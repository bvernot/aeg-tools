import pysam
from collections import Counter
import argparse


def get_base_counts(bam_file, reference_name, trim_ends, min_bqual):
    # Open the BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    # Fetch reads aligned to the given reference
    reads = bam.fetch(reference_name)
    
    # Dictionary to store base calls at each position
    position_bases = {}
    position_ref_bases = {}
    
    # Iterate over each read in the BAM file
    for read in reads:

        query_sequence = read.query_sequence
        query_qualities = read.query_qualities
        aligned_pairs = read.get_aligned_pairs(matches_only=False)

        if read.has_tag('MD'):
            ref_sequence = ''.join(x[2] if not x[2] is None else 'N' for x in read.get_aligned_pairs(matches_only=False, with_seq=True))
        else:
            ref_sequence = ''
            pass
        
        # print()
        # print(read.query_name, aligned_pairs)
        # print(read.query_name, read.get_aligned_pairs(matches_only=False, with_seq=True))
        # print(query_sequence)
        # print(ref_sequence)

        #print(read.query_name, aligned_pairs)

        # Iterate over aligned positions
        for query_pos, ref_pos in aligned_pairs:
            if ref_pos is None or query_pos is None:
                continue

            ## trim the beginning and end of each read, if necessary (this also trims clipped bases, is that what we want? I think so?)
            if query_pos < trim_ends or query_pos >= (read.query_length-trim_ends):
                # print('skipping', query_pos)
                continue

            ## skip bases that don't hit our quality threshold
            if query_qualities[query_pos] < min_bqual:
                continue
            
            base = query_sequence[query_pos]
            if read.is_reverse:
                base = base.lower()
            else:
                base = base.upper()
                pass
                
            # Store base call at this reference position
            if ref_pos not in position_bases:
                position_bases[ref_pos] = []
                pass
            
            position_bases[ref_pos].append(base)
            position_ref_bases[ref_pos] = ref_sequence[query_pos] if ref_sequence != '' else ''
            pass
        pass
    return position_bases, position_ref_bases

def report_consensus(position_bases, position_ref_bases, args):

    # Build consensus sequence from the most common base at each position
    # consensus_sequence = []
    for ref_pos in sorted(position_bases.keys()):
        base_counter = Counter(c.upper() for c in position_bases[ref_pos])
        base_counter2 = Counter(position_bases[ref_pos])
        most_common_base, most_common_base_n = base_counter.most_common(1)[0]
        depth = sum(n for _,n in base_counter.items())
        freq = most_common_base_n / depth

        ref_base = position_ref_bases[ref_pos] if position_ref_bases[ref_pos] != '' else 'no_MD_field'

        if most_common_base == 'C' and 'T' in base_counter2:
            del base_counter2['T']
        elif most_common_base == 'G' and 'a' in base_counter2:
            del base_counter2['a']
            pass

        depth2 = sum(n for _,n in base_counter2.items())
        base_counter2_upper = Counter()
        for key, count in base_counter2.items():
            base_counter2_upper[key.upper()] += count
            pass

        ####
        #somehow use args.min_depth and args.min_frequency to.. filter? create a report? flag bases?

        flags = 'PASS'
        flags2 = 'PASS'

        if depth >= args.min_depth and freq < args.min_frequency:
            flags = 'FAIL_ALL'
        elif depth < args.min_depth:
            flags = 'LOWCOV_ALL'
            pass

        if depth2 >= args.min_depth and most_common_base_n / depth2 < args.min_frequency:
            flags2 = 'FAIL_RMDMG'
        elif depth2 < args.min_depth:
            flags2 = 'LOWCOV_RMDMG'
            pass

        ## check to see if we're only reporting failed calls

        if args.report == 'failed' and 'FAIL' not in flags: continue
        if args.report == 'failed-rmdmg' and 'FAIL' not in flags2: continue

        ## print results
            
        print('CONS', ref_pos+1, ref_base,
              depth, most_common_base, most_common_base_n, freq,
              base_counter['A'], base_counter['C'], base_counter['G'], base_counter['T'],
              flags,
              depth2, most_common_base_n / depth2,
              base_counter2_upper['A'], base_counter2_upper['C'], base_counter2_upper['G'], base_counter2_upper['T'],
              flags2,
              ''.join(sorted(position_bases[ref_pos])),
              sep='\t')
        pass
    
    # Join the list of bases into a single string (the consensus sequence)
    # return ''.join(consensus_sequence)
    return


def main():
    # Set up argparse
    parser = argparse.ArgumentParser(description="Generate a consensus sequence from a BAM file.")
    parser.add_argument('-b', "--bam", help="Input BAM file")
    parser.add_argument('-c', "--chr", help="Reference sequence name (e.g., chromosome)")
    parser.add_argument('-d', "--min-depth", help="Minimum depth to report a position", type=int, default=0)
    parser.add_argument('-f', "--min-frequency", help="Minimum proportion of reads required for consensus to be called", type=float, default=.75)
    parser.add_argument('-t', "--trim-ends", help="Remove first and last X bases from each read", type=int, default=0)
    parser.add_argument('-q', "--min-bqual", help="Minimum base quality to report a base", type=int, default=25)
    parser.add_argument('--report', help="Which bases to report (default = all)", choices = ['all', 'failed', 'failed-rmdmg'], default = 'all')
    
    args = parser.parse_args()
    
    # Call the function to get the consensus sequence
    position_bases, position_ref_bases = get_base_counts(args.bam, args.chr, args.trim_ends, args.min_bqual)

    print('CONS', 'ref_pos', 'ref_base', 'depth', 'maj_base', 'n_maj_base', 'freq_maj_base',
          'n_A', 'n_C', 'n_G', 'n_T', 'flags',
          'depth_rmdmg', 'freq_maj_base_rmdmg',
          'n_A_rmdmg', 'n_C_rmdmg', 'n_G_rmdmg', 'n_T_rmdmg',
          'flags_rmdmg', 'allbases', sep='\t')
    report_consensus(position_bases, position_ref_bases, args)
    
    # Output the consensus sequence
    # print("Consensus sequence:", consensus)

if __name__ == "__main__":
    main()
