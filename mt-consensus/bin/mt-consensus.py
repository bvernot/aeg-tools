import pysam
from collections import Counter
import argparse


def get_base_counts(bam_file, reference_name, trim_ends = 0):
    # Open the BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    # Fetch reads aligned to the given reference
    reads = bam.fetch(reference_name)
    
    # Dictionary to store base calls at each position
    position_bases = {}
    
    # Iterate over each read in the BAM file
    for read in reads:
        query_sequence = read.query_sequence
        query_qualities = read.query_qualities
        aligned_pairs = read.get_aligned_pairs(matches_only=False)

        #print(read.query_name, aligned_pairs)
        #print(read.query_name, read.get_aligned_pairs(matches_only=False))

        #print(read.query_name, aligned_pairs)

        # Iterate over aligned positions
        for query_pos, ref_pos in aligned_pairs:
            if ref_pos is None or query_pos is None:
                continue

            ## trim the beginning and end of each read, if necessary (this also trims clipped bases, is that what we want? I think so?)
            if query_pos < trim_ends or query_pos >= (read.query_length-trim_ends):
                # print('skipping', query_pos)
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
            pass
        pass
    return position_bases

def report_consensus(position_bases, args):

    # Build consensus sequence from the most common base at each position
    consensus_sequence = []
    for ref_pos in sorted(position_bases.keys()):
        base_counter = Counter(c.upper() for c in position_bases[ref_pos])
        base_counter2 = Counter(position_bases[ref_pos])
        most_common_base, most_common_base_n = base_counter.most_common(1)[0]
        depth = sum(n for _,n in base_counter.items())
        freq = most_common_base_n / depth

        if most_common_base == 'C' and 'T' in base_counter2:
            del base_counter2['T']
        elif most_common_base == 'G' and 'a' in base_counter2:
            del base_counter2['a']
            pass
        depth2 = sum(n for _,n in base_counter2.items())
        
        print('CONS', ref_pos+1, depth, depth2, most_common_base, most_common_base_n, freq, most_common_base_n / depth2, base_counter, base_counter2, Counter(position_bases[ref_pos]), sep='\t')
        pass
    
    # Join the list of bases into a single string (the consensus sequence)
    # return ''.join(consensus_sequence)
    return


def main():
    # Set up argparse
    parser = argparse.ArgumentParser(description="Generate a consensus sequence from a BAM file.")
    parser.add_argument('-b', "--bam", help="Input BAM file")
    parser.add_argument('-c', "--chr", help="Reference sequence name (e.g., chromosome)")
    parser.add_argument('-d', "--depth", help="Minimum depth to report a position", type=int, default=5)
    parser.add_argument('-f', "--frequency", help="Minimum proportion of reads required for consensus to be called", type=float, default=.75)
    parser.add_argument('-t', "--trim-ends", help="Remove first and last X bases from each read", type=int, default=0)
    
    args = parser.parse_args()
    
    # Call the function to get the consensus sequence
    position_bases = get_base_counts(args.bam, args.chr, args.trim_ends)

    print('CONS', 'ref_pos', 'depth', 'depth2', 'most_common_base', 'most_common_base_n', 'freq', 'most_common_base_n / depth2', 'base_counter', 'base_counter2', 'Counter(position_bases[ref_pos])', sep='\t')
    report_consensus(position_bases, args)
    
    # Output the consensus sequence
    # print("Consensus sequence:", consensus)

if __name__ == "__main__":
    main()
