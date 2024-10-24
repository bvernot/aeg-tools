import pysam
import csv
import argparse

def load_snps(snp_file):
    """Load SNPs from a file into a dictionary with chromosome and position as keys."""
    snps = {}
    with open(snp_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        nrow = 1
        for row in reader:

            chrom, pos, ref, allele1, allele2 = row[:5]
            if nrow == 1 and not pos.isnumeric(): continue
            nrow += 1
            
            pos = int(pos)
            # Only store transition SNPs
            if is_transition(allele1, allele2):
                if chrom not in snps:
                    snps[chrom] = {}
                    pass
                snps[chrom][pos] = (allele1, allele2)
                pass
            pass
        pass
    return snps

def is_transition(base1, base2):
    """Check if the SNP is a transition (A <-> G or C <-> T)."""
    transitions = [('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')]
    return (base1, base2) in transitions

def modify_read(read, snps, chrom, X, verbose):
    """Modify the read by replacing bases overlapping with transition SNPs based on the specified conditions."""
    if chrom not in snps:
        return read, False  # No SNPs for this chromosome

    is_modified = False
    query_sequence = list(read.query_sequence)  # Convert to list to modify the sequence
    read_length = len(query_sequence)  # Total length of the read
    
    for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
        if ref_pos is None:
            continue
        
        # Check if the reference position matches a transition SNP
        if ref_pos + 1 in snps[chrom]:  # SNP positions are 1-based
            allele1, allele2 = snps[chrom][ref_pos + 1]
            read_base = query_sequence[query_pos]
            og_sequence = ''.join(query_sequence)
            
            # Check if the SNP is within the first X or last X bases of the read
            if X is None or query_pos < X or query_pos >= read_length - X:

                # Forward strand: If the read has 'T' and it's a transition SNP
                if not read.is_reverse and read_base == 'T':
                    query_sequence[query_pos] = 'N'
                    is_modified = True
                    
                # Reverse strand: If the read has 'A' and it's a transition SNP
                elif read.is_reverse and read_base == 'A':
                    query_sequence[query_pos] = 'N'
                    is_modified = True
                    pass

                pass

            if verbose: print( not read.is_reverse, allele1, allele2, og_sequence,
                               query_sequence[query_pos], query_pos, read_length - query_pos, sep='\t' )
            if verbose: print( '     ', ' ', ' ', ''.join(query_sequence), sep='\t' )
                
            pass
        pass
    
    # Convert list back to string and update the read's sequence
    read.query_sequence = ''.join(query_sequence)
    return read, is_modified

def process_bam(bam_file, snps, output_bam, X, write_mode, verbose):
    """Read BAM file, modify reads, and write to a new BAM file."""
    bam_in = pysam.AlignmentFile(bam_file, "rb")
    bam_out = pysam.AlignmentFile(output_bam, "wb", header=bam_in.header)
    
    reads_removed = 0
    n_reads = 0
    for read in bam_in.fetch():
        chrom = bam_in.get_reference_name(read.reference_id)
        modified_read, is_modified = modify_read(read, snps, chrom, X, verbose)
        n_reads += 1
        if is_modified: reads_removed += 1
        if write_mode == 'filter' and is_modified:
            if verbose: print('removing read!')
        else:
            bam_out.write(modified_read)
            pass
        pass

    print('%s %d reads out of %d' % ('Removed' if write_mode == 'filter' else 'Masked',
                                     reads_removed, n_reads))
    
    bam_in.close()
    bam_out.close()

def main():
    parser = argparse.ArgumentParser(description="Modify BAM reads by changing SNP-overlapping bases to 'N'.")
    parser.add_argument("--bam", '-b', help="Input BAM file")
    parser.add_argument("--snps", '-s', help="Input SNP file (TSV format: chrom, position, ref, allele1, allele2)")
    parser.add_argument("--output-bam", '-o', help="Output BAM file with modified reads")
    parser.add_argument("--X", '-X', type=int, default=None, help="Number of bases from the start/end of the read to check for SNP overlap. By default the entire read is processed.")
    parser.add_argument("--verbose", "-v", help="Verbose output (mostly for debugging)", action='store_true')
    parser.add_argument("--mode", '-m', choices=['mask', 'filter'], default='filter',
                        help="Mask bases vs simply removing (filter) the read from the output bam. Default is to filter the read.")
    args = parser.parse_args()
    
    # Load SNPs from file
    snps = load_snps(args.snps)
    
    # Process BAM file and write modified reads to the output BAM
    process_bam(args.bam, snps, args.output_bam, args.X, args.mode, args.verbose)

if __name__ == "__main__":
    main()
