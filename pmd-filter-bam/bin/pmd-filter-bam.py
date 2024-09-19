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

def modify_read(read, snps, chrom, X):
    """Modify the read by replacing bases overlapping with transition SNPs based on the specified conditions."""
    if chrom not in snps:
        return read  # No SNPs for this chromosome
    
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
                
                # Reverse strand: If the read has 'A' and it's a transition SNP
                elif read.is_reverse and read_base == 'A':
                    query_sequence[query_pos] = 'N'
                    pass

                pass

            print( not read.is_reverse, allele1, allele2, og_sequence, query_sequence[query_pos], query_pos, read_length - query_pos, sep='\t' )
            print( '     ', ' ', ' ', ''.join(query_sequence), sep='\t' )
                
            pass
        pass
    
    # Convert list back to string and update the read's sequence
    read.query_sequence = ''.join(query_sequence)
    return read

def process_bam(bam_file, snps, output_bam, X):
    """Read BAM file, modify reads, and write to a new BAM file."""
    bam_in = pysam.AlignmentFile(bam_file, "rb")
    bam_out = pysam.AlignmentFile(output_bam, "wb", header=bam_in.header)
    
    for read in bam_in.fetch():
        chrom = bam_in.get_reference_name(read.reference_id)
        modified_read = modify_read(read, snps, chrom, X)
        bam_out.write(modified_read)
    
    bam_in.close()
    bam_out.close()

def main():
    parser = argparse.ArgumentParser(description="Modify BAM reads by changing SNP-overlapping bases to 'N'.")
    parser.add_argument("bam_file", help="Input BAM file")
    parser.add_argument("snp_file", help="Input SNP file (TSV format: chrom, position, ref, allele1, allele2)")
    parser.add_argument("output_bam", help="Output BAM file with modified reads")
    parser.add_argument("--X", type=int, default=None, help="Number of bases from the start/end of the read to check for SNP overlap")
    args = parser.parse_args()
    
    # Load SNPs from file
    snps = load_snps(args.snp_file)
    
    # Process BAM file and write modified reads to the output BAM
    process_bam(args.bam_file, snps, args.output_bam, args.X)

if __name__ == "__main__":
    main()
