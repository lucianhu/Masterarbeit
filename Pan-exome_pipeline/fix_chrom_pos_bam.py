import argparse
import pysam
import os

def main():
    parser = argparse.ArgumentParser(description='Fix CHROM and POS columns of a BAM file to resemble a normal SAM/BAM.')

    parser.add_argument('position_report', help='FASTA position report tsv file')
    parser.add_argument('input_bam', help='Input BAM file')
    parser.add_argument('output_bam', help='Output edited BAM file')

    args = parser.parse_args()

    position_report = []

    with open(args.position_report, 'r') as report_file:
        for line in report_file:
            fields = line.strip().split()
            if len(fields) == 5:
                chrom, start, end, value1, value2 = fields
                position_report.append((chrom, int(start), int(end), int(value1), int(value2)))

    input_bam = pysam.AlignmentFile(args.input_bam, 'rb')
    output_bam = pysam.AlignmentFile(args.output_bam, 'wb', header=input_bam.header)

    for read in input_bam:
        chrom, pos = read.reference_name, read.reference_start
        min_distance = float('inf')
        selected_chrom = None

        for chrom_report, start_report, end_report, value1_report, value2_report in position_report:
            if pos < value2_report:
                distance = abs(pos - value2_report)
                if distance < min_distance:
                    min_distance = distance
                    ordinal_position = value1_report - min_distance
                    selected_chrom = chrom_report
                    selected_start_report = start_report

        if selected_chrom:
            updated_pos = selected_start_report + (ordinal_position - 1)
            read.reference_name = selected_chrom
            read.reference_start = updated_pos
        output_bam.write(read)

    input_bam.close()
    output_bam.close()

    # Sort and index the output BAM
    sorted_output_bam = args.output_bam.replace(".bam", ".sorted.bam")
    pysam.sort("-o", sorted_output_bam, args.output_bam)
    pysam.index(sorted_output_bam)

    # Clean up the unsorted output BAM
    os.remove(args.output_bam)

if __name__ == '__main__':
    main()


