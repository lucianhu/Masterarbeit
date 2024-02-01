import argparse
import gzip
import os

def main():
    parser = argparse.ArgumentParser(description='Fix CHROM and POS columns of a pangenome VCF file to resemble a normal VCF.')

    parser.add_argument('position_report', help='FASTA position report tsv file')
    parser.add_argument('input_vcf', help='Input pangenome VCF file')
    parser.add_argument('output_vcf', help='Output edited pangenome VCF file')

    args = parser.parse_args()

    # Determine whether the input file is gzip-compressed or not
    if args.input_vcf.endswith('.gz'):
        vcf_file = gzip.open(args.input_vcf, 'rt')
    else:
        vcf_file = open(args.input_vcf, 'r')

    position_report = []

    with open(args.position_report, 'r') as report_file:
        for line in report_file:
            fields = line.strip().split()
            if len(fields) == 5:
                chrom, start, end, value1, value2 = fields
                position_report.append((chrom, int(start), int(end), int(value1), int(value2)))

    updated_vcf_lines = []

    for line in vcf_file:
        if line.startswith('#'):
            updated_vcf_lines.append(line)
        else:
            fields = line.strip().split('\t')
            if len(fields) >= 5:
                chrom, pos, _, _, _ = fields[0:5]
                pos = int(pos)
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
                    fields[0] = selected_chrom
                    fields[1] = str(updated_pos)
                    updated_vcf_lines.append('\t'.join(fields) + '\n')

    vcf_file.close()

    # Open the output VCF file in plain text mode
    with open(args.output_vcf, 'w') as output_vcf_file:
        output_vcf_file.writelines(updated_vcf_lines)

if __name__ == '__main__':
    main()

