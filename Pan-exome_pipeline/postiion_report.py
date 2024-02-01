#!/usr/bin/env python3
# Position Report Generator for BED file.

import argparse

# Define a function to process a single BED file and generate TSV data
def process_bed(input_bed, output_tsv):
    # Initialize an empty list to store the results
    results = []

    # Initialize a variable to accumulate lengths
    accumulated_ordinal_length = 0

    # Open the input BED file
    with open(input_bed, 'r') as bed_file:
        # Process each line in the BED file
        for line in bed_file:
            line = line.strip()  # Remove leading/trailing whitespace
            if line.startswith('#') or not line:
                continue  # Skip comment lines and empty lines
            
            # Split the line into columns
            columns = line.split('\t')

            # Extract chromosome, start, and end positions
            chromosome = columns[0]
            start = columns[1]
            end = columns[2]

            ordinal_length = int(end) - int(start) + 1
            accumulated_ordinal_length += ordinal_length

            # Use tab characters ('\t') to separate fields in the TSV file
            info = f"{chromosome}\t{start}\t{end}\t{ordinal_length}\t{accumulated_ordinal_length}"
            results.append(info)

    # Write the results to the output TSV file
    with open(output_tsv, 'w') as tsv_file:
        tsv_data = "\n".join(results)
        tsv_file.write(tsv_data)

    print(f"Processed {len(results)} records and saved results to {output_tsv}")

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description='Position Report Generator for BED file.')

    # Add arguments for input BED file and output TSV file
    parser.add_argument('input_bed', help='Input BED file')
    parser.add_argument('output_tsv', help='Output TSV file')

    # Parse command line arguments
    args = parser.parse_args()

    # Process the input BED file and generate the output TSV file
    process_bed(args.input_bed, args.output_tsv)

if __name__ == '__main__':
    main()
