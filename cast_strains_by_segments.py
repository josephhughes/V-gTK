import csv
from collections import defaultdict
import argparse
import sys
# TODO: modify the script so it will work for different segmented viruses because they don't always have the same number of segments


# Function to pivot the data
def pivot_data(input_file, output_file, required_segments, delimiter="\t"):
    # Dictionary to hold data in a wide format, keyed by Parsed_strain
    # Each strain will have a dict where the key is the segment and value is a list of Locus values
    data_wide = defaultdict(lambda: defaultdict(list))
    
    # Set to collect all segments (which will become column names)
    segments = set()

    # Counters for complete and incomplete genomes
    complete_count = 0
    incomplete_count = 0

    try:
        # Read the input CSV file
        with open(input_file, mode='r', newline='') as csv_file:
            reader = csv.DictReader(csv_file, delimiter=delimiter)
            #print("Input file headers:", reader.fieldnames)
            
            # Check for required columns
            if not {'Parsed_strain', 'segment', 'Locus'}.issubset(reader.fieldnames):
                raise ValueError("Input file is missing one or more required columns: 'Parsed_strain', 'segment', 'Locus'")

            # Iterate over each row in the CSV
            for row in reader:
                # Get the key for grouping (Parsed_strain)
                parsed_strain = row['Parsed_strain']
                # Get the segment and locus values
                segment = row['segment']
                locus = row['Locus']

                # Append the locus value to the list for the appropriate segment
                data_wide[parsed_strain][segment].append(locus)
                
                # Track the segment for the column headers later
                segments.add(segment)
    except FileNotFoundError:
        sys.exit(f"Error: The file {input_file} was not found.")
    except ValueError as ve:
        sys.exit(f"Error: {ve}")
    except Exception as e:
        sys.exit(f"An unexpected error occurred: {e}")

    # Convert the set of segments to a sorted list for consistent column order
    segments = sorted(segments)

    # Write the pivoted data into a new CSV
    with open(output_file, mode='w', newline='') as csv_output_file:
        # Define the column names for the output file
        fieldnames = ['Parsed_strain'] + segments + ['Complete_status'] # 'Parsed_strain' + all the segments as column names + whether the genome is complete with all 8 segments

        # Create a DictWriter instance
        writer = csv.DictWriter(csv_output_file, fieldnames=fieldnames, delimiter=delimiter)

        # Write the header row
        writer.writeheader()

        # Write the rows of the wide data
        for strain, loci_dict in data_wide.items():
            # Create a row dictionary
            row = {'Parsed_strain': strain}
            
            # Populate the row with comma-separated Locus values for each segment
            for segment in segments:
                # Join multiple Locus values with a comma if they exist
                row[segment] = ','.join(loci_dict.get(segment, []))  # Leave empty if there's no locus for this segment

            # Check if all required segments have at least one Locus
            complete = all(bool(loci_dict.get(segment)) for segment in required_segments)

            # Add the "Complete" or "Incomplete" status
            row['Complete_status'] = 'Complete' if complete else 'Incomplete'

              # Update counters for complete and incomplete genomes
            if complete:
                complete_count += 1
            else:
                incomplete_count += 1
            
               
            # Write the row to the CSV
            writer.writerow(row)

    # Print the counts of Complete and Incomplete genomes
    print(f"Number of Complete genomes: {complete_count}")
    print(f"Number of Incomplete genomes: {incomplete_count}")

# Set up argument parser
parser = argparse.ArgumentParser(description='Compile the per strain accession numbers to validate a complete genome.')
parser.add_argument('input_file', help='Path to the input TSV file')
parser.add_argument('output_file', help='Path to the output TSV file')
parser.add_argument('--delimiter', default="\t", help='Delimiter used in the input and output files (default: tab)')

# Parse the command-line arguments
args = parser.parse_args()

# List of required segments (in your case, 8 segments)
required_segments = ['1', '2', '3', '4', '5', '6', '7', '8']

# Call the function to pivot the data and write to the output CSV
pivot_data(args.input_file, args.output_file, required_segments, args.delimiter)

print(f"Pivoted data has been written to {args.output_file}")
