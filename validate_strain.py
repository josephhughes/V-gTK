import csv
import re
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description='Extract influenza A strain names from the Definition field of a TSV file.')
parser.add_argument('input_file', help='Path to the input TSV file')
parser.add_argument('output_file', help='Path to the output TSV file')
args = parser.parse_args()

# Regular expression to extract the influenza A strain name
pattern = r"A\/[A-Za-z]+\/[A-Za-z0-9\-]+\/\d+"

# Open and read the TSV file
with open(args.input_file, 'r') as infile, open(args.output_file, 'w') as outfile:
    # Set up DictReader to read the input TSV file
    reader = csv.DictReader(infile, delimiter='\t')
    
    # Create fieldnames for the output TSV, adding a new 'parsed' column
    fieldnames = reader.fieldnames + ['parsed']
    
    # Set up DictWriter to write to the output TSV file
    writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
    
    # Write header to the output file
    writer.writeheader()

    # Process each row
    for row in reader:
        # Extract 'definition' column
        definition = row['Definition']
        # Ensure 'definition' is a string
        if definition is None:
            print(row['Locus'])
            parsed_text=''
        else:
            # Use regex to extract the required pattern from the 'definition' column
            match = re.search(pattern, definition)
            parsed_text = match.group(0) if match else ''  # If there's no match, return an empty string

        # Add the parsed text as a new column to the row
        row['parsed'] = parsed_text

        # Write the updated row to the output TSV file
        writer.writerow(row)