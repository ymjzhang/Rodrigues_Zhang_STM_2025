import os
import csv
import argparse
# Created on Jun 14 2024
# Created by Jason Zhang
# Last updated by Jason Zhang
# Last updated on June 14 2024

# Usage 
# python3 batch_rename.py /path/to/prefixes.csv
# if specify a different directory 
# python3 batch_rename.py --directory /path/to/your/directory /path/to/prefixes.csv


def rename_files_in_directory(directory, csv_file):
    # Read the CSV file and store the prefix mappings in a dictionary
    prefix_mapping = {}
    with open(csv_file, mode='r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip the header row if it exists
        for row in reader:
            if len(row) < 2:
                continue  # Skip rows that do not have enough columns
            current_prefix, new_prefix = row
            prefix_mapping[current_prefix] = new_prefix

    # List all files in the directory
    for filename in os.listdir(directory):
        # Check if the file name starts with any of the current prefixes
        for current_prefix, new_prefix in prefix_mapping.items():
            if filename.startswith(current_prefix):
                # Construct the new filename
                new_filename = filename.replace(current_prefix, new_prefix, 1)
                old_file_path = os.path.join(directory, filename)
                new_file_path = os.path.join(directory, new_filename)
                os.rename(old_file_path, new_file_path)
                print(f'Renamed: {filename} -> {new_filename}')
                break  # No need to check other prefixes once a match is found

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Batch rename files in a directory based on a CSV mapping file.')
    parser.add_argument('--directory', type=str, default=os.getcwd(),
                        help='The directory containing the files to be renamed (default: current directory).')
    parser.add_argument('csv_file', type=str, help='The CSV file containing the prefix mappings.')

    args = parser.parse_args()
    rename_files_in_directory(args.directory, args.csv_file)

