#!/usr/bin/env python3

import os
import sys
import time
import resource
from output.sam_output import write_sam_file
from preprocessing.preprocessor import Preprocessor
from validation.validator import Validator
from read_mapper.read_mapper import ReadMapper
import concurrent.futures


def get_memory_usage():
    """
    Get memory usage of the current process.

    Returns:
        - Memory usage in GB.
    """

    # Get memory usage in bytes and convert to GB
    mem_usage_bytes = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    mem_usage_gb = mem_usage_bytes / (1024 * 1024 * 1024)
    return mem_usage_gb


def find_unique_file(all_files, identifier, folder_path):
    """
    Find a unique file matching the identifier in the folder.

    Inputs:
        - all_files: List of all files in the folder.
        - identifier: Identifier string to match in the filenames.
        - folder_path: Path to the folder.

    Returns:
        - The path of the unique matching file.
    """

    # Filter files based on identifier
    matched_files = [f for f in all_files if identifier in f]

    # Check if there is a unique file matching the identifier
    if len(matched_files) == 1:
        return os.path.join(folder_path, matched_files[0])
    # If no files match the identifier, raise an exception
    elif not matched_files:
        raise Exception(
            f"No files with identifier '{identifier}' found in {folder_path}."
        )
    # If multiple files match the identifier, raise an exception
    else:
        matched_list = "\n".join([f"\t{x}" for x in matched_files])  # List all matches
        raise Exception(
            f"Multiple files with identifier '{identifier}' found:\n{matched_list}"
        )


def get_required_file_paths(folder_path, required_identifiers):
    """
    Retrieve required file paths based on identifiers.

    Inputs:
        - folder_path: Path to the folder containing the files.
        - required_identifiers: List of required file identifiers.

    Returns:
        - A dictionary mapping identifiers to file paths.
    """

    # Get all files in the folder
    all_files = os.listdir(folder_path)
    file_paths = {}

    # Iterate over each required identifier
    for identifier in required_identifiers:
        # Find the unique file matching the identifier
        file_paths[identifier] = find_unique_file(all_files, identifier, folder_path)
    return file_paths


def main():
    """
    Main execution function for read mapping and validation.

    Steps:
        - Parse input arguments for folder path.
        - Preprocess data and perform read mapping.
        - Validate read mappings and output metrics.
        - Write SAM file output.
    """

    # Check if folder path is provided as command-line argument
    if len(sys.argv) < 2:
        raise Exception("Folder path not provided")
    # Check if too many arguments are provided
    elif len(sys.argv) > 2:
        raise Exception("Too many arguments provided")

    # Get folder path from command-line argument
    folder_path = sys.argv[1]

    # Check if folder path is valid
    if not os.path.isdir(folder_path):
        raise Exception(f"Invalid folder path: {folder_path}")

    # List of required identifiers for input files
    required_identifiers = (
        "short_reads_ref",
        "short_reads_1",
        "short_reads_2",
        "ground_truth",
    )

    # Get file paths for required input files
    file_paths = get_required_file_paths(folder_path, required_identifiers)

    # Start timing and memory tracking for indexing process
    start_time_wall_indexing = time.time()
    start_time_cpu_indexing = time.process_time()
    mem_before_indexing = get_memory_usage()

    print("Preprocessing...")
    # Initialize Preprocessor and perform indexing
    io = Preprocessor(
        file_paths["short_reads_ref"],
        file_paths["short_reads_1"],
        file_paths["short_reads_2"],
    )

    # Measure elapsed time and memory used for indexing
    elapsed_time_wall_indexing = time.time() - start_time_wall_indexing
    elapsed_time_cpu_indexing = time.process_time() - start_time_cpu_indexing
    indexing_memory_used = get_memory_usage() - mem_before_indexing

    # Initialize ReadMapper
    read_mapper = ReadMapper(io.reference_genome, io.suffix_array)
    # Initialize Validator
    validator = Validator(file_paths["ground_truth"], io.reference_genome)
    # Get reads and quality scores
    reads_1, quals_1, reads_2, quals_2 = io.get_reads()

    # Start timing for mapping process
    start_time_wall_mapping = time.time()
    start_time_cpu_mapping = time.process_time()

    print("Mapping reads concurrently...")
    with concurrent.futures.ThreadPoolExecutor() as executor:
        future_read_1 = executor.submit(read_mapper.map_reads, reads_1, quals_1)
        future_read_2 = executor.submit(read_mapper.map_reads, reads_2, quals_2)
        read_1_alignments, read_1_memory_usage = future_read_1.result()
        read_2_alignments, read_2_memory_usage = future_read_2.result()

    # Measure elapsed time and memory used for mapping
    elapsed_time_wall_mapping = time.time() - start_time_wall_mapping
    elapsed_time_cpu_mapping = time.process_time() - start_time_cpu_mapping
    mapping_memory_used = get_memory_usage() - max(
        read_1_memory_usage, read_2_memory_usage
    )

    # Calculate performance metrics
    time_wall_mapping_min = elapsed_time_wall_mapping / 60
    total_reads = len(reads_1) + len(reads_2)
    reads_per_minute = (
        total_reads / time_wall_mapping_min if time_wall_mapping_min > 0 else 0
    )

    # Validate alignments for both sets of reads, and combine metrics
    read_1_metrics = validator.validate(read_1_alignments)
    read_2_metrics = validator.validate(read_2_alignments)
    combined_metrics = validator.combine_metrics([read_1_metrics, read_2_metrics])

    # Print metrics and performance information
    validator.print_metrics(
        combined_metrics,
        total_reads,
        elapsed_time_wall_indexing,
        elapsed_time_cpu_indexing,
        indexing_memory_used,
        elapsed_time_wall_mapping,
        elapsed_time_cpu_mapping,
        mapping_memory_used,
        reads_per_minute,
    )

    # Write SAM file output
    sam_filename = "output.sam"
    write_sam_file(sam_filename, read_1_alignments, read_2_alignments, read_mapper, reads_1, reads_2, quals_1, quals_2, io.reference_genome_name)

    print(f"SAM file written to {sam_filename}")


# Usage: ./main.py <folder_path>
if __name__ == "__main__":
    main()
