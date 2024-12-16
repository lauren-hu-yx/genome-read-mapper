class Validator:
    def __init__(self, file_path, ref):
        """
        Initialize Validator with ground truth data and reference genome.

        Inputs:
            - file_path: Path to the ground truth file.
            - ref: The reference genome sequence.
        """

        # store reference genome and ground truth data
        self.reference = ref
        self.ground_truth = self.load_ground_truth(file_path)

    def load_ground_truth(self, file_path):
        """
        Load ground truth data from a file.

        Inputs:
            - file_path: Path to the ground truth file.

        Returns:
            - ground_truth: A dictionary mapping read IDs to (start, end) positions.
        """

        ground_truth = {}

        # for each line in the file, split the line into parts
        with open(file_path, "r") as f:
            for line_num, line in enumerate(f, 1):
                parts = line.strip().split("\t")
                if len(parts) != 3:
                    print(
                        f"Warning: Line {line_num} in {file_path} does not have exactly 3 fields. Skipping."
                    )
                    continue
                # extract read ID, start, and end positions
                read_id, start, end = parts
                try:
                    # store start and end positions as integers
                    ground_truth[read_id] = (
                        int(start),
                        int(end),
                    )
                except ValueError:
                    print(
                        f"Warning: Start or end position on line {line_num} is not an integer. Skipping."
                    )
        return ground_truth

    def validate(self, predictions):
        """
        Validate predictions against the ground truth data.

        Inputs:
            - predictions: A dictionary mapping read IDs to predicted positions.

        Returns:
            - metrics: A dictionary of validation metrics (tp, fp, tn, fn, ia).
        """

        metrics = {"tp": 0, "fp": 0, "tn": 0, "fn": 0, "ia": 0}

        # iterate over each read ID and prediction
        for read_id, prediction in predictions.items():
            # check if read ID is in ground truth data
            if read_id not in self.ground_truth:
                # if not in ground truth, and prediction is none, then it is a true negative
                if prediction is None:
                    metrics["tn"] += 1  # Correctly unaligned
                # if not in ground truth, and prediction is not none, then it is a false positive
                else:
                    metrics["fp"] += 1  # Incorrectly aligned
                continue

            # extract start position from ground truth data
            truth_start, truth_end = self.ground_truth[read_id]

            # if read ID is in ground truth and prediction is none, then it is a false negative
            if prediction is None:
                metrics["fn"] += 1  # Incorrectly unaligned
            # if read ID is in ground truth and prediction is not none, then check alignment
            else:
                pred_start, pred_end, _ = prediction
                if pred_start == truth_start and pred_end == truth_end:
                    metrics["tp"] += 1  # Correctly aligned
                else:
                    # metrics["ia"] += 1  # Alignment is off, we use this metric for testing purposes
                    # We can classify incorrect alignment as false positive
                    metrics["fp"] += 1
        return metrics

    def combine_metrics(self, metrics_list):
        """
        Combine multiple sets of validation metrics.

        Inputs:
            - metrics_list: A list of metrics dictionaries to be combined.

        Returns:
            - combined: A dictionary containing the combined metrics.
        """

        combined = {
            "tp": 0,
            "fp": 0,
            "tn": 0,
            "fn": 0,
            "ia": 0,
        }

        # iterate over each set of metrics and sum the values
        for metrics in metrics_list:
            for key in combined:
                combined[key] += metrics.get(key, 0)
        return combined

    def print_metrics(
        self,
        metrics,
        total_reads,
        elapsed_time_wall_indexing,
        elapsed_time_cpu_indexing,
        indexing_memory_used,
        elapsed_time_wall_mapping,
        elaspsed_time_cpu_mapping,
        mapping_memory_used,
        reads_per_minute,
    ):
        """
        Print validation metrics and additional indexing/mapping information.

        Inputs:
            - metrics: The validation metrics dictionary.
            - total_reads: The total number of reads.
            - elapsed_time_wall_indexing: Wall clock time for indexing.
            - elapsed_time_cpu_indexing: CPU time for indexing.
            - indexing_memory_used: Memory used during indexing (in GB).
            - elapsed_time_wall_mapping: Wall clock time for mapping.
            - elaspsed_time_cpu_mapping: CPU time for mapping.
            - mapping_memory_used: Memory used during mapping (in GB).
            - reads_per_minute: Reads processed per minute.
        """

        # extract metrics values
        tp, fp, tn, fn, ia = (
            metrics["tp"],
            metrics["fp"],
            metrics["tn"],
            metrics["fn"],
            metrics["ia"],
        )

        # calculate precision and recall
        precision = (100 * tp / (tp + fp)) if (tp + fp) > 0 else 0.0
        recall = (100 * tp / (tp + fn)) if (tp + fn) > 0 else 0.0

        # format and print metrics output
        metrics_output = (
            "\nValidation Metrics:\n"
            "--------------------------------\n"
            f"Total Reads      : {total_reads}\n\n"
            "- Indexing -\n"
            f"Wall Clock Time  : {elapsed_time_wall_indexing:.2f} seconds\n"
            f"CPU Time         : {elapsed_time_cpu_indexing:.2f} seconds\n"
            f"Memory Used      : {indexing_memory_used:.2f} GB\n\n"
            "- Mapping -\n"
            f"Wall Clock Time  : {elapsed_time_wall_mapping:.2f} seconds\n"
            f"CPU Time         : {elaspsed_time_cpu_mapping:.2f} seconds\n"
            f"Memory Used      : {mapping_memory_used:.2f} GB\n"
            f"Reads per Minute : {reads_per_minute:.0f}\n\n"
            f"True Positive    : {tp}\n"
            f"False Positive   : {fp}\n"
            f"True Negative    : {tn}\n"
            f"False Negative   : {fn}\n\n"
            # f"Flawed Alignment : {ia}\n\n"
            f"Precision        : {precision:.2f}%\n"
            f"Recall           : {recall:.2f}%\n"
            "--------------------------------\n"
        )

        print(metrics_output)
