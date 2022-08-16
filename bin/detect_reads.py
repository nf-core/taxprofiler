#!/usr/bin/env python

import argparse
import csv
import sys
from typing import List, NoReturn


def parse_args(args=None) -> argparse.Namespace:
    """
    Reformatting is based on detecting whether the reads are paired or single end.
    Script appends appropriate column to samplesheet.csv file.
    """
    Description = "Reformat nf-core/taxprofiler samplesheet file."
    Epilog = "Example usage: python detect_reads.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


class ReadsModifier:
    def __init__(self):
        self.headers = None
        self.sample_index = None
        self.fastq_1_index = None
        self.fastq_2_index = None
        self.fasta_index = None

    def detect_reads_and_reformat(self, input_file_path: str, output_file_path: str) -> NoReturn:
        NEW_COLUMN_NAME = "single_end"
        new_file_rows = []

        with open(input_file_path, "r") as input_file:
            csv_reader = csv.reader(input_file, delimiter=",")
            self.headers = next(csv_reader)
            self.headers.append(NEW_COLUMN_NAME)

            self._infer_column_indexes()

            for samplesheet_row in csv_reader:

                if self._is_paired_end_short_read(samplesheet_row):
                    new_file_rows.append([*samplesheet_row, "0"])

                elif self._is_single_end_short_long_read(samplesheet_row):
                    new_file_rows.append([*samplesheet_row, "1"])

                elif self._is_single_end_long_read(samplesheet_row):
                    new_file_rows.append([*samplesheet_row, "1"])

                elif self._is_error_row(samplesheet_row):
                    self.print_error(
                        "FastQ and FastA files cannot be specified together in the same library!",
                        "Line",
                        ",".join(samplesheet_row),
                    )
                else:
                    self.print_error("Invalid combination of columns provided!", "Line", ",".join(samplesheet_row))

        self.save_reformatted_samplesheet([self.headers] + new_file_rows, output_file_path)

    def _get_row_values(self, samplesheet_row):
        """
        This method extracts data from the columns for given row of samplesheet table, based on
        previously infered column indexes.
        """
        sample = samplesheet_row[self.sample_index]
        fastq_1 = samplesheet_row[self.fastq_1_index] if self.fastq_1_index else None
        fastq_2 = samplesheet_row[self.fastq_2_index] if self.fastq_2_index else None
        fasta = samplesheet_row[self.fasta_index] if self.fasta_index else None
        return sample, fastq_1, fastq_2, fasta

    def _infer_column_indexes(self):
        """
        This method infers indexes of necessary columns from samplesheet table
        """
        self.sample_index = self.headers.index("sample")
        self.fastq_1_index = self.headers.index("fastq_1") if "fastq_1" in self.headers else None
        self.fastq_2_index = self.headers.index("fastq_2") if "fastq_2" in self.headers else None
        self.fasta_index = self.headers.index("fasta") if "fasta" in self.headers else None

    def _is_paired_end_short_read(self, samplesheet_row: List) -> bool:
        sample, fastq_1, fastq_2, _ = self._get_row_values(samplesheet_row)
        return sample and fastq_1 and fastq_2

    def _is_single_end_short_long_read(self, samplesheet_row: List) -> bool:
        sample, fastq_1, fastq_2, _ = self._get_row_values(samplesheet_row)
        return sample and fastq_1 and not fastq_2

    def _is_single_end_long_read(self, samplesheet_row: List) -> bool:
        sample, fastq_1, fastq_2, fasta = self._get_row_values(samplesheet_row)
        return sample and fasta and not fastq_1 and not fastq_2

    def _is_error_row(self, samplesheet_row: List) -> bool:
        sample, fastq_1, fastq_2, fasta = self._get_row_values(samplesheet_row)
        return fasta and (fastq_1 or fastq_2)

    @staticmethod
    def print_error(error: str, context: str = "Line", context_str: str = ""):
        error_str = "ERROR: Please check samplesheet -> {}".format(error)
        if context != "" and context_str != "":
            error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
                error, context.strip(), context_str.strip()
            )
        print(error_str)
        sys.exit(1)

    @staticmethod
    def save_reformatted_samplesheet(new_file_rows: List[List], output_file_path: str) -> NoReturn:
        """
        Write new samplesheet.
        """
        with open(output_file_path, "w") as output_file:
            csv.writer(output_file).writerows(new_file_rows)


def main(args=None):
    args = parse_args(args)
    ReadsModifier().detect_reads_and_reformat(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
