#!/usr/bin/env python
import argparse
import csv
import logging
import sys
from enum import Enum
from typing import List, NoReturn, Optional


class ColumnNames(str, Enum):
    SAMPLE = "sample"
    FASTQ_1 = "fastq_1"
    FASTQ_2 = "fastq_2"
    FASTA = "fasta"
    SINGLE_END = "single_end"


def parse_args(args: Optional[List[str]] = None) -> argparse.Namespace:
    """
    Reformatting is based on detecting whether the reads are paired or single end.
    Script appends appropriate column to samplesheet.csv file.
    """
    parser = argparse.ArgumentParser(
        description="Reformat nf-core/taxprofiler samplesheet file.",
        epilog="Example usage: python detect_reads.py <FILE_IN> <FILE_OUT>",
    )
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
        new_file_rows = []

        with open(input_file_path, "r", newline="") as input_file:
            csv_reader = csv.DictReader(input_file, delimiter=",")
            self.headers = csv_reader.fieldnames
            self.headers.append("single_end")

            for samplesheet_row in csv_reader:

                if self._is_paired_end_short_read(samplesheet_row):
                    samplesheet_row[ColumnNames.SINGLE_END] = "0"
                    new_file_rows.append(samplesheet_row.values())

                elif self._is_single_end_short_long_read(samplesheet_row):
                    samplesheet_row[ColumnNames.SINGLE_END] = "1"
                    new_file_rows.append(samplesheet_row.values())

                elif self._is_single_end_long_read(samplesheet_row):
                    samplesheet_row[ColumnNames.SINGLE_END] = "1"
                    new_file_rows.append(samplesheet_row.values())

                elif self._is_error_row(samplesheet_row):
                    logging.error(
                        "FastQ and FastA files cannot be specified together in the same library!",
                        "Line",
                        ",".join(samplesheet_row.values()),
                    )
                else:
                    logging.error(
                        "Invalid combination of columns provided!", "Line", ",".join(samplesheet_row.values())
                    )

        ReadsModifier.save_reformatted_samplesheet([self.headers] + new_file_rows, output_file_path)

    def _get_row_values(self, samplesheet_row: dict):
        """
        This method extracts data from the columns for given row of samplesheet table.
        """
        return (
            samplesheet_row.get(ColumnNames.SAMPLE),
            samplesheet_row.get(ColumnNames.FASTQ_1),
            samplesheet_row.get(ColumnNames.FASTQ_2),
            samplesheet_row.get(ColumnNames.FASTA),
        )

    def _is_paired_end_short_read(self, samplesheet_row: dict) -> bool:
        sample, fastq_1, fastq_2, _ = self._get_row_values(samplesheet_row)
        return sample and fastq_1 and fastq_2

    def _is_single_end_short_long_read(self, samplesheet_row: dict) -> bool:
        sample, fastq_1, fastq_2, _ = self._get_row_values(samplesheet_row)
        return sample and fastq_1 and not fastq_2

    def _is_single_end_long_read(self, samplesheet_row: dict) -> bool:
        sample, fastq_1, fastq_2, fasta = self._get_row_values(samplesheet_row)
        return sample and fasta and not fastq_1 and not fastq_2

    def _is_error_row(self, samplesheet_row: dict) -> bool:
        sample, fastq_1, fastq_2, fasta = self._get_row_values(samplesheet_row)
        return fasta and (fastq_1 or fastq_2)

    @classmethod
    def save_reformatted_samplesheet(cls, new_file_rows: List[List], output_file_path: str) -> NoReturn:
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
