from __future__ import annotations

import argparse

from filter_fastq import filter_fastq


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Filter FASTQ reads by GC content, length, and mean quality."
    )

    parser.add_argument("input_fastq", help="Path to input FASTQ file")
    parser.add_argument("output_fastq", help="Path to output FASTQ file")

    parser.add_argument(
        "--gc-bounds",
        nargs=2,
        type=float,
        metavar=("MIN_GC", "MAX_GC"),
        default=(0.0, 100.0),
        help="Inclusive GC percentage bounds",
    )

    parser.add_argument(
        "--length-bounds",
        nargs=2,
        type=int,
        metavar=("MIN_LEN", "MAX_LEN"),
        default=(0, 2**32),
        help="Inclusive read length bounds",
    )

    parser.add_argument(
        "--quality-threshold",
        type=float,
        default=0.0,
        help="Minimum mean Phred quality",
    )

    parser.add_argument(
        "--log-file",
        type=str,
        default="fastq_filter.log",
        help="Path to log file",
    )

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    written = filter_fastq(
        input_fastq=args.input_fastq,
        output_fastq=args.output_fastq,
        gc_bounds=tuple(args.gc_bounds),
        length_bounds=tuple(args.length_bounds),
        quality_threshold=args.quality_threshold,
        log_file=args.log_file,
    )

    print(f"Reads written: {written}")


if __name__ == "__main__":
    main()