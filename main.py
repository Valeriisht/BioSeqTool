import argparse
from loguru import logger
import sys
from pathlib import Path

from NASeqTool import filter_fastq

# Вызов данной функции приведёт к изменению свойств root logger
logs_format = "<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | <level>{level: <8}</level> | <cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>"
logger.add(  # Добавляем Handler для записи в stdout
    sys.stdout, level="DEBUG", format=logs_format, colorize=True
)


logger.add(
    "errors_logs.log",
    level="ERROR",
    format=logs_format,
    colorize=False,
    rotation="500 KB",  # ограничиваем по памяти
)


def parser_tuple(val):
    """Parsing the gc_bounds/lengh_bound argument into tuple or float."""
    try:
        if len(val) > 1:
            min, max = float(val[0]), float(val[1])
            return min, max
        return 0, val
    except (TypeError, ValueError) as e:
        raise ValueError(f"Invalid value with {e}")


def parse_args():
    """Function for parsing commands arguments"""

    parser = argparse.ArgumentParser(
        prog="filter_fastq",
        description="The program provides filtering of input files by quality, GC-composition and sequence length.",
        epilog="Example: python script.py -I input.fq -O output.fq -gc 20 80 -l 50 150 -q 30",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "-I", "--input", type=Path, required=True, help="Input FASTQ file path"
    )  # positional argume

    parser.add_argument(
        "-O",
        "--output",
        type=Path,
        default=Path("output_filtered.fastq"),
        help="Output FASTQ file path  (default: <input>_filtered.fq)",
    )

    parser.add_argument(
        "-gc",
        type=int,
        nargs=2,
        metavar=("MIN", "MAX"),
        default=[0, 100],
        help="param for filtering by gc-composition",
    )
    parser.add_argument(
        "-l",
        "--length",
        type=int,
        nargs=2,
        metavar=("MIN", "MAX"),
        default=[0, 2**32],
        help="param for filtering by length_bounds",
    )
    parser.add_argument(
        "-q", "--quality", type=int, default=0, help="Parametr for filtering by quality"
    )

    return parser.parse_args()


def main():

    try:

        logger.info("Starting main function")

        args = parse_args()

        gc_bounds = parser_tuple(args.gc)
        length_bounds = parser_tuple(args.length)

        logger.info(f"Process start with {args.input}")
        logger.info(
            f"Parameters: GC={gc_bounds}, Length={length_bounds}, Quality≥{args.quality}"
        )

        filter_fastq(
            input_fastq=args.input,
            output_fastq=args.output,
            gc_bounds=gc_bounds,
            length_bounds=length_bounds,
            quality_threshold=args.quality,
        )

    except Exception as e:
        logger.critical(f"Failed with {e}")
        sys.exit(1)  # 1 - общепринятый код ошибки


if __name__ == "__main__":
    main()
