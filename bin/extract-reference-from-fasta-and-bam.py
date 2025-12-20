import argparse
from pathlib import Path

from dark.extract_reference import extract_fasta, extract_sam
from dark.process import Executor


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract a reference sequence from a FASTA file and its matches "
            "from a BAM file."
        ),
    )

    parser.add_argument(
        "--reference-id",
        required=True,
        metavar="ID",
        help="The ID of the reference sequence.",
    )

    parser.add_argument(
        "--out-dir",
        type=Path,
        required=True,
        metavar="DIR",
        help="The directory into which to put the output files.",
    )

    parser.add_argument(
        "--bam", "--sam", type=Path, required=True, help="The BAM or SAM file."
    )

    parser.add_argument("--fasta", type=Path, required=True, help="The FASTA file.")

    parser.add_argument("--quiet", action="store_true", help="Do not print a summary.")

    parser.add_argument(
        "--keep-description",
        action="store_true",
        help=(
            "Do not ignore the description (if any) of the sequence in the FASTA file. "
            "Normally, anything after the first space in the ID in the FASTA file will "
            "be ignored (for finding the ID in the FASTA file, in the output FASTA, "
            "and will not be included in the output BAM (where it is illegal). If you "
            "use this option, the description will not be ignored, so the "
            "--reference-id you give will be used in full. This may result in not "
            "finding what you are looking for (in the input FASTA)."
        ),
    )

    parser.add_argument(
        "--force", action="store_true", help="Overwrite pre-existing output files."
    )

    parser.add_argument(
        "--index-fasta",
        action="store_true",
        help="Index the extracted FASTA (using samtools faidx).",
    )

    parser.add_argument(
        "--keep-sam",
        action="store_true",
        help=(
            "Do not remove the intermediate SAM file with the matching records "
            "(produced en route to making a BAM file)"
        ),
    )

    return parser.parse_args()


def main() -> None:
    args = get_args()

    ref_id = args.reference_id
    out_dir = args.out_dir

    if not out_dir.exists():
        out_dir.mkdir()

    e = Executor()

    ref_len = extract_fasta(
        args.fasta,
        ref_id,
        out_dir,
        args.index_fasta,
        args.quiet,
        args.force,
        args.keep_description,
        e,
    )

    matches = extract_sam(
        args.bam, ref_id, ref_len, out_dir, args.keep_sam, args.force, e
    )

    if not args.quiet:
        s = "" if matches == 1 else "es"
        print(f"Found {matches} match{s} for {ref_id!r} in {str(args.bam)!r}.")


if __name__ == "__main__":
    main()
