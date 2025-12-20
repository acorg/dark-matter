import argparse
import sys
from pathlib import Path

from dark.fasta import FastaReads
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


def extract_fasta(
    fasta: Path,
    ref_id: str,
    out_dir: Path,
    index_fasta: bool,
    quiet: bool,
    force: bool,
    keep_description: bool,
    e: Executor,
) -> int:
    new_reference = out_dir / f"{ref_id}.fasta"

    if new_reference.exists() and not force:
        sys.exit(
            f"Output FASTA file {str(new_reference)!r} already exists. Exiting. "
            "Re-run with --force to overwrite it."
        )

    # Find the reference in the FASTA and write it to the new_reference file.  In
    # writing it, make sure it has the same id as the reference we're going to look for
    # in the BAM file. We drop anything after the first space in the ID in the FASTA
    # file because that cannot appear in the BAM/SAM file we'll search (and write).
    for read in FastaReads(fasta):
        found = truncated = False
        if read.id == ref_id:
            found = True
            truncated = False
        elif not keep_description and read.id.split()[0] == ref_id:
            found = truncated = True

        if found:
            if truncated:
                if not quiet and not keep_description:
                    print(
                        f"Reference ID {ref_id!r} was found in {str(fasta)!r} but with "
                        f"a description ({read.id!r}). The description will be dropped "
                        "in the output FASTA so that it matches the reference name in "
                        "the output BAM (since spaces in an ID in a BAM file are "
                        "illegal). This is done so the IDs in both files match, as "
                        "required by some tools. To instead keep the description in "
                        "the output FASTA, re-run with --keep-description."
                    )
                read.id = ref_id

            with open(new_reference, "w") as fp:
                print(read.toString(), end="", file=fp)

            if index_fasta:
                e.execute(f"samtools faidx {str(new_reference)!r}")

            return len(read)

    if keep_description:
        # Have another look, dropping the description. Even though we were told not to,
        # pointing out that the ID is found but with a description is likely to be
        # helpful.
        for read in FastaReads(fasta):
            if read.id.split()[0] == ref_id:
                sys.exit(
                    f"A sequence with ID {ref_id!r} (exactly) does not appear in "
                    f"{str(fasta)!r}. However, that FASTA file has at least one "
                    f"sequence with an ID ({read.id!r}) whose first part matches your "
                    "reference ID. But, you specified --keep-description, so this ID "
                    "with its description has not been considered as matching. Maybe "
                    "try re-running without --keep-description?"
                )

    sys.exit(f"Could not find reference {ref_id!r} in {str(fasta)!r}.")


def extract_sam(
    bam: Path,
    ref_id: str,
    reference_len: int,
    out_dir: Path,
    keep_sam: bool,
    force: bool,
    e: Executor,
) -> int:
    new_sam = out_dir / f"{ref_id}.sam"
    new_bam = new_sam.with_suffix(".bam")

    for f in new_sam, new_bam:
        if f.exists() and not force:
            sys.exit(
                f"Output file {str(f)!r} already exists. Exiting. "
                "Re-run with --force to overwrite it."
            )

    if bam.name.lower().endswith(".bam"):
        # Convert the given BAM to SAM so we can easily read it (we could instead use
        # the pysam module).
        sam = out_dir / bam.with_suffix(".sam").name
        e.execute(f"samtools view -h -O SAM {str(bam)!r} --output {str(sam)!r}")
        our_sam = True
    else:
        if not bam.name.lower().endswith(".sam"):
            sys.exit(
                f"Input BAM/SAM file {str(bam)!r} does not have a .bam or .sam "
                "suffix so I don't know how to read it."
            )
        sam = bam
        our_sam = False

    read_count = 0

    with open(sam) as sam_in, open(new_sam, "w") as sam_out:
        for line in sam_in:
            line = line.strip()
            fields = line.split("\t", maxsplit=3)
            if line.startswith("@"):
                if fields[0] == "@SQ":
                    if fields[1] == f"SN:{ref_id}":
                        if fields[2] == f"LN:{reference_len}":
                            print(line, file=sam_out)
                        else:
                            sys.exit(
                                f"Found reference {ref_id!r} in "
                                f"{str(sam)!r} but the length field in the SAM "
                                f"file is {fields[2]!r} whereas the length of the "
                                f"reference in {str(ref_id)!r} is "
                                f"{reference_len}."
                            )
                else:
                    # A non-@SQ header line.
                    print(line, file=sam_out)
            else:
                # A non-header line. Print it if this is a match on the reference.
                if fields[2] == ref_id:
                    read_count += 1
                    print(line, file=sam_out)

    e.execute(f"samtools sort -O BAM {str(new_sam)!r} -o {str(new_bam)!r}")
    e.execute(f"samtools index {str(new_bam)!r}")

    # Remove the intermediate SAM, if we made it. This was just a direct conversion of
    # the input BAM that we made so we could read it easily to find the matching records
    # (which we wrote to new_sam).
    if our_sam:
        sam.unlink()

    if not keep_sam:
        new_sam.unlink()

    return read_count


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
