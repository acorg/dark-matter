#!/usr/bin/env python
"""
Extract FASTA sequences from BWA index files (.pac, .ann, .amb).

Usage:
    python bwa_index_to_fasta.py <index_prefix> [-o output.fa]

The index prefix is the same value you'd pass to `bwa mem` — typically the
original reference filename (e.g. ref.fa), with .ann/.pac/.amb alongside it.
"""

import argparse
import sys
from pathlib import Path
from typing import TypedDict


class SequenceInfo(TypedDict):
    name: str
    anno: str
    offset: int
    length: int


def readAnn(annPath: str) -> tuple[int, list[SequenceInfo]]:
    """Parse .ann file to get sequence names, annotations, and lengths.

    Format (from bwa bntseq.c bns_dump):
        Line 1: <total_pac_length> <num_sequences> <seed>
        Per sequence:
            <gi> <name> [<annotation>]
            <offset> <length> <n_ambs>
    """
    sequences = []
    with open(annPath) as fp:
        header = fp.readline().split()
        totalLength = int(header[0])
        numSeqs = int(header[1])

        for _ in range(numSeqs):
            parts = fp.readline().rstrip("\n").split(None, 2)  # gi, name, [annotation]
            name = parts[1]
            # BWA writes "(null)" when there's no annotation in the FASTA header
            anno = parts[2] if len(parts) > 2 else ""
            if anno == "(null)":
                anno = ""

            parts2 = fp.readline().split()
            offset = int(parts2[0])
            length = int(parts2[1])

            sequences.append(
                {
                    "name": name,
                    "anno": anno,
                    "offset": offset,
                    "length": length,
                }
            )

    return totalLength, sequences


def readAmb(ambPath: str) -> list[tuple[int, int, str]]:
    """Parse .amb file to get positions of ambiguous (N) bases.

    Format:
        Line 1: <total_length> <num_sequences> <num_holes>
        Per hole:
            <offset> <length> <character>
    """
    holes: list[tuple[int, int, str]] = []
    with open(ambPath) as fp:
        header = fp.readline().split()
        numHoles = int(header[2])

        for _ in range(numHoles):
            parts = fp.readline().split()
            offset = int(parts[0])
            length = int(parts[1])
            char = parts[2]
            holes.append((offset, length, char))

    return holes


def readPac(pacPath: str) -> list[str]:
    """Decode the .pac file (2-bit packed sequence).

    Each byte stores 4 bases in the high-to-low bits:
        bits 7-6 = first base, bits 5-4 = second, etc.
        0=A, 1=C, 2=G, 3=T

    The very last byte of the file stores how many bases are
    encoded in the penultimate byte (0 means 4).
    """
    decode = ["A", "C", "G", "T"]

    with open(pacPath, "rb") as fp:
        data = fp.read()

    # Last byte = number of valid bases in the final data byte
    remainder = data[-1]
    pacData = data[:-1]

    seq: list[str] = []
    nBases = 4
    for i, byte in enumerate(pacData):
        # If on last byte
        if i == len(pacData) - 1:
            nBases = remainder if remainder > 0 else 4

        for j in range(nBases):
            # Shift the base we're trying to read to the rightmost two bits
            # and extract them with binary and b00000011
            base = (byte >> (6 - 2 * j)) & 0x3
            seq.append(decode[base])

    return seq


def applyAmbiguities(seq: list[str], holes: list[tuple[int, int, str]]) -> None:
    """Replace positions marked in .amb with their ambiguity character (usually N)."""
    for offset, length, char in holes:
        for i in range(length):
            pos = offset + i
            if pos < len(seq):
                seq[pos] = char


def reconstructFasta(
    sequences: list[SequenceInfo],
    holes: list[tuple[int, int, str]],
    seqData: list[str],
    lineWidth: int = 60,
) -> str:
    """Reconstruct FASTA from BWA index files and return as a string.

    Args:
        prefix: BWA index prefix (e.g. ref.fa — expects ref.fa.ann, ref.fa.pac,
                                  ref.fa.amb)
        lineWidth: Number of bases per line in output (default 60)

    Returns:
        The reconstructed FASTA as a string.
    """
    applyAmbiguities(seqData, holes)

    fastaContent: list[str] = []
    for s in sequences:
        header = s["name"]
        if s["anno"]:
            header += " " + s["anno"]
        fastaContent.append(f">{header}\n")

        start = s["offset"]
        end = start + s["length"]
        subseq = "".join(seqData[start:end])

        for i in range(0, len(subseq), lineWidth):
            fastaContent.append(subseq[i : i + lineWidth] + "\n")

    return "".join(fastaContent)


def createParser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract FASTA from BWA index files (.pac/.ann/.amb)"
    )
    parser.add_argument(
        "prefix",
        help=(
            "BWA index prefix (e.g. ref.fa — expects ref.fa.ann, ref.fa.pac, "
            "ref.fa.amb)"
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output FASTA file (default: stdout)",
        default=None,
    )
    parser.add_argument(
        "--lineWidth", default=60, help="Line width to use for FASTA file."
    )

    return parser


def main() -> None:
    parser = createParser()
    args = parser.parse_args()

    prefix = args.prefix
    annPath = f"{prefix}.ann"
    pacPath = f"{prefix}.pac"
    ambPath = f"{prefix}.amb"

    # Check files exist
    missingFiles: list[str] = []
    for path in [annPath, pacPath, ambPath]:
        if not Path(path).exists():
            missingFiles.append(path)

    if missingFiles:
        for path in missingFiles:
            print(f"Error: {path} not found", file=sys.stderr)
        sys.exit(1)

    # Read index components
    totalLength, sequences = readAnn(annPath)
    holes = readAmb(ambPath)
    seqData = readPac(pacPath)

    # Sanity check
    if len(seqData) < totalLength:
        print(
            f"Warning: .pac contains {len(seqData)} bases but .ann says {totalLength}",
            file=sys.stderr,
        )

    fastaStr = reconstructFasta(sequences, holes, seqData, args.lineWidth)

    # Write output
    if args.output:
        with open(args.output, "w") as fp:
            fp.write(fastaStr)
        print(
            f"Wrote {len(sequences)} sequence(s) ({totalLength} bp) to {args.output}",
            file=sys.stderr,
        )
    else:
        print(fastaStr)


if __name__ == "__main__":
    main()
