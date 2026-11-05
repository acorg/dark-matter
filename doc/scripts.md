# dark-matter command-line scripts

Every script described below accepts a `--help` option, which prints a
full description of what it does along with its complete list of
command-line options. The one-line summaries below are necessarily
brief — if you're not sure exactly how a script behaves,
`script.py --help` is the place to check first.

## Running a script

**From a clone of this repository**, use [uv](https://docs.astral.sh/uv/)
from the repository's top-level directory:

```sh
uv run bin/fasta-lengths.py --help
```

This works the same way on Linux, macOS, and Windows.

**If you've installed dark-matter** (e.g., via `pip install dark-matter` or
`uv add dark-matter`, as described in the top-level `README.md`), the
scripts are installed into your Python environment's `bin` directory
alongside the library itself. If you're using uv, you don't need to find
that directory yourself — from within a project that depends on
dark-matter, `uv run` will find and run an installed script by name, on
any platform:

```sh
uv run fasta-lengths.py --help
```

Without uv, on Linux and macOS, that `bin` directory is normally on your
shell's `PATH` already (e.g., once you've activated a virtual environment
with `source .venv/bin/activate`), so you can just run a script by name:

```sh
fasta-lengths.py --help
```

On Windows, the installed scripts are plain `.py` files rather than `.exe`
launchers, so the reliable way to run one — regardless of your `PATH`
setup — is to invoke it explicitly with Python:

```bat
python path\to\your\venv\bin\fasta-lengths.py --help
```

(use `py` in place of `python` if that's what's set up on your system).

**If you don't want to install anything**, not even into a project's
environment, use `uvx` (part of uv) to run a script straight from PyPI in
a disposable, one-off environment:

```sh
uvx --from dark-matter fasta-lengths.py --help
```

This downloads dark-matter and its dependencies into a temporary
environment just for this command and leaves nothing behind. The
`--from dark-matter` part is required — `uvx` can only infer the package
to install when a script's name matches its package's name, which isn't
the case for any of dark-matter's scripts.

---


## 1. FASTA/FASTQ: general read inspection, filtering, and manipulation

<details><summary><strong><code>fasta-ids.py</code></strong> — Given FASTA on stdin, write the sequence ids to stdout. <span title="Show --help">ℹ️</span></summary>

```
usage: fasta-ids.py [-h] [--fastaFile FILENAME] [--readClass CLASSNAME]
                    [--fasta | --fastq | --fasta-ss]
                    [FILENAME]

Given FASTA on stdin, write the sequence ids to stdout.

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given.

options:
  -h, --help            show this help message and exit
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead.
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ).
  --fastq               If specified, input will be treated as FASTQ.
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly.
```

</details>

<details><summary><strong><code>fasta-lengths.py</code></strong> — Given FASTA on stdin, write the sequence ids and lengths to stdout. <span title="Show --help">ℹ️</span></summary>

```
usage: fasta-lengths.py [-h] [--fastaFile FILENAME] [--readClass CLASSNAME]
                        [--fasta | --fastq | --fasta-ss]
                        [FILENAME]

Given FASTA on stdin, write the sequence ids and lengths to stdout.

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given.

options:
  -h, --help            show this help message and exit
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead.
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ).
  --fastq               If specified, input will be treated as FASTQ.
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly.
```

</details>

<details><summary><strong><code>fasta-count.py</code></strong> — Given FASTA on stdin, print the number of sequences to stdout. <span title="Show --help">ℹ️</span></summary>

```
usage: fasta-count.py [-h] [--fastaFile FILENAME] [--readClass CLASSNAME]
                      [--fasta | --fastq | --fasta-ss]
                      [FILENAME]

Given FASTA on stdin, print the number of sequences to stdout.

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given. (default: None)

options:
  -h, --help            show this help message and exit
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
                        (default: None)
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead. (default: DNARead)
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ). (default: False)
  --fastq               If specified, input will be treated as FASTQ.
                        (default: False)
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly. (default: False)
```

</details>

<details><summary><strong><code>fasta-count-chars.py</code></strong> — Given FASTA on stdin, write information about the counts of characters in
the sequences to stdout. <span title="Show --help">ℹ️</span></summary>

```
usage: fasta-count-chars.py [-h] [--ignoreCase] [--sep SEP] [--omitIds]
                            [--omitHeader] [--ambiguous] [--omitZeroes]
                            [--locations] [--sort] [--python] [--chars CHARS]
                            [--fastaFile FILENAME] [--readClass CLASSNAME]
                            [--fasta | --fastq | --fasta-ss]
                            [FILENAME]

Given FASTA on stdin, write information about the counts of characters in the
sequences to stdout.

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given.

options:
  -h, --help            show this help message and exit
  --ignoreCase, -i      Map all characters to uppercase.
  --sep SEP             The output field separator character.
  --omitIds             Do not print sequence ids as the first field.
  --omitHeader          Do not print a header.
  --ambiguous           If --chars is not given, print counts for all
                        ambiguous DNA codes.
  --omitZeroes          Do not print chars whose counts are zero.
  --locations           Also print the 1-based genome locations of the matched
                        characters.
  --sort                Sort the characters for each sequence by decreasing
                        count (with draws broken in favour of alphabetically
                        increasing character).
  --python              Print 0-based Python closed-open location ranges
                        (i.e., the upper bound is not included) when
                        --locations is used.
  --chars, -c CHARS     The characters to count. Default is all IUPAC
                        nucleotide codes.
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead.
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ).
  --fastq               If specified, input will be treated as FASTQ.
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly.
```

</details>

<details><summary><strong><code>fasta-join.py</code></strong> — Given FASTA on stdin, write the ids, sequences, and quality strings on a
single TAB-separated line to stdout (used internally by <code>fasta-diff.sh</code>). <span title="Show --help">ℹ️</span></summary>

```
usage: fasta-join.py [-h] [--separator SEPARATOR] [--removeIds]
                     [--fastaFile FILENAME] [--readClass CLASSNAME] [--fasta |
                     --fastq | --fasta-ss]
                     [FILENAME]

Given FASTA on stdin, write the ids, sequences, and quality strings on a
single TAB-separated line to stdout.

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given.

options:
  -h, --help            show this help message and exit
  --separator SEPARATOR
                        The character string to separate ids from sequences
                        and quality strings (if any)
  --removeIds           Do not print sequence ids
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead.
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ).
  --fastq               If specified, input will be treated as FASTQ.
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly.
```

</details>

<details><summary><strong><code>fasta-sort.py</code></strong> — Write sorted FASTA/FASTQ to stdout. Sorting is by sequence id, then by
sequence, then by quality (if FASTQ). <span title="Show --help">ℹ️</span></summary>

```
usage: fasta-sort.py [-h] [--regex REGEX] [--reverse] [--length]
                     [--fastaFile FILENAME] [--readClass CLASSNAME] [--fasta |
                     --fastq | --fasta-ss]
                     [FILENAME]

Write sorted FASTA/Q to stdout. Sorting is by sequence id, then by sequence,
then by quality (if FASTQ)

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given.

options:
  -h, --help            show this help message and exit
  --regex REGEX         A regular expression to specify a region (or regions)
                        to extract from sequence ids to sort on. Regions will
                        be converted to integers if possible.
  --reverse, -r         Sort by decreasing value.
  --length, -l          Sort by sequence length.
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead.
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ).
  --fastq               If specified, input will be treated as FASTQ.
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly.
```

</details>

<details><summary><strong><code>fasta-split.py</code></strong> — Split sequences in a FASTA/FASTQ file into numerically named files, each
containing a given number of sequences. <span title="Show --help">ℹ️</span></summary>

```
usage: fasta-split.py [-h] [--outDir OUTDIR] [--verbose] [--count COUNT]
                      [--noLeadingZeroes] [--force]
                      [--saveAs {fasta,fastq,fasta-ss}] [--fastaFile FILENAME]
                      [--readClass CLASSNAME] [--fasta | --fastq | --fasta-ss]
                      [FILENAME]

Split sequences in a FASTA/Q file into numerically named files, each
containing a given number of sequences.

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given. (default: None)

options:
  -h, --help            show this help message and exit
  --outDir OUTDIR       The directory to make the files in. (default: .)
  --verbose             If given, print sequence ids as they are processed.
                        (default: False)
  --count COUNT         The number of sequences to put in each file. (default:
                        100)
  --noLeadingZeroes     If given, numeric filenames will not have leading
                        zeroes. (default: False)
  --force               If given, overwrite pre-existing files. (default:
                        False)
  --saveAs {fasta,fastq,fasta-ss}
                        The output format. The default is to match the input
                        format, so there is usually no need to specify this
                        option. It can be used to force conversion from FASTQ
                        to FASTA (default: None)
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
                        (default: None)
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead. (default: DNARead)
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ). (default: False)
  --fastq               If specified, input will be treated as FASTQ.
                        (default: False)
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly. (default: False)
```

</details>

<details><summary><strong><code>fasta-split-by-id.py</code></strong> — Split sequences in a FASTA file into separate files, named either by their
sequence id or numerically. <span title="Show --help">ℹ️</span></summary>

```
usage: fasta-split-by-id.py [-h] [--outDir OUTDIR] [--verbose] [--numeric]
                            [--noLeadingZeroes]
                            [--filenamePrefix FILENAMEPREFIX]
                            [--filenameSuffix FILENAMESUFFIX]
                            [--idPrefix IDPREFIX] [--idSuffix IDSUFFIX]
                            [--force] [--dryRun]
                            [--saveAs {fasta,fastq,fasta-ss}]
                            [--fastaFile FILENAME] [--readClass CLASSNAME]
                            [--fasta | --fastq | --fasta-ss]
                            [FILENAME]

Split sequences in a FASTA file into separate files, named either by their
sequence id or numerically.

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given. (default: None)

options:
  -h, --help            show this help message and exit
  --outDir OUTDIR       The directory to make the files in. (default: .)
  --verbose, -v         If given, print sequence ids as they are processed.
                        (default: False)
  --numeric             If given, use numeric filenames, like 1.fasta,
                        2.fasta, etc. (default: False)
  --noLeadingZeroes     If given, numeric filenames will not have leading
                        zeroes. (default: False)
  --filenamePrefix FILENAMEPREFIX
                        Use this as a prefix in output filenames. (default: )
  --filenameSuffix FILENAMESUFFIX
                        Use this as a suffix (coming before '.fasta') in
                        output filenames. (default: )
  --idPrefix IDPREFIX   Use this as a prefix for sequence IDs. (default: )
  --idSuffix IDSUFFIX   Use this as a suffix for sequence IDs. (default: )
  --force, -f           If given, overwrite pre-existing files. (default:
                        False)
  --dryRun, -n          Don't create files, just print what would be done.
                        Implies --verbose. (default: False)
  --saveAs {fasta,fastq,fasta-ss}
                        The output format. The default is to match the input
                        format, so there is usually no need to specify this
                        option. It can be used to force conversion from FASTQ
                        to FASTA. (default: None)
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
                        (default: None)
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead. (default: DNARead)
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ). (default: False)
  --fastq               If specified, input will be treated as FASTQ.
                        (default: False)
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly. (default: False)
```

</details>

<details><summary><strong><code>fasta-split-by-length.py</code></strong> — Split FASTA/FASTQ sequences longer than a given length; print the
beginning and end sequences of sufficiently long sequences separately. <span title="Show --help">ℹ️</span></summary>

```
usage: fasta-split-by-length.py [-h] --length LENGTH [--leftPrefix LEFTPREFIX]
                                [--leftSuffix LEFTSUFFIX]
                                [--rightPrefix RIGHTPREFIX]
                                [--rightSuffix RIGHTSUFFIX] [--noWarn]
                                [--fastaFile FILENAME] [--readClass CLASSNAME]
                                [--fasta | --fastq | --fasta-ss]
                                [FILENAME]

Split FASTA/Q sequences longer than a given length. Print the beginning and
end sequences of sufficiently long sequences separately.

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given. (default: None)

options:
  -h, --help            show this help message and exit
  --length LENGTH       Sequences longer than this will be split in two.
  --leftPrefix LEFTPREFIX
                        The string to put before the id of the starting (left)
                        subsequence. (default: )
  --leftSuffix LEFTSUFFIX
                        The string to put after the id of the starting (left)
                        subsequence. (default: )
  --rightPrefix RIGHTPREFIX
                        The string to put before the id of the ending (right)
                        subsequence. (default: )
  --rightSuffix RIGHTSUFFIX
                        The string to put after the id of the ending (right)
                        subsequence. (default: )
  --noWarn              Do not warn the user if no prefixes are given (which
                        can result in FASTA output with duplicated IDs.
                        (default: True)
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
                        (default: None)
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead. (default: DNARead)
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ). (default: False)
  --fastq               If specified, input will be treated as FASTQ.
                        (default: False)
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly. (default: False)
```

</details>

<details><summary><strong><code>fasta-subset.py</code></strong> — Extract a subset of FASTA reads by id. <span title="Show --help">ℹ️</span></summary>

```
usage: fasta-subset.py [-h] [--readIdFile READIDFILE] [ids ...]

Extract a subset of FASTA reads by id

positional arguments:
  ids                   Wanted read ids.

options:
  -h, --help            show this help message and exit
  --readIdFile READIDFILE
                        A file of wanted read ids.

Given a set of FASTA sequence identifiers from sys.argv or in a file, read
FASTA from stdin, and print FASTA to stdout for the given sequence ids.
```

</details>

<details><summary><strong><code>fasta-subtraction.py</code></strong> — Given two or more FASTA file names as arguments, remove from the first
file any read whose id also appears in a later file, and print what's
left. <span title="Show --help">ℹ️</span></summary>

```
usage: fasta-subtraction.py [-h] FASTA-file [FASTA-file ...]

Given two or more FASTA file names as arguments, remove from the first file
any read whose id also appears in a later file, and print what's left.

positional arguments:
  FASTA-file  A FASTA file name. Reads are removed from the first file named
              if their id also appears in any of the later files.

options:
  -h, --help  show this help message and exit
```

</details>

<details><summary><strong><code>filter-fasta.py</code></strong> — Given FASTA on stdin and a set of filtering criteria, write filtered
FASTA to stdout. Has many options (see <code>--help</code>); the <code>bin/README.md</code>
documents a couple of non-obvious uses (e.g., reformatting multi-line
FASTA to one-line-per-sequence). <span title="Show --help">ℹ️</span></summary>

```
usage: filter-fasta.py [-h] [--quiet] [--saveAs {fasta,fastq,fasta-ss}]
                       [--checkResultCount CHECKRESULTCOUNT]
                       [--fastaFile FILENAME] [--readClass CLASSNAME]
                       [--fasta | --fastq | --fasta-ss] [--minLength N]
                       [--maxLength N] [--maxNFraction N]
                       [--whitelist SEQUENCE-ID] [--blacklist SEQUENCE-ID]
                       [--whitelistFile SEQUENCE-ID-FILE]
                       [--blacklistFile SEQUENCE-ID-FILE] [--titleRegex REGEX]
                       [--negativeTitleRegex REGEX]
                       [--sequenceWhitelist SEQUENCE]
                       [--sequenceBlacklist SEQUENCE]
                       [--sequenceWhitelistFile SEQUENCE-FILE]
                       [--sequenceBlacklistFile SEQUENCE-FILE]
                       [--sequenceRegex REGEX] [--sequenceNegativeRegex REGEX]
                       [--keepSequences NUMBER,RANGE,... |
                       --removeSequences NUMBER,RANGE,...] [--head N]
                       [--removeDuplicates | --removeDuplicatesById]
                       [--removeDuplicatesUseMD5] [--randomSubset N]
                       [--trueLength N] [--sampleFraction FRACTION]
                       [--sequenceNumbersFile FILENAME]
                       [--keepSites KEEPSITES |
                       --keepSitesFile KEEPSITESFILE |
                       --removeSites REMOVESITES |
                       --removeSitesFile REMOVESITESFILE] [--removeGaps]
                       [--truncateTitlesAfter STR] [--removeDescriptions]
                       [--idLambda LAMBDA-FUNCTION]
                       [--readLambda LAMBDA-FUNCTION] [--reverse]
                       [--reverseComplement] [--rotate N]
                       [--setSite site:base[:quality] [site:base[:quality] ...]]
                       [--setSiteZeroBased] [--upper | --lower] [--upperId |
                       --lowerId] [--setQuality Q]
                       [FILENAME]

Given FASTA on stdin and a set of filtering criteria write filtered FASTA to
stdout.

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given.

options:
  -h, --help            show this help message and exit
  --quiet               If True, do not print the final sequence summary.
  --saveAs {fasta,fastq,fasta-ss}
                        The output format. The default is to match the input
                        format, so there is usually no need to specify this
                        option. It can be used to force conversion from FASTQ
                        to FASTA
  --checkResultCount CHECKRESULTCOUNT
                        The number of reads expected in the output. If this
                        number is not seen, the script exits with status 1 and
                        an error message is printed unless --quiet was used.
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead.
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ).
  --fastq               If specified, input will be treated as FASTQ.
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly.
  --minLength N         The minimum sequence length.
  --maxLength N         The maximum sequence length.
  --maxNFraction N      The maximum fraction of Ns that can be present in the
                        sequence.
  --whitelist SEQUENCE-ID
                        Sequence titles (ids) that should be whitelisted. May
                        be repeated.
  --blacklist SEQUENCE-ID
                        Sequence titles (ids) that should be blacklisted. May
                        be repeated.
  --whitelistFile SEQUENCE-ID-FILE
                        The name of a file that contains sequence titles (ids)
                        that should be whitelisted, one per line
  --blacklistFile SEQUENCE-ID-FILE
                        The name of a file that contains sequence titles (ids)
                        that should be blacklisted, one per line.
  --titleRegex REGEX    A regex that sequence titles (ids) must match.
  --negativeTitleRegex REGEX
                        A regex that sequence titles (ids) must not match.
  --sequenceWhitelist SEQUENCE
                        Sequences that should be whitelisted. May be repeated.
  --sequenceBlacklist SEQUENCE
                        Sequences that should be blacklisted. May be repeated.
  --sequenceWhitelistFile SEQUENCE-FILE
                        The name of a file that contains sequences that should
                        be whitelisted, one per line.
  --sequenceBlacklistFile SEQUENCE-FILE
                        The name of a file that contains sequences that should
                        be blacklisted, one per line.
  --sequenceRegex REGEX
                        A regex that sequences must match.
  --sequenceNegativeRegex REGEX
                        A regex that sequences must not match.
  --keepSequences NUMBER,RANGE,...
                        Specify (1-based) ranges of sequence numbers that
                        should be kept. E.g., --keepSequences 1-3,5 will
                        output just the 1st, 2nd, 3rd, and 5th sequences. All
                        others will be omitted.
  --removeSequences NUMBER,RANGE,...
                        Specify (1-based) ranges of sequence numbers that
                        should be removed. E.g., --removeSequences 1-3,5 will
                        output all but the 1st, 2nd, 3rd, and 5th sequences.
                        All others will be ouput.
  --head N              Only the first N sequences will be printed.
  --removeDuplicates    Duplicate reads will be removed, based only on
                        sequence identity. The first occurrence is kept.
  --removeDuplicatesById
                        Duplicate reads will be removed, based only on read
                        id. The first occurrence is kept.
  --removeDuplicatesUseMD5
                        MD5 sums will be stored instead of the full sequence
                        or read id when either --removeDuplicates or
                        removeDuplicatesById are given. One of
                        --removeDuplicates and --removeDuplicatesById must
                        also be given. Note that this makes duplicate removal
                        probabilistic. This option can be used to reduce the
                        amount of RAM consumed during duplicate removal.
  --randomSubset N      An integer giving the number of sequences that should
                        be kept. These will be selected at random.
  --trueLength N        The number of reads in the FASTA input. Only to be
                        used with randomSubset
  --sampleFraction FRACTION
                        A [0.0, 1.0] C{float} indicating a fraction of the
                        reads that should be allowed to pass through the
                        filter. The sample size will only be approximately the
                        product of the sample fraction and the number of
                        reads. The sample is taken at random.
  --sequenceNumbersFile FILENAME
                        A file of (1-based) sequence numbers to retain.
                        Numbers must be one per line.
  --keepSites, --keepIndices KEEPSITES
                        Specify 1-based sequence sites to keep. All other
                        sites will be removed. The sites must be given in the
                        form e.g., 24,100-200,260. Note that the requested
                        sites will be taken from the input sequences in order,
                        not in the order given by --keepSites. I.e.,
                        --keepSites 5,8-10 will get you the same result as
                        --keepSites 8-10,5.
  --keepSitesFile, --keepIndicesFile KEEPSITESFILE
                        Specify a file containing 1-based sites to keep. All
                        other sequence sites will be removed. Lines in the
                        file must be given in the form e.g., 24,100-200,260.
                        See --keepSites for more detail.
  --removeSites, --removeIndices REMOVESITES
                        Specify 1-based sites to remove. All other sequence
                        sites will be kept. The sites must be given in the
                        form e.g., 24,100-200,260. See --keepSites for more
                        detail.
  --removeSitesFile, --removeIndicesFile REMOVESITESFILE
                        Specify a file containing 1-based sites to remove. All
                        other sequence sites will be kept. Lines in the file
                        must be given in the form e.g., 24,100-200,260. See
                        --keepSites for more detail.
  --removeGaps          If True, gap ('-') characters in sequences will be
                        removed.
  --truncateTitlesAfter STR
                        A string that sequence titles (ids) will be truncated
                        beyond. If the truncated version of a title has
                        already been seen, that title will be skipped.
  --removeDescriptions  Read id descriptions will be removed. The description
                        is the part of a sequence id after the first
                        whitespace (if any).
  --idLambda LAMBDA-FUNCTION
                        A one-argument function taking and returning a read
                        id. E.g., --idLambda "lambda id: id.split('_')[0]" or
                        --idLambda "lambda id: id[:10]". If the function
                        returns None, the read will be filtered out.
  --readLambda LAMBDA-FUNCTION
                        A one-argument function taking and returning a read.
                        E.g., --readLambda "lambda r: Read(r.id.split('_')[0],
                        r.sequence.strip('-'))". Make sure to also modify the
                        quality string if you change the length of a FASTQ
                        sequence. If the function returns None, the read will
                        be filtered out. The function will be passed to eval
                        with the dark.reads classes Read, DNARead, AARead,
                        etc. all in scope.
  --reverse             Reverse the sequences. Note that this is NOT reverse
                        complementing.
  --reverseComplement   Reverse complement the sequences.
  --rotate N            Rotate sequences. A positive integer rotates the
                        sequence to the right, negative to the left
  --setSite site:base[:quality] [site:base[:quality] ...]
                        Set a (1-based) site to a certain base (with a
                        quality, if FASTQ). May be given with multiple
                        specifications to change several sites, e.g.:
                        --setSite 4:A 20:G Several things should be noted: 1)
                        because the sites are 1-based, you cannot use -1 as a
                        site to change the final site, you need to use 0
                        instead (since 1 is subtracted from the site indices
                        to get Python's 0-based offsets). If you find this
                        weird, use --setSiteZeroBased and give 0-based sites.
                        2) To pass a negative site you should give a leading
                        space so the argparse parser doesn't think your offset
                        is the start of another argument. Therefore: --setSite
                        ' -1:A' can be used. 3) These changes are applied
                        before any other read-altering operations (e.g.,
                        sequence rotation, applying a read lambda, reversing
                        or reverse complementing).
  --setSiteZeroBased    Indicate that the sites given in --setSite are zero
                        based.
  --upper               Convert sequences to uppercase.
  --lower               Convert sequences to lowercase.
  --upperId             Convert sequence IDs to uppercase.
  --lowerId             Convert sequence IDs to lowercase.
  --setQuality Q        Set all quality characters to 'Q'. Note that this is a
                        string value, not an integer. As well as changing
                        FASTQ qualities, this can also be used to convert
                        FASTA to FASTQ (with a constant quality for all bases)
                        when used with --saveAs fastq. A PHRED quality of 20
                        can be given as '5', of 30 as '?', and 32 as 'A'. See
                        man ascii for other values.
```

</details>

<details><summary><strong><code>filter-fasta-by-complexity.py</code></strong> — Filter FASTA for low-complexity sequences. <span title="Show --help">ℹ️</span></summary>

```
usage: filter-fasta-by-complexity.py [-h] [--maxLowComplexity FRAC]
                                     [--rejectFile FILE]
                                     fasta

Filter FASTA for low complexity sequences

positional arguments:
  fasta                 The FASTA file to read.

options:
  -h, --help            show this help message and exit
  --maxLowComplexity FRAC
                        Sequences with a fraction of low-complexity bases up
                        to this level will be filtered out.
  --rejectFile FILE     The name of a file to save rejected sequences to.

Read DNA FASTA from stdin and print FASTA to stdout, excluding sequences of
low complexity.
```

</details>

<details><summary><strong><code>filter-fasta-by-taxonomy.py</code></strong> — Filter FASTA based on taxonomy, matched against a taxonomic regular
expression. <strong>Note:</strong> only works on FASTA with old-style NCBI <code>gi|...|</code>
description lines, a numbering scheme NCBI retired in 2016 — likely a
candidate for the "probably obsolete" pile (see bottom of this document). <span title="Show --help">ℹ️</span></summary>

```
usage: filter-fasta-by-taxonomy.py [-h] --taxonomy TAXONOMY [--invert]
                                   [--detailsFile FILE]

Filter FASTA based on taxonomy

options:
  -h, --help           show this help message and exit
  --taxonomy TAXONOMY  The regex to match the taxonomy on. Case is ignored.
  --invert             If True, only write sequences whose taxonomy does not
                       match.
  --detailsFile FILE   The name of a file to save taxonomy details to

Read DNA FASTA from stdin and print FASTA to stdout, only including sequences
that match the requested taxonomy regular expression at any level.
```

</details>

<details><summary><strong><code>format-fasta.py</code></strong> — Given FASTA on stdin, write it to stdout according to a given output format. <span title="Show --help">ℹ️</span></summary>

```
usage: format-fasta.py [-h] [--prefix PREFIX] [--start START] [--end END]
                       [--format FORMAT] [--fastaFile FILENAME]
                       [--readClass CLASSNAME] [--fasta | --fastq |
                       --fasta-ss]
                       [FILENAME]

Given FASTA on stdin, write it to stdout according to a given format.

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given.

options:
  -h, --help            show this help message and exit
  --prefix PREFIX       The string prefix to be used. E.g., use --prefix r for
                        raw strings. The default is "f" (i.e., use Python
                        f-strings) for Python 3.6 and later, or else "" (i.e.,
                        regular strings).
  --start START         The integer to start counting from, for the "count"
                        variable.
  --end END             A string to use as the "end" keyword argument for
                        print(). This will be passed to eval, allowing
                        convenient use of Python string constructs such as \t
                        on the command line. If not specified, a newline will
                        be printed.
  --format FORMAT       The output format. The default is to reproduce the
                        input. The format can be that used by Python to
                        specify items in a dictionary. Available dictionary
                        items are "id", "sequence", "quality", "length",
                        "count", and "md5" (of the sequence). A
                        dark.reads.Read instance (named 'read') and an 'md5'
                        function (that takes and returns a string) are both
                        also in scope for each FASTA record and these can be
                        accessed using f-string format (with Python 3.6 or
                        greater). The passed format string will be eval'd in a
                        double quoted string, as with the --end string,
                        allowing the use of Python string conveniences such as
                        \t on the command line.
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead.
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ).
  --fastq               If specified, input will be treated as FASTQ.
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly.
```

</details>

<details><summary><strong><code>fasta-diff.py</code></strong> — Diff two FASTA sequences with visual highlighting (this is the
Python-based, within-file diff; compare with the shell-based
<code>fasta-diff.sh</code> below, which diffs whole FASTA/FASTQ <em>files</em>). <span title="Show --help">ℹ️</span></summary>

```
usage: fasta-diff.py [-h] [--name1 NAME1] [--name2 NAME2] [-c CONTEXT]
                     [-w WIDTH] [--noColor] [--noTable] [--maxDiffs MAXDIFFS]
                     [--full] [--aligner {edlib,mafft,needle}]
                     [--alignerOptions ALIGNEROPTIONS] [--threads THREADS]
                     [--verbose]
                     file1 file2

Diff two FASTA sequences with visual highlighting.

positional arguments:
  file1                 First FASTA file.
  file2                 Second FASTA file.

options:
  -h, --help            show this help message and exit
  --name1 NAME1         Sequence name to use from file1 (default: first).
                        (default: None)
  --name2 NAME2         Sequence name to use from file2 (default: first).
                        (default: None)
  -c, --context CONTEXT
                        Bases of context around each difference. (default: 10)
  -w, --width WIDTH     Alignment display width. (default: 80)
  --noColor             Disable colored output. (default: False)
  --noTable             Don't print the per-position difference table.
                        (default: False)
  --maxDiffs MAXDIFFS   The maximum number of differences for which the diff
                        table and alignment are printed. If the number of
                        differences exceeds this value, only the summary is
                        printed. Use -1 (the default) to always print
                        differences. (default: 20)
  --full                Show full alignment (no context windowing). (default:
                        False)
  --aligner {edlib,mafft,needle}
                        The alignment algorithm to use for sequences of
                        different length. (default: edlib)
  --alignerOptions ALIGNEROPTIONS
                        Optional arguments to pass to the alignment algorithm.
                        If the aligner is 'mafft', the default options are '--
                        auto --maxiterate 1000 --preservecase'. If 'needle',
                        the default is 'auto'. Do not try to set the number of
                        threads here; use the --threads argument instead. If
                        you are using mafft, see https://mafft.cbrc.jp/alignme
                        nt/software/algorithms/algorithms.html for some
                        possible option combinations. (default: None)
  --threads THREADS     The number of threads to use when running the aligner
                        (if the alignment algorithm can make use of multiple
                        threads; mafft can, needle and edlib cannot).
                        (default: 12)
  --verbose             Print alignment progress to standard error. (default:
                        False)
```

</details>

**`fasta-diff.sh`** — A different, cruder FASTA/FASTQ diff: converts both files' records to
single sorted TAB-separated lines and runs `diff`, then converts the
result back into FASTA/FASTQ-like output. Good for a quick "are these
basically the same" check, not for pinpointing exact differences.

<details><summary><strong><code>fasta-find.py</code></strong> — Given FASTA on stdin and a sequence (or genome region) to find, write
information about (potentially overlapping) matches, optionally including
reversed/complemented/reverse-complemented matches. <span title="Show --help">ℹ️</span></summary>

```
usage: fasta-find.py [-h] (--pattern PATTERN | --region REGION)
                     [--reverseComplement] [--reverse] [--complement] [--end]
                     [--offsetsOnly] [--checkMatchCount CHECKMATCHCOUNT]
                     [--caseSensitive] [--printSequences]
                     [--printTrivialMatch] [--sep SEP] [--zeroBased]
                     [--gapCharacter GAPCHARACTER] [--ignoreGaps]
                     [--fastaFile FILENAME] [--readClass CLASSNAME] [--fasta |
                     --fastq | --fasta-ss]
                     [FILENAME]

Given FASTA on stdin and a sequence (or genome region) to find, write
information about (potentially overlapping) matches (and, optionally,
reversed, complemented, or reversed complement matches), to standard output.
Offsets are inclusive and 1-based unless --zeroBased is used (in which case
they are 0-based closed/open offsets as used in Python).

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given.

options:
  -h, --help            show this help message and exit
  --pattern, --sequence PATTERN
                        The sequence to look for.
  --region REGION       A region of the form A-B, where A and B are numeric
                        sites, with A <= B. A and B will be treated as 1-based
                        offsets unless --zeroBased is used. This can be used
                        to find repeat copies of an input sequence region.
                        I.e., for each input sequence, the specified region of
                        the sequence will be searched for in the rest of the
                        sequence.
  --reverseComplement, --rc
                        Also search for the reverse complement. Reverse
                        complement matches are indicated with an (rc) suffix.
  --reverse             Also search for the reversed sequence. Reverse matches
                        are indicated with an (r) suffix.
  --complement          Also search for the complemented sequence. Complement
                        matches are indicated with a (c) suffix.
  --end                 Print the index of the end of the match, not the
                        start.
  --offsetsOnly         Only print the matching offset (or offsets, separated
                        by the --sep value).
  --checkMatchCount CHECKMATCHCOUNT
                        The expected match count. This can be used to force an
                        error if the number of occurrences of the search
                        pattern found is not as expected.
  --caseSensitive       Do not ignore upper/lower case when matching.
  --printSequences      Print the sequences searched for. Reversed,
                        complemented, and reverse-complemented sequences are
                        indicated by (r), (c), and (rc), respectively.
  --printTrivialMatch   If --region is used, there is a guaranteed match at
                        that location. Do not include that match in the
                        printed matches unless this option is given.
  --sep SEP             The output field separator string (default=TAB).
  --zeroBased           Return Python-style zero-based closed/open match
                        offsets. When --end is used, the output offset will be
                        for the character that follows the end of the match
                        (i.e., it is not part of the match).
  --gapCharacter GAPCHARACTER
                        The gap character (only used if --ignoreGaps is used).
  --ignoreGaps          Ignore gaps when matching. This is useful for finding
                        a subsequence of an aligned sequence. Note that the
                        reported index into input sequences will be for the
                        (possibly gapped) original input sequence.
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead.
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ).
  --fastq               If specified, input will be treated as FASTQ.
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly.
```

</details>

<details><summary><strong><code>fasta-match-offsets.py</code></strong> — Given FASTA on stdin and a sequence pattern, write information about
which sequences the pattern matches, and where. <span title="Show --help">ℹ️</span></summary>

```
usage: fasta-match-offsets.py [-h] --regex REGEX [--noUpper] [--overlapping]
                              [--noId] [--noEnd] [--noSequence]
                              [--noOverlapping] [--zeroBased] [--count]
                              [--fastaFile FILENAME] [--readClass CLASSNAME]
                              [--fasta | --fastq | --fasta-ss]
                              [FILENAME]

Given FASTA on stdin and a sequence pattern, write information about which
sequences and where the pattern matches to stdout.

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given.

options:
  -h, --help            show this help message and exit
  --regex REGEX         The regular expression pattern to match. For matching,
                        input sequences will be converted to uppercase unless
                        --noUpper is used.
  --noUpper             Do not convert sequences to uppercase for the
                        matching.
  --overlapping         Include overlapping matches.
  --noId                Do not print sequence ids.
  --noEnd               Do not print the end offset of matches.
  --noSequence          Do not print sequence match regions.
  --noOverlapping       Do not print a '+' indicator to identify matches that
                        overlap the previously printed one.
  --zeroBased           Print half-open zero-based offsets, as in Python.
  --count, -c           Just print the count of the number of matches, if the
                        count is non-zero.
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead.
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ).
  --fastq               If specified, input will be treated as FASTQ.
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly.
```

</details>

<details><summary><strong><code>fasta-base-indices.py</code></strong> — Given FASTA on stdin and a set of bases to look for, write the 1-based
indices of where those bases occur in each sequence (or, with <code>--any</code>,
where they occur in any sequence). <span title="Show --help">ℹ️</span></summary>

```
usage: fasta-base-indices.py [-h] [--bases BASES] [--matchCase] [--any]
                             [--fastaFile FILENAME] [--readClass CLASSNAME]
                             [--fasta | --fastq | --fasta-ss]
                             [FILENAME]

Given FASTA on standard input, and bases to look for, write the 1-based
indices of where the bases occur in all sequences to standard output. If
standard output is a terminal, these will be space separated, else newline
separated. Use --any to print indices that match any sequence.

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given.

options:
  -h, --help            show this help message and exit
  --bases BASES         The sequence bases whose indices should be printed. If
                        not specified, this will be the defined set of bases
                        for the input sequence type (i.e., "ACGT" for DNA).
                        This will have the effect of printing the indices for
                        which any sequence has an ambiguous or missing base.
  --matchCase           If specified, sequence case will be considered in
                        matching.
  --any                 If specified, print indices of bases that match any
                        sequence, otherwise (the default) indices are only
                        printed if they match all sequences.
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead.
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ).
  --fastq               If specified, input will be treated as FASTQ.
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly.

This can be used to find all columns of a FASTA multiple sequence alignment
that contain gaps, ambiguous nucleotides or AAs, or to find all columns that
do not contain such things. Note that if sequences are of uneven lengths and
--any is not specified, only indices up to the length of the shortest input
sequence can be printed (i.e., there is a strict interpretation of all
sequences needing to have a matching base at an index: if the index does not
exist in even one sequence, then no base can occur in all sequences at that
index).
```

</details>

<details><summary><strong><code>fasta-identity-table.py</code></strong> — Print a FASTA sequence identity table (pairwise percent identity between
all sequences). <span title="Show --help">ℹ️</span></summary>

```
usage: fasta-identity-table.py [-h] [--text] [--strict] [--numberedColumns]
                               [--omitMatchDetails] [--digits DIGITS]
                               [--addZeroes] [--showLengths] [--showGaps]
                               [--showNoCoverage] [--showNs] [--footer]
                               [--div] [--fastaFile2 FILENAME]
                               [--gapChars CHARS] [--noCoverageChars CHARS]
                               [--defaultColor DEFAULTCOLOR] [--color COLOR]
                               [--align] [--upperOnly] [--concise] [--sort]
                               [--aligner {edlib,mafft,needle}]
                               [--emptyCellMinSide EMPTYCELLMINSIDE]
                               [--alignerOptions ALIGNEROPTIONS] [--verbose]
                               [--highlightBest] [--fastaFile FILENAME]
                               [--readClass CLASSNAME] [--fasta | --fastq |
                               --fasta-ss] [--minLength N] [--maxLength N]
                               [--maxNFraction N] [--whitelist SEQUENCE-ID]
                               [--blacklist SEQUENCE-ID]
                               [--whitelistFile SEQUENCE-ID-FILE]
                               [--blacklistFile SEQUENCE-ID-FILE]
                               [--titleRegex REGEX]
                               [--negativeTitleRegex REGEX]
                               [--sequenceWhitelist SEQUENCE]
                               [--sequenceBlacklist SEQUENCE]
                               [--sequenceWhitelistFile SEQUENCE-FILE]
                               [--sequenceBlacklistFile SEQUENCE-FILE]
                               [--sequenceRegex REGEX]
                               [--sequenceNegativeRegex REGEX]
                               [--keepSequences NUMBER,RANGE,... |
                               --removeSequences NUMBER,RANGE,...] [--head N]
                               [--removeDuplicates | --removeDuplicatesById]
                               [--removeDuplicatesUseMD5] [--randomSubset N]
                               [--trueLength N] [--sampleFraction FRACTION]
                               [--sequenceNumbersFile FILENAME]
                               [--keepSites KEEPSITES |
                               --keepSitesFile KEEPSITESFILE |
                               --removeSites REMOVESITES |
                               --removeSitesFile REMOVESITESFILE]
                               [--removeGaps] [--truncateTitlesAfter STR]
                               [--removeDescriptions]
                               [--idLambda LAMBDA-FUNCTION]
                               [--readLambda LAMBDA-FUNCTION] [--reverse]
                               [--reverseComplement] [--rotate N]
                               [--setSite site:base[:quality] [site:base[:quality] ...]]
                               [--setSiteZeroBased] [--upper | --lower]
                               [--upperId | --lowerId] [--setQuality Q]
                               [FILENAME]

Print a FASTA sequence identity table.

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given.

options:
  -h, --help            show this help message and exit
  --text                If specified, just print a simple text table.
  --strict              If given, do not allow ambiguous nucleotide symbols to
                        match.
  --numberedColumns     Use a sequence (row) number as the header of each
                        column instead of the sequence id.
  --omitMatchDetails    Do not show match details.
  --digits DIGITS       The number of digits to round identities to.
  --addZeroes           If given, make sure all identities have the same
                        number of digits shown by adding zeroes when needed.
  --showLengths         If given, show the lengths of sequences.
  --showGaps            If given, show the number of gaps in sequences.
  --showNoCoverage      If given, show the number of no-coverage characters in
                        sequences.
  --showNs              If given, show the number of fully ambiguous N
                        characters in sequences.
  --footer              If given, also show sequence ids at the bottom of the
                        table.
  --div                 If given, print an HTML <div> fragment only, not a
                        full HTML document (ignored if --text is used).
  --fastaFile2 FILENAME
                        The name of a second FASTA input file. If no second
                        FASTA file is given, sequences in the first FASTA file
                        will be compared with each other.
  --gapChars CHARS      The sequence characters that should be considered to
                        be gaps. These characters will be ignored in computing
                        sequence lengths and identity fractions.
  --noCoverageChars CHARS
                        The sequence characters that indicate lack of
                        coverage. These characters will be ignored in identity
                        fractions.
  --defaultColor DEFAULTCOLOR
                        The (background) color for cells. This will be used
                        for all cells that do not otherwise have a color due
                        to use of --color. This option is ignored if --text is
                        given.
  --color COLOR         Specify cell background coloring. This option must be
                        given as a space-separated "value color" pair. The
                        value is an identity fraction in [0..1] and the color
                        is any color specification that can be given to CSS.
                        This argument can be repeated. E.g., --color "0.9 red"
                        --color "0.75 rgb(23, 190, 207)" --color "0.1
                        #CF3CF3". Cells will be colored using the color of the
                        highest identity fraction they satisfy. The default is
                        to color all cells with the --defaultColor color. This
                        option is ignored if --text is given.
  --align               Do pairwise alignment of all sequences.
  --upperOnly           Only show the upper diagonal (only valid when a single
                        FASTA input file is given, resulting in a square all-
                        against-all identity table).
  --concise             When only one file is given (and the table is
                        therefore square with symmetric entries), Make the
                        table more concise by omitting the first column and
                        last row. These can be omitted because the first
                        column has the sequence that is in the first row
                        (where it is matched against everything else) and the
                        last row is the same sequence as the final column
                        (where it is matched against everything else). When
                        these two are omitted, the table is slightly smaller
                        and arguably slightly harder to understand.
  --sort                Sort the input sequences by id.
  --aligner {edlib,mafft,needle}
                        The alignment algorithm to use.
  --emptyCellMinSide EMPTYCELLMINSIDE
                        The CSS minimum width and height for empty HTML table
                        cells.
  --alignerOptions ALIGNEROPTIONS
                        Optional arguments to pass to the alignment algorithm.
                        If the aligner is mafft, the default options are '--
                        globalpair --maxiterate 1000 --preservecase'. If
                        needle, "auto". Do not try to set the number of
                        threads here - use the --threads argument instead. If
                        you are using mafft, see https://mafft.cbrc.jp/alignme
                        nt/software/algorithms/algorithms.html for some
                        possible option combinations.
  --verbose             Print progress to standard error.
  --highlightBest       Highlight the best identity in each row.
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead.
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ).
  --fastq               If specified, input will be treated as FASTQ.
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly.
  --minLength N         The minimum sequence length.
  --maxLength N         The maximum sequence length.
  --maxNFraction N      The maximum fraction of Ns that can be present in the
                        sequence.
  --whitelist SEQUENCE-ID
                        Sequence titles (ids) that should be whitelisted. May
                        be repeated.
  --blacklist SEQUENCE-ID
                        Sequence titles (ids) that should be blacklisted. May
                        be repeated.
  --whitelistFile SEQUENCE-ID-FILE
                        The name of a file that contains sequence titles (ids)
                        that should be whitelisted, one per line
  --blacklistFile SEQUENCE-ID-FILE
                        The name of a file that contains sequence titles (ids)
                        that should be blacklisted, one per line.
  --titleRegex REGEX    A regex that sequence titles (ids) must match.
  --negativeTitleRegex REGEX
                        A regex that sequence titles (ids) must not match.
  --sequenceWhitelist SEQUENCE
                        Sequences that should be whitelisted. May be repeated.
  --sequenceBlacklist SEQUENCE
                        Sequences that should be blacklisted. May be repeated.
  --sequenceWhitelistFile SEQUENCE-FILE
                        The name of a file that contains sequences that should
                        be whitelisted, one per line.
  --sequenceBlacklistFile SEQUENCE-FILE
                        The name of a file that contains sequences that should
                        be blacklisted, one per line.
  --sequenceRegex REGEX
                        A regex that sequences must match.
  --sequenceNegativeRegex REGEX
                        A regex that sequences must not match.
  --keepSequences NUMBER,RANGE,...
                        Specify (1-based) ranges of sequence numbers that
                        should be kept. E.g., --keepSequences 1-3,5 will
                        output just the 1st, 2nd, 3rd, and 5th sequences. All
                        others will be omitted.
  --removeSequences NUMBER,RANGE,...
                        Specify (1-based) ranges of sequence numbers that
                        should be removed. E.g., --removeSequences 1-3,5 will
                        output all but the 1st, 2nd, 3rd, and 5th sequences.
                        All others will be ouput.
  --head N              Only the first N sequences will be printed.
  --removeDuplicates    Duplicate reads will be removed, based only on
                        sequence identity. The first occurrence is kept.
  --removeDuplicatesById
                        Duplicate reads will be removed, based only on read
                        id. The first occurrence is kept.
  --removeDuplicatesUseMD5
                        MD5 sums will be stored instead of the full sequence
                        or read id when either --removeDuplicates or
                        removeDuplicatesById are given. One of
                        --removeDuplicates and --removeDuplicatesById must
                        also be given. Note that this makes duplicate removal
                        probabilistic. This option can be used to reduce the
                        amount of RAM consumed during duplicate removal.
  --randomSubset N      An integer giving the number of sequences that should
                        be kept. These will be selected at random.
  --trueLength N        The number of reads in the FASTA input. Only to be
                        used with randomSubset
  --sampleFraction FRACTION
                        A [0.0, 1.0] C{float} indicating a fraction of the
                        reads that should be allowed to pass through the
                        filter. The sample size will only be approximately the
                        product of the sample fraction and the number of
                        reads. The sample is taken at random.
  --sequenceNumbersFile FILENAME
                        A file of (1-based) sequence numbers to retain.
                        Numbers must be one per line.
  --keepSites, --keepIndices KEEPSITES
                        Specify 1-based sequence sites to keep. All other
                        sites will be removed. The sites must be given in the
                        form e.g., 24,100-200,260. Note that the requested
                        sites will be taken from the input sequences in order,
                        not in the order given by --keepSites. I.e.,
                        --keepSites 5,8-10 will get you the same result as
                        --keepSites 8-10,5.
  --keepSitesFile, --keepIndicesFile KEEPSITESFILE
                        Specify a file containing 1-based sites to keep. All
                        other sequence sites will be removed. Lines in the
                        file must be given in the form e.g., 24,100-200,260.
                        See --keepSites for more detail.
  --removeSites, --removeIndices REMOVESITES
                        Specify 1-based sites to remove. All other sequence
                        sites will be kept. The sites must be given in the
                        form e.g., 24,100-200,260. See --keepSites for more
                        detail.
  --removeSitesFile, --removeIndicesFile REMOVESITESFILE
                        Specify a file containing 1-based sites to remove. All
                        other sequence sites will be kept. Lines in the file
                        must be given in the form e.g., 24,100-200,260. See
                        --keepSites for more detail.
  --removeGaps          If True, gap ('-') characters in sequences will be
                        removed.
  --truncateTitlesAfter STR
                        A string that sequence titles (ids) will be truncated
                        beyond. If the truncated version of a title has
                        already been seen, that title will be skipped.
  --removeDescriptions  Read id descriptions will be removed. The description
                        is the part of a sequence id after the first
                        whitespace (if any).
  --idLambda LAMBDA-FUNCTION
                        A one-argument function taking and returning a read
                        id. E.g., --idLambda "lambda id: id.split('_')[0]" or
                        --idLambda "lambda id: id[:10]". If the function
                        returns None, the read will be filtered out.
  --readLambda LAMBDA-FUNCTION
                        A one-argument function taking and returning a read.
                        E.g., --readLambda "lambda r: Read(r.id.split('_')[0],
                        r.sequence.strip('-'))". Make sure to also modify the
                        quality string if you change the length of a FASTQ
                        sequence. If the function returns None, the read will
                        be filtered out. The function will be passed to eval
                        with the dark.reads classes Read, DNARead, AARead,
                        etc. all in scope.
  --reverse             Reverse the sequences. Note that this is NOT reverse
                        complementing.
  --reverseComplement   Reverse complement the sequences.
  --rotate N            Rotate sequences. A positive integer rotates the
                        sequence to the right, negative to the left
  --setSite site:base[:quality] [site:base[:quality] ...]
                        Set a (1-based) site to a certain base (with a
                        quality, if FASTQ). May be given with multiple
                        specifications to change several sites, e.g.:
                        --setSite 4:A 20:G Several things should be noted: 1)
                        because the sites are 1-based, you cannot use -1 as a
                        site to change the final site, you need to use 0
                        instead (since 1 is subtracted from the site indices
                        to get Python's 0-based offsets). If you find this
                        weird, use --setSiteZeroBased and give 0-based sites.
                        2) To pass a negative site you should give a leading
                        space so the argparse parser doesn't think your offset
                        is the start of another argument. Therefore: --setSite
                        ' -1:A' can be used. 3) These changes are applied
                        before any other read-altering operations (e.g.,
                        sequence rotation, applying a read lambda, reversing
                        or reverse complementing).
  --setSiteZeroBased    Indicate that the sites given in --setSite are zero
                        based.
  --upper               Convert sequences to uppercase.
  --lower               Convert sequences to lowercase.
  --upperId             Convert sequence IDs to uppercase.
  --lowerId             Convert sequence IDs to lowercase.
  --setQuality Q        Set all quality characters to 'Q'. Note that this is a
                        string value, not an integer. As well as changing
                        FASTQ qualities, this can also be used to convert
                        FASTA to FASTQ (with a constant quality for all bases)
                        when used with --saveAs fastq. A PHRED quality of 20
                        can be given as '5', of 30 as '?', and 32 as 'A'. See
                        man ascii for other values.
```

</details>

**`randomize-fasta.py`** — Read DNA FASTA from stdin and print FASTA to stdout, with the sequences
randomized (ids and lengths preserved).

<details><summary><strong><code>fastq-set-quality.py</code></strong> — Given FASTA/FASTQ on stdin, write FASTQ to stdout with all quality scores
set to a given value. <span title="Show --help">ℹ️</span></summary>

```
usage: fastq-set-quality.py [-h] [--quality QUALITY] [--fastaFile FILENAME]
                            [--readClass CLASSNAME] [--fasta | --fastq |
                            --fasta-ss]
                            [FILENAME]

Given FASTA/Q on stdin, write FASTQ to stdout with all nucleotide quality
scores set to a given value.

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given.

options:
  -h, --help            show this help message and exit
  --quality, -q QUALITY
                        The integer PHRED quality to set. Sanger (32-based
                        PHRED is used).
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead.
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ).
  --fastq               If specified, input will be treated as FASTQ.
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly.
```

</details>

<details><summary><strong><code>combine-sequences.py</code></strong> — Combine multiple sequences into a single sequence. <span title="Show --help">ℹ️</span></summary>

```
usage: combine-sequences.py [-h] --fastaFiles FASTAFILES [FASTAFILES ...]
                            [--threads THREADS]
                            [--alignerOptions ALIGNEROPTIONS]

Combine multiple sequences into a single sequence.

options:
  -h, --help            show this help message and exit
  --fastaFiles FASTAFILES [FASTAFILES ...]
                        The name of the FASTA files containing the sequences
                        to be combined.
  --threads THREADS     The number of threads to use when running the aligner
                        (if --align is used and the alignment algorithm can
                        make use of multiple threads). (default: 12)
  --alignerOptions ALIGNEROPTIONS
                        Optional arguments to pass to the alignment algorithm.
                        The default options are '--globalpair --maxiterate
                        1000 --preservecase'. See https://mafft.cbrc.jp/alignm
                        ent/software/algorithms/algorithms.html for some
                        possible option combinations. (default: --globalpair
                        --maxiterate 1000 --preservecase)
```

</details>

<details><summary><strong><code>summarize-fasta-bases.py</code></strong> — Given FASTA on stdin, write a summary of sequence base categories to
stdout (categories are currently hard-coded, not selectable on the
command line). <span title="Show --help">ℹ️</span></summary>

```
usage: summarize-fasta-bases.py [-h] [--baseType {nucl,prot}]
                                [--minLength MINLENGTH] [--concise]
                                [--fastaFile FILENAME] [--readClass CLASSNAME]
                                [--fasta | --fastq | --fasta-ss]
                                [FILENAME]

Given FASTA on stdin, write a summary of sequence base categories to stdout.
It is currently not possible to specify the categories on the command line.

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given.

options:
  -h, --help            show this help message and exit
  --baseType {nucl,prot}
                        The type of the bases in the input.
  --minLength MINLENGTH
                        If specified, stretches of reads that are less than
                        this length will not be reported but will be
                        summarized by an ellipsis.
  --concise             If specified, do not show the individual sequence
                        regions.
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead.
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ).
  --fastq               If specified, input will be treated as FASTQ.
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly.
```

</details>

<details><summary><strong><code>summarize-reads.py</code></strong> — Given a FASTA or FASTQ filename argument, print summary statistics: read
count, total/average/median/min/max length, and per-base counts. <span title="Show --help">ℹ️</span></summary>

```
usage: summarize-reads.py [-h] filename

Given a FASTA or FASTQ file, print summary statistics: read count,
total/average/median/min/max length, and per-base counts.

positional arguments:
  filename    The FASTA or FASTQ file to summarize.

options:
  -h, --help  show this help message and exit
```

</details>

<details><summary><strong><code>fasta-coverage.py</code></strong> — Given a set of (presumably same-length or aligned) FASTA sequences, print
per-character coverage/composition information (e.g., GC fraction) to
stdout. <span title="Show --help">ℹ️</span></summary>

```
usage: fasta-coverage.py [-h] [--digits DIGITS] [--noNames] [--strict] [--gc]
                         [--reverse] [--sortBy SORTBY] [--header] [--sep SEP]
                         [--printCounts] [--sortChars] [--fastaFile FILENAME]
                         [--readClass CLASSNAME] [--fasta | --fastq |
                         --fasta-ss]
                         [FILENAME]

Print coverage information to standard output.

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given.

options:
  -h, --help            show this help message and exit
  --digits DIGITS       The number of digits of precision to use for floating
                        point outputs (coverage, gc fraction).
  --noNames             Do not print the ids of input sequences.
  --strict              Treat all ambiguous nucleotides as not covered.
                        Additional characters specified with --chars will also
                        be counted as not covered.
  --gc                  Include the GC percentage in the output.
  --reverse, -r         Sort lines by decreasing value.
  --sortBy SORTBY       What to sort lines by. Either one of id, coverage,
                        length, gc, not-covered or else a single nucleotide
                        character.
  --header              Write a header line and in cells print just the
                        values.
  --sep SEP             The output separator character.
  --printCounts         Print the count for each character found in each read.
  --sortChars           Sort the output columns.
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead.
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ).
  --fastq               If specified, input will be treated as FASTQ.
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly.
```

</details>

<details><summary><strong><code>trim-reads.py</code></strong> — Given a FASTA filename argument, detect a prefix/suffix common to all
reads (e.g., leftover primer/adaptor) and print the reads with that
trimmed off. <span title="Show --help">ℹ️</span></summary>

```
usage: trim-reads.py [-h] filename

Given a FASTA file, detect a prefix/suffix common to all reads (e.g., leftover
primer or adaptor sequence) and print the reads with that trimmed off.

positional arguments:
  filename    The FASTA file to trim.

options:
  -h, --help  show this help message and exit
```

</details>

<details><summary><strong><code>trim-primers.py</code></strong> — Given FASTA on stdin, look for a primer sequence and write trimmed FASTA
(the part after the primer) to stdout. <span title="Show --help">ℹ️</span></summary>

```
usage: trim-primers.py [-h] [--verbose VERBOSE] primer

Given FASTA on stdin, look for a primer sequence and write trimmed FASTA
(after the primer) to stdout.

positional arguments:
  primer             The primer sequence.

options:
  -h, --help         show this help message and exit
  --verbose VERBOSE  Print information on found primers.
```

</details>

<details><summary><strong><code>split-fasta-by-adaptors.py</code></strong> — Given a set of known adaptors and FASTA on stdin, write separate FASTA
files split by adaptor. <span title="Show --help">ℹ️</span></summary>

```
usage: split-fasta-by-adaptors.py [-h] [--adaptor-offset N]
                                  [--maximum-distance dist]
                                  [--outputPrefix prefix] [--dry-run DRYRUN]
                                  [--verbose VERBOSE]
                                  adaptor [adaptor ...]

From a set of known adaptors and a FASTA on stdin, write FASTA files split by
adaptor

positional arguments:
  adaptor               the set of adaptors that were used in sequencing

options:
  -h, --help            show this help message and exit
  --adaptor-offset N    the offset of the adaptor in each sequence
  --maximum-distance dist
                        The maximum distance an unknown adaptor will be mapped
                        to in an attempt to find its nearest known adaptor.
  --outputPrefix prefix
                        the prefix to use for newly created FASTA files
  --dry-run DRYRUN      If True, do not write new FASTA files, just show what
                        would be done
  --verbose VERBOSE     If True, print information about adaptor classes
```

</details>

<details><summary><strong><code>adaptor-distances.py</code></strong> — Print a Levenshtein distance matrix for a set of known adaptors (helper
for choosing/checking adaptors used by the two scripts above). <span title="Show --help">ℹ️</span></summary>

```
usage: adaptor-distances.py [-h] adaptor [adaptor ...]

Print a (Levenshtein) distance matrix for a set of known adaptors

positional arguments:
  adaptor     The set of adaptors used in sequencing.

options:
  -h, --help  show this help message and exit
```

</details>

<details><summary><strong><code>reverse-complement.py</code></strong> — Print the reverse complement of one or more DNA sequences given directly
as command-line arguments (not stdin). <span title="Show --help">ℹ️</span></summary>

```
usage: reverse-complement.py [-h] SEQUENCE [SEQUENCE ...]

Print the reverse complement of one or more DNA sequences.

positional arguments:
  SEQUENCE    A DNA sequence.

options:
  -h, --help  show this help message and exit
```

</details>

<details><summary><strong><code>local-align.py</code></strong> — Given two FASTA files, perform a local Smith-Waterman alignment of all
sequences in one against all those in the other. <span title="Show --help">ℹ️</span></summary>

```
usage: local-align.py [-h] --fastaFile1 FASTAFILE1 --fastaFile2 FASTAFILE2
                      [--matchScore MATCHSCORE]
                      [--mismatchScore MISMATCHSCORE]
                      [--gapOpenScore GAPOPENSCORE]
                      [--gapExtendScore GAPEXTENDSCORE]
                      [--gapExtendDecay GAPEXTENDDECAY]

Given two FASTA files, perform a local Smith Waterman alignment on all the
sequences in one file against all those in the other.

options:
  -h, --help            show this help message and exit
  --fastaFile1 FASTAFILE1
                        The first FASTA file.
  --fastaFile2 FASTAFILE2
                        The second FASTA file.
  --matchScore MATCHSCORE
                        The match score. (default: 1)
  --mismatchScore MISMATCHSCORE
                        The mismatch score. (default: -1)
  --gapOpenScore GAPOPENSCORE
                        The gap open score. (default: -1)
  --gapExtendScore GAPEXTENDSCORE
                        The gap extend score. (default: -1)
  --gapExtendDecay GAPEXTENDDECAY
                        The gap extend decay. (default: 0.0)
```

</details>

<details><summary><strong><code>fasta-to-phylip.py</code></strong> — Given FASTA sequences of equal length on stdin, write Phylip format to
stdout (for feeding into Phylip-family or other tree-building tools —
could equally live in the "phylogenetic trees" category below). <span title="Show --help">ℹ️</span></summary>

```
usage: fasta-to-phylip.py [-h] [--addIdSpace] [--fastaFile FILENAME]
                          [--readClass CLASSNAME] [--fasta | --fastq |
                          --fasta-ss]
                          [FILENAME]

Given FASTA sequences (of equal length) on stdin write Phylip to stdout.

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given.

options:
  -h, --help            show this help message and exit
  --addIdSpace          If True, print an extra space after each sequence id
                        (the space is expected by some programs that process
                        Phylip files, e.g., baseml).
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead.
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ).
  --fastq               If specified, input will be treated as FASTQ.
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly.
```

</details>

## 2. DNA/amino-acid translation and codon tools

<details><summary><strong><code>dna-to-aa.py</code></strong> — Convert DNA to amino acids. <span title="Show --help">ℹ️</span></summary>

```
usage: dna-to-aa.py [-h] [--minORFLength LEN] [--fastaFile FILENAME]
                    [--readClass CLASSNAME] [--fasta | --fastq | --fasta-ss]
                    [FILENAME]

Convert DNA to AA

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given.

options:
  -h, --help            show this help message and exit
  --minORFLength LEN    Translations to AA that do not contain an ORF of at
                        least this length will not be produced.
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead.
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ).
  --fastq               If specified, input will be treated as FASTQ.
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly.

Given DNA FASTA on stdin, output AA FASTA to stdout. Optionally, filter by
minimum required ORF length.
```

</details>

<details><summary><strong><code>aa-to-dna.py</code></strong> — Given amino acid FASTA on stdin, write FASTA with DNA sequences that
would translate to the given amino acid sequences (codons chosen fixed or
at random). <span title="Show --help">ℹ️</span></summary>

```
usage: aa-to-dna.py [-h] [--random] [--fastaFile FILENAME]
                    [--readClass CLASSNAME] [--fasta | --fastq | --fasta-ss]
                    [FILENAME]

Given amino acid FASTA on stdin, write FASTA with DNA sequences to stdout,
where the DNA sequences would translate to the given amino acid seqeunces. The
DNA bases used for each amino acid is either fixed or is chosen at random from
the corresponding codons.

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given.

options:
  -h, --help            show this help message and exit
  --random              Use codons chosen at random for each amino acid.
                        Otherwise a fixed codon will be used for each amino
                        acid so that the output for a given amino acid
                        sequence is always identical.
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead.
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ).
  --fastq               If specified, input will be treated as FASTQ.
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly.
```

</details>

<details><summary><strong><code>fasta-translate.py</code></strong> — Translate FASTA and write the result to standard output. <span title="Show --help">ℹ️</span></summary>

```
usage: fasta-translate.py [-h] [--fastaFile FILENAME] [--readClass CLASSNAME]
                          [--fasta | --fastq | --fasta-ss]
                          [FILENAME]

Translate FASTA and write the result to standard output.

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given.

options:
  -h, --help            show this help message and exit
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead.
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ).
  --fastq               If specified, input will be treated as FASTQ.
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly.
```

</details>

<details><summary><strong><code>translate.py</code></strong> — Print the amino acid translation of one or more DNA sequences given as
command-line arguments. <span title="Show --help">ℹ️</span></summary>

```
usage: translate.py [-h] DNA-SEQUENCE [DNA-SEQUENCE ...]

Print the amino acid translation of one or more DNA sequences.

positional arguments:
  DNA-SEQUENCE  A DNA sequence to translate.

options:
  -h, --help    show this help message and exit
```

</details>

<details><summary><strong><code>extract-ORFs.py</code></strong> — Convert DNA or amino-acid FASTA to amino-acid, ORF-only FASTA. <span title="Show --help">ℹ️</span></summary>

```
usage: extract-ORFs.py [-h] [--allowOpenORFs] [--minORFLength LEN]
                       [--maxORFLength LEN] [--kozakOnly]
                       [--kozakInfoFile KOZAKINFOFILE] [--fastaFile FILENAME]
                       [--readClass CLASSNAME] [--fasta | --fastq |
                       --fasta-ss]
                       [FILENAME]

Convert DNA or AA FASTA to AA ORF-only FASTA

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given.

options:
  -h, --help            show this help message and exit
  --allowOpenORFs       If True, ORFs that do not meet the length requirement
                        (as specified by --minORFLength) will be output as
                        long as they are open.
  --minORFLength LEN    Only ORFs of at least this length will be written to
                        stdout.
  --maxORFLength LEN    Only ORFs of a maximum of this length will be written
                        to stdout.
  --kozakOnly           Only ORFs that also have a Kozak consensus will be
                        written to stdout. Only applicable if DNA reads are
                        given.
  --kozakInfoFile KOZAKINFOFILE
                        Filename to which all Kozak consensus information is
                        written. Only applicable if DNA reads are given.
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead.
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ).
  --fastq               If specified, input will be treated as FASTQ.
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly.

Given DNA or AA FASTA on stdin, output AA ORF-only FASTA to stdout.
Optionally, filter by minimum required ORF length and allow the output of
short ORFs that are open.
```

</details>

<details><summary><strong><code>aa-info.py</code></strong> — Print detailed information (name, 1- and 3-letter abbreviations, codons,
properties) for one or more amino acids or codons named on the command
line, or for all amino acids if none are given. <span title="Show --help">ℹ️</span></summary>

```
usage: aa-info.py [-h] [--details] [targets ...]

Print detailed information (name, 1- and 3-letter abbreviations, codons,
properties) for one or more amino acids or codons named on the command line,
or for all amino acids if none are given.

positional arguments:
  targets     The items to look up.

options:
  -h, --help  show this help message and exit
  --details   Display numeric amino acid property values.
```

</details>

<details><summary><strong><code>aa-compare.py</code></strong> — Compare two amino acids (or codons), given as command-line arguments:
look each one up by name/abbreviation/codon, then print their shared and
different properties and the possible codon-to-codon paths (with
transition counts) between them, with colored terminal output. <span title="Show --help">ℹ️</span></summary>

```
usage: aa-compare.py [-h] AMINO_ACID_OR_CODON_1 AMINO_ACID_OR_CODON_2

Compare two amino acids or codons: look each one up by name, 3- or 1-letter
abbreviation, or codon, then print the properties they have in common and
those that differ, their codons, and the possible codon-to-codon paths (with
transition/transversion counts) between them.

positional arguments:
  AMINO_ACID_OR_CODON_1
                        The first amino acid or codon.
  AMINO_ACID_OR_CODON_2
                        The second amino acid or codon.

options:
  -h, --help            show this help message and exit
```

</details>

## 3. Pairwise sequence comparison

<details><summary><strong><code>compare-sequences.py</code></strong> — Compare two DNA sequences: align them (MAFFT by default, or another
aligner), then print a detailed match report (exact/ambiguous matches,
mismatches, gaps, etc., via <code>matchToString</code>). <span title="Show --help">ℹ️</span></summary>

```
usage: compare-sequences.py [-h] [--index1 INDEX1] [--index2 INDEX2] [--align]
                            [--aligner {edlib,mafft,needle}]
                            [--alignerOptions ALIGNEROPTIONS]
                            [--threads THREADS]
                            [--alignmentFile ALIGNMENTFILE] [--strict]
                            [--quiet] [--noGapLocations]
                            [--noNoCoverageLocations] [--sites SITES]
                            [--showDiffs] [--showAmbiguous]
                            [--includeAmbiguousMatches]
                            [--includeNonGapMismatches]
                            [--includeDifferenceCounts]
                            [--includeDifferenceLocations] [--gapChars CHARS]
                            [--noCoverageChars CHARS] [--fastaFile FILENAME]
                            [--readClass CLASSNAME] [--fasta | --fastq |
                            --fasta-ss]
                            [FILENAME]

Compare two DNA sequences.

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given. (default: None)

options:
  -h, --help            show this help message and exit
  --index1 INDEX1       The (1-based) index in the input of the first
                        sequence. (default: 1)
  --index2 INDEX2       The (1-based) index in the input of the second
                        sequence. (default: 2)
  --align               If given, use MAFFT (the default), edlib, or needle
                        (according to the algorithm selected by --aligner) to
                        align the two sequences. (default: False)
  --aligner {edlib,mafft,needle}
                        The alignment algorithm to use when --align is given.
                        Default 'mafft'. Implies --align. (default: None)
  --alignerOptions ALIGNEROPTIONS
                        Optional arguments to pass to the alignment algorithm.
                        If the aligner is 'mafft', the default options are '--
                        auto --maxiterate 1000 --preservecase'. If 'needle',
                        the default is 'auto'. Do not try to set the number of
                        threads here; use the --threads argument instead. If
                        you are using mafft, see https://mafft.cbrc.jp/alignme
                        nt/software/algorithms/algorithms.html for some
                        possible option combinations. Implies --align.
                        (default: None)
  --threads THREADS     The number of threads to use when running the aligner
                        (if --align is used and the alignment algorithm can
                        make use of multiple threads (mafft can, needle
                        cannot)). (default: 12)
  --alignmentFile ALIGNMENTFILE
                        The file to save the alignment to (implies --align).
                        (default: None)
  --strict              If given, do not allow ambiguous nucleotide symbols to
                        match. (default: False)
  --quiet               Do not print information about aligning. (default:
                        True)
  --noGapLocations      Do not indicate the (1-based) locations of sequence
                        gaps. (default: True)
  --noNoCoverageLocations
                        Do not indicate the (1-based) locations of no
                        coverage. (default: True)
  --sites SITES         Specify (1-based) sequence sites to keep. All other
                        sites will be ignored. The sites must be given in the
                        form e.g., 24,100-200,260. (default: None)
  --showDiffs           Print (1-based) sites where the sequence nucleotides
                        differ. (default: False)
  --showAmbiguous       Print (1-based) sites where either sequence has an
                        ambiguous nucleotide code. (default: False)
  --includeAmbiguousMatches
                        Print (1-based) sites of ambiguous matches. The output
                        gives the site, the base from the first sequence, then
                        the base from the second sequence. (default: False)
  --includeNonGapMismatches
                        Print (1-based) sites of mismatches that do not
                        involve a gap. The output gives the site, the base
                        from the first sequence, then the base from the second
                        sequence. (default: False)
  --includeDifferenceCounts
                        Include the counts (i.e., number of occurrences) of
                        all differences. Note that mismatches involving areas
                        of no coverage or where one read is shorter than
                        another do not count as differences. (default: False)
  --includeDifferenceLocations
                        Include the locations (and counts) (i.e., number of
                        occurrences) of all differences. Note that mismatches
                        involving areas of no coverage or where one read is
                        shorter than another do not count as differences.
                        (default: False)
  --gapChars CHARS      The sequence characters that should be considered to
                        be gaps. These characters will be ignored in computing
                        sequence lengths and identity fractions. (default: -)
  --noCoverageChars CHARS
                        The sequence characters that indicate lack of
                        coverage. These characters will be ignored in identity
                        fractions. (default: None)
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
                        (default: None)
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead. (default: DNARead)
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ). (default: False)
  --fastq               If specified, input will be treated as FASTQ.
                        (default: False)
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly. (default: False)
```

</details>

<details><summary><strong><code>compare-aa-sequences.py</code></strong> — The amino-acid equivalent of <code>compare-sequences.py</code>, aligning with
<code>needle</code> (EMBOSS) and reporting the comparison via <code>dark.aa.matchToString</code>. <span title="Show --help">ℹ️</span></summary>

```
usage: compare-aa-sequences.py [-h] [--index1 INDEX1] [--index2 INDEX2]
                               [--align] [--alignmentFile ALIGNMENTFILE]
                               [--sites SITES] [--showDiffs]
                               [--fastaFile FILENAME] [--readClass CLASSNAME]
                               [--fasta | --fastq | --fasta-ss]
                               [FILENAME]

Compare two amino acid sequences.

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given. (default: None)

options:
  -h, --help            show this help message and exit
  --index1 INDEX1       The (1-based) index in the input of the first
                        sequence. (default: 1)
  --index2 INDEX2       The (1-based) index in the input of the second
                        sequence. (default: 2)
  --align               If given, use needle to do a global alignment of the
                        two sequences. (default: False)
  --alignmentFile ALIGNMENTFILE
                        The file to save the alignment to (implies --align).
                        (default: None)
  --sites SITES         Specify (1-based) sequence sites to keep. All other
                        sites will be ignored. The sites must be given in the
                        form e.g., 24,100-200,260. (default: None)
  --showDiffs           Print (1-based) sites where the sequence amino acids
                        differ. (default: False)
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
                        (default: None)
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead. (default: DNARead)
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ). (default: False)
  --fastq               If specified, input will be treated as FASTQ.
                        (default: False)
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly. (default: False)
```

</details>

## 4. Multiple-sequence-alignment (MSA) and alignment-position analysis

<details><summary><strong><code>msa-find-and-extract.py</code></strong> — Find nucleotide sequences in an MSA and extract the matched region. <span title="Show --help">ℹ️</span></summary>

```
usage: msa-find-and-extract.py [-h] [--fastaFile FILENAME]
                               [--readClass CLASSNAME] [--fasta | --fastq |
                               --fasta-ss] [--prefix PREFIX] [--suffix SUFFIX]
                               [--id ID] [--allowUnequalLengths]
                               [--allowDistributedMatch] [--adjustIds]
                               [--verbose] [--zeroBased]
                               [--gapCharacter GAPCHARACTER] [--tsv TSV]
                               [FILENAME]

Find nucleotide sequences in an MSA and extract the matched region.

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given.

options:
  -h, --help            show this help message and exit
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead.
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ).
  --fastq               If specified, input will be treated as FASTQ.
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly.
  --prefix PREFIX       The start of the nucleotide sequence.
  --suffix SUFFIX       The end of the nucleotide sequence.
  --id ID               The id of the input sequence in which the prefix and
                        suffix must be found. If not specified, all input
                        sequences will be examined and those in which the
                        prefix/suffix are found must all agree.
  --allowUnequalLengths
                        Do not exit if all sequences are not the same length.
                        I.e., if this is not an MSA.
  --allowDistributedMatch
                        If --prefix and --suffix are both used, the default
                        behaviour is to make sure that they are found together
                        in at least one of the input sequences. If you don't
                        mind if the prefix is found in one sequence and the
                        suffix in another, use this option to allow it.
  --adjustIds           Add ' (position X-Y)' to the end of the output
                        sequence ids to indicate the subsequence range of the
                        result.
  --verbose             Print (to standard error) a summary of total prefixes
                        and suffixes found.
  --zeroBased           If --tsv or --verbose are used, the reported prefix
                        and suffix positions will refer to zero-based Python
                        style [closed, open) regions.
  --gapCharacter GAPCHARACTER
                        The 1-letter gap character.
  --tsv TSV             The (optional) name of a TSV file to write match
                        offsets to.
```

</details>

<details><summary><strong><code>msa-merge-gapped-region.py</code></strong> — Find nucleotide sequences bounding a region in an MSA and merge the gaps
within that region into one gap, placed at a specified offset, moving the
non-gap characters to one side (or around it). <span title="Show --help">ℹ️</span></summary>

```
usage: msa-merge-gapped-region.py [-h] [--fastaFile FILENAME]
                                  [--readClass CLASSNAME] [--fasta | --fastq |
                                  --fasta-ss] [--prefix PREFIX]
                                  [--suffix SUFFIX] [--id ID]
                                  [--gapOffset GAPOFFSET]
                                  [--allowDistributedMatch] [--verbose]
                                  [--zeroBased] [--gapCharacter GAPCHARACTER]
                                  [--tsv TSV]
                                  [FILENAME]

Find nucleotide sequences bounding a region in an MSA and merge the gaps
within that region. I.e., merge all the gaps into a single gap, place it at a
specific offset within the region, and put the non-gap characters (now
contiguous) to the left or right (or around) the gap.

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given.

options:
  -h, --help            show this help message and exit
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead.
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ).
  --fastq               If specified, input will be treated as FASTQ.
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly.
  --prefix PREFIX       The start of the match to find in sequences.
  --suffix SUFFIX       The end of the match to find in sequences.
  --id ID               The id of the input sequence in which the prefix and
                        suffix must be found. If not specified, all input
                        sequences will be examined and the offsets in which
                        the prefix/suffix are found must all agree.
  --gapOffset GAPOFFSET
                        The offset at which the merged gaps should be put. The
                        default (0) will put the merged gaps on the left
                        (i.e., with zero non-gap chars to the left of them.
                        This value cannot be larger than the number of non-gap
                        chars in the sequence region (otherwise there are not
                        enough non-gap chars to put before the gap). Use -1 to
                        put the merged gaps on the right.
  --allowDistributedMatch
                        If --prefix and --suffix are both used, the default
                        behaviour is to make sure that they are both found in
                        at least one of the input sequences. If you don't mind
                        if the prefix is found in one sequence and the suffix
                        in another, use this option to allow it.
  --verbose             Print a summary of prefixes/suffixes found to standard
                        error.
  --zeroBased           If --tsv or --verbose are used, the reported prefix
                        and suffix positions will refer to zero-based Python
                        style [closed, open) regions.
  --gapCharacter GAPCHARACTER
                        The 1-letter alignment gap character.
  --tsv TSV             The (optional) name of a TSV file to write match
                        offsets to.
```

</details>

<details><summary><strong><code>window-split-alignment.py</code></strong> — Split a FASTA alignment into equal-sized windows. <span title="Show --help">ℹ️</span></summary>

```
usage: window-split-alignment.py [-h] [--fastaFile FILENAME]
                                 [--readClass CLASSNAME] [--fasta | --fastq |
                                 --fasta-ss] [--base BASE] [--window WINDOW]
                                 [--verbose] [--zeroBased]
                                 [--minWindow MINWINDOW]
                                 [--maxWindows MAXWINDOWS] [--rotate N]
                                 [--startOffset STARTOFFSET]
                                 [--invalidChars INVALIDCHARS]
                                 [FILENAME]

Split a FASTA alignment into equal-sized windows.

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given. (default: None)

options:
  -h, --help            show this help message and exit
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
                        (default: None)
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead. (default: DNARead)
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ). (default: False)
  --fastq               If specified, input will be treated as FASTQ.
                        (default: False)
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly. (default: False)
  --base BASE           The base output filename. If you want output to be
                        placed in a sub-directory, specify the sub-directory
                        in the value. Required non-existent intermediate
                        directories (if any) will be created. (default:
                        window-)
  --window WINDOW       The window size. (default: 100)
  --verbose             Write processing information to standard error.
                        (default: False)
  --zeroBased           The --startOffset and offsets used in filenames will
                        be zero-based. If not given, offsets are one-based.
                        (default: True)
  --minWindow MINWINDOW
                        If the length of the final windowed sequence is less
                        than this value, do not write a file for it. (default:
                        None)
  --maxWindows MAXWINDOWS
                        The maximum number of windowed subsequences to
                        extract. (default: None)
  --rotate N            Rotate sequences left (if N is negative) or right (if
                        positive) before extracting sub-sequence windows.
                        Output filenames will be adjusted to match the
                        rotation amount and the sequence length. For example,
                        if the window size is 200, a --rotate value of -100 is
                        given, and the sequences have length 2000, the first
                        output file will have a name like
                        window-1901-0100.fasta (or window-1900-0100.fasta if
                        --zeroBased is used) and the final name will be
                        window-1701-1900.fasta. This option is useful for
                        making windows that match the location of genome
                        features, especially in the case of a circular genome
                        (e.g., for hepatitis B virus). (default: 0)
  --startOffset STARTOFFSET
                        The zero-based genome offset to start windowing from
                        (default value 0). Use --oneBased if you would prefer
                        to use a one-based offset (in which case the default
                        value is 1 and indices in the output files will also
                        be one-based). (default: None)
  --invalidChars INVALIDCHARS
                        Sequence characters (case insensitive) to consider
                        invalid for the purposes of tree-making. If a windowed
                        sequence for a sample consists entirely of these
                        characters, no sequence will be written to the output
                        for that window of that sample. (default: N?-)
```

</details>

<details><summary><strong><code>plot-windowed-identity.py</code></strong> — Make a plot showing windowed nucleotide identity across an alignment. <span title="Show --help">ℹ️</span></summary>

```
usage: plot-windowed-identity.py [-h] [--fastaFile FILENAME]
                                 [--readClass CLASSNAME] [--fasta | --fastq |
                                 --fasta-ss] --reference REFERENCE
                                 [--include INCLUDE] [--window WINDOW]
                                 [--minWindow MINWINDOW] [--jump JUMP]
                                 [--sort] [--startOffset STARTOFFSET]
                                 [--stopOffset STOPOFFSET] [--html HTML]
                                 [--out OUT] [--title TITLE]
                                 [--legendTitle LEGENDTITLE] [--dryRun]
                                 [--bestOnly] [--matplotlib] [--y01]
                                 [--ambiguous] [--markerSize MARKERSIZE]
                                 [FILENAME]

Make a plot showing windowed nucleotide identity across an alignment.

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given. (default: None)

options:
  -h, --help            show this help message and exit
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
                        (default: None)
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead. (default: DNARead)
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ). (default: False)
  --fastq               If specified, input will be treated as FASTQ.
                        (default: False)
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly. (default: False)
  --reference REFERENCE
                        A regular expression indicating the reference sequence
                        id (i.e., the sequence that others should be compared
                        to. This must obviously match exactly one input
                        sequence.
  --include INCLUDE     A regular expression indicating the sequence ids to
                        include in the plot. If not specified, all sequences
                        (other than the reference itself) will be compared to
                        the reference. (default: .)
  --window WINDOW       The width (number of nucleotides) of the window.
                        (default: 200)
  --minWindow MINWINDOW
                        The smallest window size to examine. (default: 50)
  --jump JUMP           The number of nucleotides to move across to the start
                        of the next window. (default: 50)
  --sort                Sort the ids of the sequences that are compared. If
                        not given, the order that matching ids are found in
                        the input FASTA will be used. (default: False)
  --startOffset STARTOFFSET
                        The zero-based genome offset to start the first window
                        at. (default: 0)
  --stopOffset STOPOFFSET
                        The zero-based genome offset to stop the last window
                        at. (default: None)
  --html HTML           The name of an HTML file to write a plotly image to.
                        (default: None)
  --out OUT             The name of the file to write the image to. The file
                        suffix will determine the output format. The image
                        will be created with plotly unless --matplotlib is
                        specified. (default: None)
  --title TITLE         The plot title. (default: None)
  --legendTitle LEGENDTITLE
                        The legend title. (default: Sample)
  --dryRun, -n          Just print what would be compared without actually
                        doing it. (default: False)
  --bestOnly            Only show the best matching sequences for each window.
                        (default: False)
  --matplotlib          Use matplotlib to make a static image. Otherwise,
                        plotly will be used. (default: False)
  --y01                 Set the y-axis limits to (0.0, 1.0). Only applicable
                        for matplotlib plots. (default: False)
  --ambiguous           Count compatible ambiguous nucleotide codes as
                        matching when computing identity. (default: True)
  --markerSize MARKERSIZE
                        The size of the markers in plot lines. (default: None)
```

</details>

<details><summary><strong><code>fasta-variable-sites.py</code></strong> — Given FASTA on stdin, write reads with only the variable sites to
standard output (nothing is written if no sites are variable). <span title="Show --help">ℹ️</span></summary>

```
usage: fasta-variable-sites.py [-h] [--printSites] [--sitesOnly]
                               [--homogeneous HOMOGENEOUS] [--confirm]
                               [--reference REFERENCE] [--unknownAreAmbiguous]
                               [--fastaFile FILENAME] [--readClass CLASSNAME]
                               [--fasta | --fastq | --fasta-ss]
                               [FILENAME]

Given FASTA on standard input, write reads with only the variable sites to
standard output. If no sites are variable, nothing is written to standard
output.

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given.

options:
  -h, --help            show this help message and exit
  --printSites          Print the variable site locations to standard error
                        along with the nucleotides found at the site, the
                        homogeneity fraction of the site and the reference
                        base (if a reference is given).
  --sitesOnly           Only print the (comma separated) variable sites to
                        standard output. This output can then be used directly
                        as an argument for --keepSites to filter-fasta.py
  --homogeneous HOMOGENEOUS
                        If the most-common nucleotide frequency at a site is
                        at least this value, the site will be considered
                        homogeneous.
  --confirm             Only keep sites where there is confirm variation
                        (i.e., ambiguous sites that are compatible with there
                        being no variation are not included).
  --reference REFERENCE
                        A sequence id. Site offsets will be reported relative
                        to this sequence. That means the number of gap
                        (hyphen) chars at the start of this sequence will be
                        counted and offsets will be incremented by that amount
                        (and lower offsets in sequences that start before the
                        specified sequence will be ignored).
  --unknownAreAmbiguous
                        Any unknown character (e.g., a '-' gap or '?' unknown
                        base) will be treated as being fully ambiguous (i.e.,
                        could be any of ACGT). Otherwise, all unknown
                        characters are counted as '-' characters.
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead.
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ).
  --fastq               If specified, input will be treated as FASTQ.
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly.
```

</details>

<details><summary><strong><code>position-summary.py</code></strong> — For a FASTA file of (presumably aligned) sequences, summarize what's
happening at one specific (1-based) position. <span title="Show --help">ℹ️</span></summary>

```
usage: position-summary.py [-h] --fastaFile FASTAFILE --position POSITION

For a fasta file with sequences, summarize what is happening at a specific
position.

options:
  -h, --help            show this help message and exit
  --fastaFile FASTAFILE
                        The name of the FASTA file to read.
  --position POSITION   The (one-based) position to examine.
```

</details>

## 5. Consensus-sequence building

<details><summary><strong><code>make-consensus.py</code></strong> — Make a consensus sequence from a FASTA reference and a BAM file, using external tools. <span title="Show --help">ℹ️</span></summary>

```
usage: make-consensus.py [-h] --reference REFERENCE [--bam BAM]
                         [--vcfFile VCFFILE] [--id ID |
                         --idLambda LAMBDA-FUNCTION] [--sample SAMPLE]
                         [--dryRun] [--maxDepth MAXDEPTH] [--ploidy PLOIDY]
                         [--maskLowCoverage MASKLOWCOVERAGE] [--log]
                         [--noClean] [--callHaplotypesGATK]
                         [--picardJar PICARDJAR] [--ivar] [--quiet]
                         [--ivarFrequencyThreshold IVARFREQUENCYTHRESHOLD]
                         [--ivarBedFile IVARBEDFILE]

Make a consensus sequence.

options:
  -h, --help            show this help message and exit
  --reference REFERENCE
                        The reference FASTA file.
  --bam BAM             The BAM file from which the consensus should be made.
                        Required if --maskLowCoverage is used. If no BAM file
                        is given, a VCF file must be provided. If both a BAM
                        and a VCF file are given, the VCF file will take
                        precedence. (default: None)
  --vcfFile VCFFILE     The VCF file. If omitted, bcftools will be used to
                        make a VCF file from the BAM file. (default: None)
  --id ID               The id to use in the consensus sequence in the output
                        FASTA. If not given, the reference sequence id will be
                        used. (default: None)
  --idLambda LAMBDA-FUNCTION
                        A one-argument function taking and returning a read
                        id. This can be used to set the id of the consensus
                        sequence based on the id of the reference sequence
                        (the function will be called with the id of the
                        reference sequence). E.g., --idLambda "lambda id:
                        id.split('_')[0]" or --idLambda "lambda id: id[:10] +
                        '-consensus'". (default: None)
  --sample SAMPLE       The name of the sample (from the @RG SM tag in the
                        original alignment BAM file) for which a consensus
                        should be made. If not given, the first sample name
                        (from the #CHROM header) in the VCF file will be used.
                        (default: None)
  --dryRun              Do not run commands, just print what would be done.
                        (default: False)
  --maxDepth MAXDEPTH   The maximum read depth to consider when calling the
                        consensus via bcftools consensus. (default: 5000)
  --ploidy PLOIDY       The ploidy to use when calling the consensus via
                        bcftools consensus. Use 1 for viruses, but be aware
                        that bcftools will NOT emit ambiguous nucleotide
                        characters in that case. (default: 1)
  --maskLowCoverage, --minimumCoverageDepth, --minCoverageDepth MASKLOWCOVERAGE
                        Put an N into sites where the coverage is below the
                        specified cutoff. If you specify a negative numer,
                        masking will be turned off. Note that specifying 0 is
                        likely to cause problems, either with the bcftools or
                        the samtools mpileup | ivar calling. If you're using
                        bcftools and you set --maskLowCoverage to 0, bcftools
                        will put the consensus base into all non-covered sites
                        (because no low-coverage sites will be masked out by
                        bcftools filter). If you're using ivar, the piped
                        output from samtools mpileup -aa will include sites
                        with zero coverage and the -m 0 given to ivar will
                        result in ivar emitting nothing at all for those sites
                        (so your consensus will end up being shorter than your
                        reference, which may well not be what you want (i.e.,
                        might be totally wrong)). Requires --bam. (default: 1)
  --log                 Show a log of commands that were (or would be, if
                        --dryRun is used) executed. (default: False)
  --noClean             Do not remove intermediate files or the temporary
                        directory. (default: True)
  --callHaplotypesGATK  Use GATK to call haplotypes. See
                        https://gatk.broadinstitute.org for details on GATK.
                        (default: False)
  --picardJar PICARDJAR
                        The path to the Picard jar file. See
                        https://github.com/broadinstitute/picard for details
                        on Picard. (default: None)
  --ivar, --iVar        If given, iVar will be used to call the consensus.
                        (default: False)
  --quiet               Do not emit a warning about --maskLowCoverage being
                        set to zero. (default: False)
  --ivarFrequencyThreshold IVARFREQUENCYTHRESHOLD
                        The frequency threshold used by iVar when calling the
                        consensus. If the frequency of the most-common
                        nucleotide at a site meets this threshold, the
                        nucleotide will be called. Otherwise, an ambiguous
                        nucleotide code will be produced, based on the
                        smallest set of most-frequent nucleotides whose summed
                        frequencies meet the threshold. See https://andersen-
                        lab.github.io/ivar/html/manualpage.html#autotoc_md19
                        for more information. If not given, 0.6 is used. Only
                        used if --ivar is specified. (default: None)
  --ivarBedFile IVARBEDFILE
                        If iVar should trim primers, a BED file of the primer
                        positions. (default: None)
```

</details>

## 6. Read mapping and SAM/BAM analysis

<details><summary><strong><code>run-bwa.py</code></strong> — Run <code>bwa</code> on a FASTA file; optionally convert the result to BAM, sort,
and index it. <span title="Show --help">ℹ️</span></summary>

```
usage: run-bwa.py [-h] --bwaIndex BWAINDEX [--bwaArgs BWAARGS]
                  --fastaFile FASTAFILE [--base BASE] [--verbose] [--noBAM]
                  [--noSort] [--noIndex] [--noClean] [--force] [--dryRun]

Run bwa on a FASTA file. Optionally convert the result to BAM, sorting, and
indexing.

options:
  -h, --help            show this help message and exit
  --bwaIndex BWAINDEX   The name of the BWA index (created with bwa index).
  --bwaArgs BWAARGS     All arguments to be passed to BWA. (default: mem)
  --fastaFile FASTAFILE
                        The name of the FASTA file whose sequences should be
                        aligned with the BWA index.
  --base BASE           The base name of files to create. Suffixes such as
                        .sam and .bam will be added. If not given, the
                        basename of the --fastaFile name will be used,
                        stripped of its final suffix. (default: None)
  --verbose             If specified, show the commands that were (or would
                        be) executed. (default: False)
  --noBAM               If specified, do not convert SAM to BAM. (default:
                        False)
  --noSort              If specified, do not sort the BAM. (default: False)
  --noIndex             If specified, do not index the sorted BAM. (default:
                        False)
  --noClean             If specified, do not remove intermediate .sam and .bam
                        files. (default: False)
  --force               If specified, overwrite pre-existing output files.
                        (default: False)
  --dryRun              If specified, do not run commands, just print what
                        would be done. (default: False)
```

</details>

<details><summary><strong><code>run-bowtie2.py</code></strong> — Run <code>bowtie2</code> on a FASTA file; optionally convert the result to BAM,
sort, and index it. <span title="Show --help">ℹ️</span></summary>

```
usage: run-bowtie2.py [-h] [--index INDEX] [--ignoreIndex IGNOREDINDICES]
                      [--fastq1 FASTQ1] [--fastq2 FASTQ2]
                      [--bowtie2Args BOWTIE2ARGS]
                      [--samtoolsViewArgs SAMTOOLSVIEWARGS]
                      [--tempdir TEMPDIR] [--out OUT] [--reference REFERENCE]
                      [--vcfFile VCFFILE] [--markDuplicatesGATK]
                      [--markDuplicatesPicard] [--picardJar PICARDJAR]
                      [--tmpChmod TMPCHMOD] [--removeDuplicates] [--verbose]
                      [--log] [--threads THREADS] [--noAlign] [--noBAM]
                      [--noSort] [--noIndexBAM] [--noClean] [--force]
                      [--dryRun] [--callHaplotypesGATK |
                      --callHaplotypesBcftools]
                      [--removePrimersFromBedFile REMOVEPRIMERSFROMBEDFILE]

Run bowtie2 on a FASTA file. Optionally convert the result to BAM, sorting,
and indexing.

options:
  -h, --help            show this help message and exit
  --index INDEX         Either: an accession number, a filename or the name of
                        a pre-existing bowtie2 index (created with
                        bowtie2-build). If not given and --reference is used,
                        the reference will be used to build a bowtie2 index.
                        (default: None)
  --ignoreIndex IGNOREDINDICES
                        Either: an accession number, a filename or the name of
                        a pre-existing bowtie2 index (created with
                        bowtie2-build). Reads matching this index will be
                        ignored. May be repeated. (default: None)
  --fastq1, -1 FASTQ1   The FASTQ reads to match against the bowtie2 index
                        given by --index. Also use --fast2 if you have paired
                        reads. If not given, single-end FASTQ reads will be
                        read from standard input. (default: None)
  --fastq2, -2 FASTQ2   The FASTQ reads to match against the bowtie2 index
                        given by --index. Use this with --fastq1 to specify
                        the mate file for paired-end reads. (default: None)
  --bowtie2Args BOWTIE2ARGS
                        Extra arguments to be passed to Bowtie2 (use --threads
                        to specify a thread count). (default: --no-unal)
  --samtoolsViewArgs SAMTOOLSVIEWARGS
                        Arguments to be passed to samtools view to create the
                        BAM file. (default: -F 3332 -q 30)
  --tempdir TEMPDIR     The temporary directory to use. If not specified, the
                        value of the TMPDIR environment variable (if any) is
                        used, or else /tmp. (default: None)
  --out, -o OUT         The output file name. If not given, the resulting SAM
                        or BAM will be written to standard output will be
                        used. (default: None)
  --reference REFERENCE
                        The reference FASTA file for use with
                        --callHaplotypesGATK and --callHaplotypesBcftools.
                        This will be used to build a Bowtie2 index if --index
                        is not given. (default: None)
  --vcfFile VCFFILE     The file to write VCF info to if --callHaplotypesGATK
                        or --callHaplotypesBcftools are used. (default: None)
  --markDuplicatesGATK  Use GATK to mark duplicates. See
                        https://gatk.broadinstitute.org for details on GATK.
                        (default: False)
  --markDuplicatesPicard
                        Use Picard to mark duplicates. See
                        https://github.com/broadinstitute/picard for details
                        on Picard. (default: False)
  --picardJar PICARDJAR
                        The path to the Picard jar file. See
                        https://github.com/broadinstitute/picard for details
                        on Picard. (default: None)
  --tmpChmod TMPCHMOD   A chmod string for setting the permission on the
                        temporary directory created. This will be passed to
                        chmod(1), so you could specify "g+rwx", for example.
                        (default: None)
  --removeDuplicates    Remove duplicates from the resulting SAM/BAM file.
                        Best used in combination with an option that marks
                        duplicates, such as --markDuplicatesGATK. (default:
                        False)
  --verbose             Print a description of commands as they are (or would
                        be, if --dryRun is used) executed. (default: False)
  --log                 Show a log of commands that were (or would be, if
                        --dryRun is used) executed. (default: False)
  --threads THREADS     The number of threads to use when running bowtie2
                        commands. (default: None)
  --noAlign             Do not align with Bowtie2, just build an index.
                        (default: True)
  --noBAM               Do not convert SAM to BAM. (default: True)
  --noSort              Do not sort the BAM. (default: True)
  --noIndexBAM          Do not index the BAM file. (default: True)
  --noClean             Do not remove intermediate files or the temporary
                        directory. (default: True)
  --force               Overwrite pre-existing output file. (default: False)
  --dryRun              Do not run commands, just print what would be done.
                        (default: False)
  --callHaplotypesGATK  Use GATK to call haplotypes. See
                        https://gatk.broadinstitute.org for details on GATK.
                        (default: False)
  --callHaplotypesBcftools
                        Use bcftools call to call haplotypes. (default: False)
  --removePrimersFromBedFile REMOVEPRIMERSFROMBEDFILE
                        If a bed file with Primers is specified, the Primers
                        will be soft-clipped from the bam file using iVar
                        (default: None)
```

</details>

<details><summary><strong><code>bwa-index-to-fasta.py</code></strong> — Extract FASTA from BWA index files (<code>.pac</code>/<code>.ann</code>/<code>.amb</code>). <span title="Show --help">ℹ️</span></summary>

```
usage: bwa-index-to-fasta.py [-h] [-o OUTPUT] [--lineWidth LINEWIDTH] prefix

Extract FASTA from BWA index files (.pac/.ann/.amb)

positional arguments:
  prefix                BWA index prefix (e.g. ref.fa — expects ref.fa.ann,
                        ref.fa.pac, ref.fa.amb)

options:
  -h, --help            show this help message and exit
  -o, --output OUTPUT   Output FASTA file (default: stdout)
  --lineWidth LINEWIDTH
                        Line width to use for FASTA file.
```

</details>

<details><summary><strong><code>filter-sam.py</code></strong> — Given a SAM/BAM file and a set of filtering criteria, write filtered SAM/BAM to stdout. <span title="Show --help">ℹ️</span></summary>

```
usage: filter-sam.py [-h] [--quiet] [--bam]
                     [--checkResultCount CHECKRESULTCOUNT] [--minLength N]
                     [--maxLength N] [--maxNFraction N]
                     [--whitelist SEQUENCE-ID] [--blacklist SEQUENCE-ID]
                     [--whitelistFile SEQUENCE-ID-FILE]
                     [--blacklistFile SEQUENCE-ID-FILE] [--titleRegex REGEX]
                     [--negativeTitleRegex REGEX]
                     [--sequenceWhitelist SEQUENCE]
                     [--sequenceBlacklist SEQUENCE]
                     [--sequenceWhitelistFile SEQUENCE-FILE]
                     [--sequenceBlacklistFile SEQUENCE-FILE]
                     [--sequenceRegex REGEX] [--sequenceNegativeRegex REGEX]
                     [--keepSequences NUMBER,RANGE,... |
                     --removeSequences NUMBER,RANGE,...] [--head N]
                     [--removeDuplicates | --removeDuplicatesById]
                     [--removeDuplicatesUseMD5] [--randomSubset N]
                     [--trueLength N] [--sampleFraction FRACTION]
                     [--sequenceNumbersFile FILENAME] --samfile SAMFILE
                     [--referenceId ID] [--dropUnmapped] [--dropSecondary]
                     [--dropSupplementary] [--dropDuplicates]
                     [--keepQCFailures] [--minScore FLOAT] [--maxScore FLOAT]
                     [--scoreTag TAG]

Given a SAM/BAM file and a set of filtering criteria write filtered SAM/BAM to
stdout.

options:
  -h, --help            show this help message and exit
  --quiet               If True, do not print the final summary. (default:
                        False)
  --bam                 If given, write (gzip compressed) BAM output.
                        (default: )
  --checkResultCount CHECKRESULTCOUNT
                        The number of alignments expected in the output. If
                        this number is not seen, the script exits with status
                        1 (and an error message is also printed, unless
                        --quiet was used). (default: None)
  --minLength N         The minimum sequence length. (default: None)
  --maxLength N         The maximum sequence length. (default: None)
  --maxNFraction N      The maximum fraction of Ns that can be present in the
                        sequence. (default: None)
  --whitelist SEQUENCE-ID
                        Sequence titles (ids) that should be whitelisted. May
                        be repeated. (default: None)
  --blacklist SEQUENCE-ID
                        Sequence titles (ids) that should be blacklisted. May
                        be repeated. (default: None)
  --whitelistFile SEQUENCE-ID-FILE
                        The name of a file that contains sequence titles (ids)
                        that should be whitelisted, one per line (default:
                        None)
  --blacklistFile SEQUENCE-ID-FILE
                        The name of a file that contains sequence titles (ids)
                        that should be blacklisted, one per line. (default:
                        None)
  --titleRegex REGEX    A regex that sequence titles (ids) must match.
                        (default: None)
  --negativeTitleRegex REGEX
                        A regex that sequence titles (ids) must not match.
                        (default: None)
  --sequenceWhitelist SEQUENCE
                        Sequences that should be whitelisted. May be repeated.
                        (default: None)
  --sequenceBlacklist SEQUENCE
                        Sequences that should be blacklisted. May be repeated.
                        (default: None)
  --sequenceWhitelistFile SEQUENCE-FILE
                        The name of a file that contains sequences that should
                        be whitelisted, one per line. (default: None)
  --sequenceBlacklistFile SEQUENCE-FILE
                        The name of a file that contains sequences that should
                        be blacklisted, one per line. (default: None)
  --sequenceRegex REGEX
                        A regex that sequences must match. (default: None)
  --sequenceNegativeRegex REGEX
                        A regex that sequences must not match. (default: None)
  --keepSequences NUMBER,RANGE,...
                        Specify (1-based) ranges of sequence numbers that
                        should be kept. E.g., --keepSequences 1-3,5 will
                        output just the 1st, 2nd, 3rd, and 5th sequences. All
                        others will be omitted. (default: None)
  --removeSequences NUMBER,RANGE,...
                        Specify (1-based) ranges of sequence numbers that
                        should be removed. E.g., --removeSequences 1-3,5 will
                        output all but the 1st, 2nd, 3rd, and 5th sequences.
                        All others will be ouput. (default: None)
  --head N              Only the first N sequences will be printed. (default:
                        None)
  --removeDuplicates    Duplicate reads will be removed, based only on
                        sequence identity. The first occurrence is kept.
                        (default: False)
  --removeDuplicatesById
                        Duplicate reads will be removed, based only on read
                        id. The first occurrence is kept. (default: False)
  --removeDuplicatesUseMD5
                        MD5 sums will be stored instead of the full sequence
                        or read id when either --removeDuplicates or
                        removeDuplicatesById are given. One of
                        --removeDuplicates and --removeDuplicatesById must
                        also be given. Note that this makes duplicate removal
                        probabilistic. This option can be used to reduce the
                        amount of RAM consumed during duplicate removal.
                        (default: False)
  --randomSubset N      An integer giving the number of sequences that should
                        be kept. These will be selected at random. (default:
                        None)
  --trueLength N        The number of reads in the FASTA input. Only to be
                        used with randomSubset (default: None)
  --sampleFraction FRACTION
                        A [0.0, 1.0] C{float} indicating a fraction of the
                        reads that should be allowed to pass through the
                        filter. The sample size will only be approximately the
                        product of the sample fraction and the number of
                        reads. The sample is taken at random. (default: None)
  --sequenceNumbersFile FILENAME
                        A file of (1-based) sequence numbers to retain.
                        Numbers must be one per line. (default: None)
  --samfile SAMFILE     The SAM/BAM file to filter.
  --referenceId ID      A reference sequence id whose alignments should be
                        kept (alignments against other references will be
                        dropped). If omitted, alignments against all
                        references will be kept. May be repeated. (default:
                        None)
  --dropUnmapped        If given, unmapped matches will not be output.
                        (default: False)
  --dropSecondary       If given, secondary matches will not be output.
                        (default: False)
  --dropSupplementary   If given, supplementary matches will not be output.
                        (default: False)
  --dropDuplicates      If given, matches flagged as optical or PCR duplicates
                        will not be output. (default: False)
  --keepQCFailures      If given, reads that are considered quality control
                        failures will be included in the output. (default:
                        False)
  --minScore FLOAT      If given, alignments with --scoreTag (default AS)
                        values less than this value will not be output. If
                        given, alignments that do not have a score will not be
                        output. (default: None)
  --maxScore FLOAT      If given, alignments with --scoreTag (default AS)
                        values greater than this value will not be output. If
                        given, alignments that do not have a score will not be
                        output. (default: None)
  --scoreTag TAG        The alignment tag to extract for --minScore and
                        --maxScore comparisons. (default: AS)
```

</details>

<details><summary><strong><code>sam-merge.py</code></strong> — Concatenate, sort, and index SAM/BAM files with <code>samtools</code>. <span title="Show --help">ℹ️</span></summary>

```
usage: sam-merge.py [-h] --bam BAM [-f] [-@ THREADS] FILE [FILE ...]

Concatenate, sort, and index SAM/BAM files with samtools.

positional arguments:
  FILE                  Input SAM/BAM files to merge.

options:
  -h, --help            show this help message and exit
  --bam BAM             Output BAM filename.
  -f, --force           Overwrite the output file if it already exists.
  -@, --threads THREADS
                        Number of threads to use for sorting/indexing
                        (default: 1).
```

</details>

<details><summary><strong><code>sam-references.py</code></strong> — Print SAM/BAM file reference names and lengths. <span title="Show --help">ℹ️</span></summary>

```
usage: sam-references.py [-h] FILENAME

Print SAM/BAM file reference names and lengths.

positional arguments:
  FILENAME    The name of a SAM/BAM alignment file.

options:
  -h, --help  show this help message and exit
```

</details>

<details><summary><strong><code>sam-reference-read-counts.py</code></strong> — Print SAM/BAM file reference-sequence and unmapped-read counts. <span title="Show --help">ℹ️</span></summary>

```
usage: sam-reference-read-counts.py [-h] [--sortBy {count,coverage,name}]
                                    [--excludeZeroes]
                                    [--excludeIfNoAdditional]
                                    [--topReferenceIdsFile FILE]
                                    FILENAME

Print SAM/BAM file reference sequence and unmapped read counts.

positional arguments:
  FILENAME              The name of a SAM/BAM alignment file.

options:
  -h, --help            show this help message and exit
  --sortBy {count,coverage,name}
                        Use 'count' to sort the output by decreasing read
                        count (i.e., the reference with the highest number of
                        matching reads is printed first), or 'coverage' to
                        sort by decreasing genome coverage, or 'name' to sort
                        by reference name. (default: count)
  --excludeZeroes       Do not print references unmatched by any read.
                        (default: False)
  --excludeIfNoAdditional
                        Do not print references that are mapped by no
                        additional reads (i.e., reads other than those
                        matching already-printed references). This is useful
                        when --sortBy count is used, making it possible to not
                        print references only matched by reads matching
                        references that have already been printed. (default:
                        False)
  --topReferenceIdsFile FILE
                        The file to write the (sorted) ids of the reads for
                        the reference with the highest number of matching
                        reads. Only valid when --sortBy count is used.
                        (default: None)
```

</details>

<details><summary><strong><code>sam-coverage.py</code></strong> — Print SAM/BAM file coverage statistics. <span title="Show --help">ℹ️</span></summary>

```
usage: sam-coverage.py [-h] [--noFilter] [--includeZeroes] [--filter]
                       [--regex] [--match] [--listReferences]
                       [--minReferenceLength MINREFERENCELENGTH] [--tsv TSV]
                       [--excel EXCEL] [--minLength N] [--maxLength N]
                       [--maxNFraction N] [--whitelist SEQUENCE-ID]
                       [--blacklist SEQUENCE-ID]
                       [--whitelistFile SEQUENCE-ID-FILE]
                       [--blacklistFile SEQUENCE-ID-FILE] [--titleRegex REGEX]
                       [--negativeTitleRegex REGEX]
                       [--sequenceWhitelist SEQUENCE]
                       [--sequenceBlacklist SEQUENCE]
                       [--sequenceWhitelistFile SEQUENCE-FILE]
                       [--sequenceBlacklistFile SEQUENCE-FILE]
                       [--sequenceRegex REGEX] [--sequenceNegativeRegex REGEX]
                       [--keepSequences NUMBER,RANGE,... |
                       --removeSequences NUMBER,RANGE,...] [--head N]
                       [--removeDuplicates | --removeDuplicatesById]
                       [--removeDuplicatesUseMD5] [--randomSubset N]
                       [--trueLength N] [--sampleFraction FRACTION]
                       [--sequenceNumbersFile FILENAME] [--referenceId ID]
                       [--dropUnmapped] [--dropSecondary]
                       [--dropSupplementary] [--dropDuplicates]
                       [--keepQCFailures] [--minScore FLOAT]
                       [--maxScore FLOAT] [--scoreTag TAG]
                       samfile

Print SAM/BAM file coverage statistics.

positional arguments:
  samfile               The SAM/BAM file to filter.

options:
  -h, --help            show this help message and exit
  --noFilter            Do not use our SAM filtering. This now the default
                        since filtering can be extremely slow. Use --filter to
                        turn filtering on. Note that if you use this option,
                        any filtering option (other than --referenceId) you
                        also specify that is provided by the
                        SAMFilter.addFilteringOptions will be silently
                        ignored! (default: False)
  --includeZeroes       Also write output for references that have no covering
                        reads. (default: False)
  --filter              Use our SAM filtering. Note that this can be extremely
                        slow. (default: False)
  --regex               Interpret the --referenceId arguments as regular
                        expressions. (default: False)
  --match               Use 'match' when checking if references match regular
                        expressions. This will cause the regular expression to
                        be checked against the start of the reference name. If
                        not given, the regular expression must match the
                        entire reference name not just its prefix. (default:
                        False)
  --listReferences      Print the list of matching reference IDs and then
                        exit. This makes it easy to adjust reference regular
                        expressions to see what work would be done without
                        doing it. (default: False)
  --minReferenceLength MINREFERENCELENGTH
                        The minimum length reference to use. This can be
                        useful if you give a reference regular expression that
                        matches too many references and you want to also limit
                        them by length. Or you can just give a minimal length
                        and that way pick up all chromosomes without needing
                        to name them. (default: None)
  --tsv TSV             The name of a TSV file to save results to. (default:
                        None)
  --excel EXCEL         The name of an Excel file to save results to.
                        (default: None)
  --minLength N         The minimum sequence length. (default: None)
  --maxLength N         The maximum sequence length. (default: None)
  --maxNFraction N      The maximum fraction of Ns that can be present in the
                        sequence. (default: None)
  --whitelist SEQUENCE-ID
                        Sequence titles (ids) that should be whitelisted. May
                        be repeated. (default: None)
  --blacklist SEQUENCE-ID
                        Sequence titles (ids) that should be blacklisted. May
                        be repeated. (default: None)
  --whitelistFile SEQUENCE-ID-FILE
                        The name of a file that contains sequence titles (ids)
                        that should be whitelisted, one per line (default:
                        None)
  --blacklistFile SEQUENCE-ID-FILE
                        The name of a file that contains sequence titles (ids)
                        that should be blacklisted, one per line. (default:
                        None)
  --titleRegex REGEX    A regex that sequence titles (ids) must match.
                        (default: None)
  --negativeTitleRegex REGEX
                        A regex that sequence titles (ids) must not match.
                        (default: None)
  --sequenceWhitelist SEQUENCE
                        Sequences that should be whitelisted. May be repeated.
                        (default: None)
  --sequenceBlacklist SEQUENCE
                        Sequences that should be blacklisted. May be repeated.
                        (default: None)
  --sequenceWhitelistFile SEQUENCE-FILE
                        The name of a file that contains sequences that should
                        be whitelisted, one per line. (default: None)
  --sequenceBlacklistFile SEQUENCE-FILE
                        The name of a file that contains sequences that should
                        be blacklisted, one per line. (default: None)
  --sequenceRegex REGEX
                        A regex that sequences must match. (default: None)
  --sequenceNegativeRegex REGEX
                        A regex that sequences must not match. (default: None)
  --keepSequences NUMBER,RANGE,...
                        Specify (1-based) ranges of sequence numbers that
                        should be kept. E.g., --keepSequences 1-3,5 will
                        output just the 1st, 2nd, 3rd, and 5th sequences. All
                        others will be omitted. (default: None)
  --removeSequences NUMBER,RANGE,...
                        Specify (1-based) ranges of sequence numbers that
                        should be removed. E.g., --removeSequences 1-3,5 will
                        output all but the 1st, 2nd, 3rd, and 5th sequences.
                        All others will be ouput. (default: None)
  --head N              Only the first N sequences will be printed. (default:
                        None)
  --removeDuplicates    Duplicate reads will be removed, based only on
                        sequence identity. The first occurrence is kept.
                        (default: False)
  --removeDuplicatesById
                        Duplicate reads will be removed, based only on read
                        id. The first occurrence is kept. (default: False)
  --removeDuplicatesUseMD5
                        MD5 sums will be stored instead of the full sequence
                        or read id when either --removeDuplicates or
                        removeDuplicatesById are given. One of
                        --removeDuplicates and --removeDuplicatesById must
                        also be given. Note that this makes duplicate removal
                        probabilistic. This option can be used to reduce the
                        amount of RAM consumed during duplicate removal.
                        (default: False)
  --randomSubset N      An integer giving the number of sequences that should
                        be kept. These will be selected at random. (default:
                        None)
  --trueLength N        The number of reads in the FASTA input. Only to be
                        used with randomSubset (default: None)
  --sampleFraction FRACTION
                        A [0.0, 1.0] C{float} indicating a fraction of the
                        reads that should be allowed to pass through the
                        filter. The sample size will only be approximately the
                        product of the sample fraction and the number of
                        reads. The sample is taken at random. (default: None)
  --sequenceNumbersFile FILENAME
                        A file of (1-based) sequence numbers to retain.
                        Numbers must be one per line. (default: None)
  --referenceId ID      A reference sequence id whose alignments should be
                        kept (alignments against other references will be
                        dropped). If omitted, alignments against all
                        references will be kept. May be repeated. (default:
                        None)
  --dropUnmapped        If given, unmapped matches will not be output.
                        (default: False)
  --dropSecondary       If given, secondary matches will not be output.
                        (default: False)
  --dropSupplementary   If given, supplementary matches will not be output.
                        (default: False)
  --dropDuplicates      If given, matches flagged as optical or PCR duplicates
                        will not be output. (default: False)
  --keepQCFailures      If given, reads that are considered quality control
                        failures will be included in the output. (default:
                        False)
  --minScore FLOAT      If given, alignments with --scoreTag (default AS)
                        values less than this value will not be output. If
                        given, alignments that do not have a score will not be
                        output. (default: None)
  --maxScore FLOAT      If given, alignments with --scoreTag (default AS)
                        values greater than this value will not be output. If
                        given, alignments that do not have a score will not be
                        output. (default: None)
  --scoreTag TAG        The alignment tag to extract for --minScore and
                        --maxScore comparisons. (default: AS)
```

</details>

<details><summary><strong><code>sam-coverage-depth.py</code></strong> — Print SAM/BAM per-site nucleotide composition and coverage statistics
(sites with no coverage produce no output). <span title="Show --help">ℹ️</span></summary>

```
usage: sam-coverage-depth.py [-h] [--referenceId ID] [--fasta FASTA]
                             [--diffsFrom {reference,commonest}]
                             [--window WINDOW] [--workers WORKERS]
                             [--noOffsets] [--noStats] [--statsFile STATSFILE]
                             [--verbosity VERBOSITY] [--minDepth MINDEPTH]
                             [--maxDepth MAXDEPTH]
                             [--minSiteMutationRate MINSITEMUTATIONRATE]
                             [--maxSiteMutationRate MAXSITEMUTATIONRATE]
                             [--minSiteMutationCount MINSITEMUTATIONCOUNT]
                             [--minEntropy MINENTROPY]
                             [--maxEntropy MAXENTROPY] [--noFilter] [--tsv]
                             [--sep SEP] [--minBases {1,2,3,4}]
                             [--maxIdenticalReadIdsPerSite MAXIDENTICALREADIDSPERSITE]
                             [--minBaseQuality MINBASEQUALITY]
                             [--minMapQuality MINMAPQUALITY]
                             samfile

Print SAM/BAM per-site nucleotide composition and coverage statistics. Note
that no output is produced for sites with no read coverage.

positional arguments:
  samfile               The SAM/BAM file to read.

options:
  -h, --help            show this help message and exit
  --referenceId ID      A reference sequence id. If not given and the SAM/BAM
                        file has just one reference, that one will be used.
                        Otherwise, a reference id must be given. (default:
                        None)
  --fasta FASTA         Optionally give a FASTA file in which the reference
                        can be found. If you do not provide a reference
                        sequence, mutations will be counted as differences
                        from the most-common nucleotide at each genome site.
                        (default: None)
  --diffsFrom {reference,commonest}
                        The base from which to assess differences. If
                        'reference', differences from the reference base (at
                        each genome site) will be treated as mutations. If
                        'commonest', differences from the most common base
                        will be treated as mutations (with ties broken in
                        favour of the reference base, if a reference is
                        given). It is important that you understand the
                        difference between these two and choose the right
                        option. If not specified, mutations will not be
                        counted and the mutation-related columns in the --tsv
                        output will be zeroes. (default: None)
  --window, -n WINDOW   The number of genome sites to have each subprocess
                        handle. Large values will have a major performance
                        impact when --maxDepth is also high (e.g., >10,000).
                        This seems to be due to an unknown pysam/htslib
                        slowdown (based on samtools v1.21, as of 2025-01-31).
                        (default: 250)
  --workers WORKERS     The number of subprocesses to run in parallel. If not
                        specified, the process pool will typically use all
                        available cores, but see https://docs.python.org/3/lib
                        rary/concurrent.futures.html#processpoolexecutor for
                        details. (default: None)
  --noOffsets           Do not print per-site details of base counts.
                        (default: True)
  --noStats             Do not print final average, standard deviation etc.
                        statistics. (default: True)
  --statsFile STATSFILE
                        The file to print overall statistics to. If not
                        specified, standard output will be used. (default:
                        None)
  --verbosity VERBOSITY
                        Print progress information (to standard error). 0 =
                        don't print progress information. 1 = print
                        information after each genome region is processed. 2 =
                        also print information during the processing of each
                        region. (default: 0)
  --minDepth MINDEPTH   The minimum read depth to report on. Sites with a
                        lower depth will not be reported. (default: None)
  --maxDepth MAXDEPTH   The maximum read depth to request from pysam for each
                        site. (default: 1000000.0)
  --minSiteMutationRate MINSITEMUTATIONRATE
                        The minimum site mutation rate to report on. Sites
                        with a lower mutation rate will not be reported. The
                        mutation rate is the fraction of bases at the site
                        that are not the reference or commonest base
                        (depending on the --diffsFrom option, when given).
                        (default: None)
  --maxSiteMutationRate MAXSITEMUTATIONRATE
                        The maximum site mutation rate to report on. Sites
                        with a higher mutation rate will not be reported. The
                        mutation rate is the fraction of bases at the site
                        that are not the reference or commonest base
                        (depending on the --diffsFrom option, when given).
                        (default: None)
  --minSiteMutationCount MINSITEMUTATIONCOUNT
                        The minimum number of mutant nucleotides that must be
                        present at a site in order for it to be reported.
                        E.g., if this were set to 4, there would need to be at
                        least four nucleotides (of any mixture of bases) that
                        were not the consensus or reference (depending on the
                        --diffsFrom setting, when given). (default: 0)
  --minEntropy MINENTROPY
                        The minimum site entropy to report on. Sites with
                        lower entropy will not be reported. (default: 0.0)
  --maxEntropy MAXENTROPY
                        The maximum site entropy to report on. Sites with
                        higher entropy will not be reported. (default: 2.0)
  --noFilter            This behaviour is now the default and this option is
                        ignored. (default: True)
  --tsv                 Produce TSV output (use --sep , if you want CSV). Note
                        that overall SAM/BAM file statistics will be printed
                        at the bottom of the output unless you also specify
                        --noStats or --statsFile. Columns will hopefully be
                        self-explanatory, but a few clarifying comments will
                        probably be helpful. 1) The 'Base' column indicates
                        the nucleotide base that was used to count mutations.
                        The base used will depend on the value of --diffsFrom
                        (when given). See the help text for that option. 2)
                        The 'Reference' column will have an empty value if no
                        reference sequence was provided. 3) The 'Site' column
                        gives 1-based indexes into the reference genome.
                        (default: False)
  --sep SEP             The separator to use when --tsv is used. Ignored
                        otherwise. (default: )
  --minBases {1,2,3,4}  The minimum number of different bases that must occur
                        at a site for the site to be reported. E.g., if this
                        were set to 3, sites are only counted if at most one
                        nucleotide is not present. The default value of 1 will
                        do no filtering because all (covered) sites will have
                        at least one nucleotide present. Higher values of this
                        parameter can be used to add a diversity filter, by
                        only counting sites at which there is greater
                        nucleotide diversity. (default: 1)
  --maxIdenticalReadIdsPerSite MAXIDENTICALREADIDSPERSITE
                        The number of reads with identical IDs that are
                        permitted at each site. This should be 2 for paired-
                        end reads (Illumina sequencing uses the identical read
                        ID for both members of the pair) because the read pair
                        may overlap at some sites. If you have unpaired reads,
                        you could set this to 1. A factor that needs to be
                        considered is whether the mapping algorithm was told
                        to only produce the best match for each read (or read
                        pair). If so, setting this to 1 or 2 (as just
                        described) makes sense because it will catch
                        inadvertent errors. But if reads are allowed to match
                        multiple times, there could be many matches of the
                        same read in the same reference section. In that case,
                        you can use a value of -1 to disable the checking.
                        (default: 2)
  --minBaseQuality MINBASEQUALITY
                        The minimum nucleotide base quality. Note that
                        samtools depth defaults to 13. (default: 13)
  --minMapQuality MINMAPQUALITY
                        The minimum mapping quality (SAM MAPQ field).
                        (default: 0)
```

</details>

<details><summary><strong><code>plot-coverage-depth.py</code></strong> — Plot per-base coverage depth from a SAM/BAM file. <span title="Show --help">ℹ️</span></summary>

```
usage: plot-coverage-depth.py [-h] --out OUT [--reference REFERENCE]
                              [--title TITLE] [--xTitle XTITLE]
                              [--yTitle YTITLE] [--legendTitle LEGENDTITLE]
                              [--places PLACES]
                              [--titleFontSize TITLEFONTSIZE]
                              [--axisFontSize AXISFONTSIZE]
                              [--legendSquareSize LEGENDSQUARESIZE]
                              [--region START:END[:LABEL[:COLOR]]]
                              SAM-OR-BAM-FILE

Plot per-base coverage depth from a SAM/BAM file.

positional arguments:
  SAM-OR-BAM-FILE       The SAM or BAM file to read.

options:
  -h, --help            show this help message and exit
  --out OUT             Output filename. The file suffix determines the
                        format: .html for interactive HTML, or
                        .png/.svg/.pdf/.jpg for a static image.
  --reference REFERENCE
                        The reference sequence name to plot coverage for. If
                        not given, the reference is auto-detected from the
                        aligned reads (this fails if reads map to more than
                        one reference or no mapped reads are found). (default:
                        None)
  --title TITLE         The plot title. You can embed {median} {mean} {sd}
                        {min} and {max} to incorporate coverage statistics
                        into the title. The --places argument determines the
                        number of decimal places that will be used in the mean
                        and standard deviation. {reference} will embed the
                        name of the reference. (default: None)
  --xTitle XTITLE       The x-axis title. (default: Genome position)
  --yTitle YTITLE       The y-axis title. (default: Coverage depth)
  --legendTitle LEGENDTITLE
                        The legend title. (default: Regions)
  --places PLACES       Number of decimal places for floating-point coverage
                        statistics (median, mean, sd) that can be substituted
                        into the plot title. See the --title option for how to
                        do this. (default: 2)
  --titleFontSize TITLEFONTSIZE
                        The font size for the plot title. (default: 16)
  --axisFontSize AXISFONTSIZE
                        The font size for axis labels. (default: 14)
  --legendSquareSize LEGENDSQUARESIZE
                        The length of a legend colour square (for --region
                        arguments, if any). (default: 8)
  --region START:END[:LABEL[:COLOR]]
                        Highlight a genome region with a coloured bar below
                        y=0. START and END should be 1-based inclusive
                        offsets. LABEL is optional text shown in the legend
                        (default: region-1, region-2, …). COLOR is any
                        CSS/plotly colour name or hex code; if omitted a
                        colour is assigned automatically. Overlapping regions
                        will be stacked vertically. This argument may be given
                        multiple times. (default: None)
```

</details>

<details><summary><strong><code>sam-match-lengths.py</code></strong> — Print information about the lengths of query matches in a SAM/BAM file
(query length, matched-region length excluding soft clips, soft-clip
lengths on each end); can print per-query TSV and/or summary statistics. <span title="Show --help">ℹ️</span></summary>

```
usage: sam-match-lengths.py [-h] [--tsv] [--stats] samfile

Print information about the lengths of query matches in a SAM/BAM file. If
--tsv is given, print (to standard output) tab-separated fields containing:
the query id, query lengths, and matched region length (i.e., excluding soft-
clipped bases). The number of soft-clipped bases on the left and right are
also given in the last two TSV fields. Queries that did not align produce no
output. If --stats is given, summary statistics are printed to standard output
(or to standard error if --tsv is also given).

positional arguments:
  samfile     The SAM/BAM file to read.

options:
  -h, --help  show this help message and exit
  --tsv       Print TSV to standard output.
  --stats     Print summary statistics to standard output (or standard error
              if --tsv is given).
```

</details>

<details><summary><strong><code>genome-protein-summary.py</code></strong> — Print SAM/BAM file protein-match statistics. <span title="Show --help">ℹ️</span></summary>

```
usage: genome-protein-summary.py [-h]
                                 --proteinGenomeDatabase PROTEINGENOMEDATABASE
                                 [--progress]
                                 [--sortOn {coverage,depth,name,offset,readCount}]
                                 [--minReadOffsetCount MINREADOFFSETCOUNT]
                                 [--skipTranslationChecks] [--minLength N]
                                 [--maxLength N] [--maxNFraction N]
                                 [--whitelist SEQUENCE-ID]
                                 [--blacklist SEQUENCE-ID]
                                 [--whitelistFile SEQUENCE-ID-FILE]
                                 [--blacklistFile SEQUENCE-ID-FILE]
                                 [--titleRegex REGEX]
                                 [--negativeTitleRegex REGEX]
                                 [--sequenceWhitelist SEQUENCE]
                                 [--sequenceBlacklist SEQUENCE]
                                 [--sequenceWhitelistFile SEQUENCE-FILE]
                                 [--sequenceBlacklistFile SEQUENCE-FILE]
                                 [--sequenceRegex REGEX]
                                 [--sequenceNegativeRegex REGEX]
                                 [--keepSequences NUMBER,RANGE,... |
                                 --removeSequences NUMBER,RANGE,...]
                                 [--head N] [--removeDuplicates |
                                 --removeDuplicatesById]
                                 [--removeDuplicatesUseMD5] [--randomSubset N]
                                 [--trueLength N] [--sampleFraction FRACTION]
                                 [--sequenceNumbersFile FILENAME]
                                 [--samfile [SAMFILE ...]] [--referenceId ID]
                                 [--dropUnmapped] [--dropSecondary]
                                 [--dropSupplementary] [--dropDuplicates]
                                 [--keepQCFailures] [--minScore FLOAT]
                                 [--maxScore FLOAT] [--scoreTag TAG]

Print SAM/BAM file protein match statistics.

options:
  -h, --help            show this help message and exit
  --proteinGenomeDatabase PROTEINGENOMEDATABASE
                        The filename of an Sqlite3 database holding protein
                        and genome information, as built by make-protein-
                        database.py
  --progress            Print progress info to standard error. (default:
                        False)
  --sortOn {coverage,depth,name,offset,readCount}
                        How to sort proteins for output. (default: readCount)
  --minReadOffsetCount MINREADOFFSETCOUNT
                        The minimum number of reads offsets that must overlap
                        a protein for the read to be considered as
                        sufficiently intersecting the protein. Use this to
                        prevent reads that just overlap the protein in a very
                        small number offsets from being counted. (default:
                        None)
  --skipTranslationChecks
                        Skip the sanity check that database protein sequences
                        can all be translated from the database genome
                        sequence. (default: True)
  --minLength N         The minimum sequence length. (default: None)
  --maxLength N         The maximum sequence length. (default: None)
  --maxNFraction N      The maximum fraction of Ns that can be present in the
                        sequence. (default: None)
  --whitelist SEQUENCE-ID
                        Sequence titles (ids) that should be whitelisted. May
                        be repeated. (default: None)
  --blacklist SEQUENCE-ID
                        Sequence titles (ids) that should be blacklisted. May
                        be repeated. (default: None)
  --whitelistFile SEQUENCE-ID-FILE
                        The name of a file that contains sequence titles (ids)
                        that should be whitelisted, one per line (default:
                        None)
  --blacklistFile SEQUENCE-ID-FILE
                        The name of a file that contains sequence titles (ids)
                        that should be blacklisted, one per line. (default:
                        None)
  --titleRegex REGEX    A regex that sequence titles (ids) must match.
                        (default: None)
  --negativeTitleRegex REGEX
                        A regex that sequence titles (ids) must not match.
                        (default: None)
  --sequenceWhitelist SEQUENCE
                        Sequences that should be whitelisted. May be repeated.
                        (default: None)
  --sequenceBlacklist SEQUENCE
                        Sequences that should be blacklisted. May be repeated.
                        (default: None)
  --sequenceWhitelistFile SEQUENCE-FILE
                        The name of a file that contains sequences that should
                        be whitelisted, one per line. (default: None)
  --sequenceBlacklistFile SEQUENCE-FILE
                        The name of a file that contains sequences that should
                        be blacklisted, one per line. (default: None)
  --sequenceRegex REGEX
                        A regex that sequences must match. (default: None)
  --sequenceNegativeRegex REGEX
                        A regex that sequences must not match. (default: None)
  --keepSequences NUMBER,RANGE,...
                        Specify (1-based) ranges of sequence numbers that
                        should be kept. E.g., --keepSequences 1-3,5 will
                        output just the 1st, 2nd, 3rd, and 5th sequences. All
                        others will be omitted. (default: None)
  --removeSequences NUMBER,RANGE,...
                        Specify (1-based) ranges of sequence numbers that
                        should be removed. E.g., --removeSequences 1-3,5 will
                        output all but the 1st, 2nd, 3rd, and 5th sequences.
                        All others will be ouput. (default: None)
  --head N              Only the first N sequences will be printed. (default:
                        None)
  --removeDuplicates    Duplicate reads will be removed, based only on
                        sequence identity. The first occurrence is kept.
                        (default: False)
  --removeDuplicatesById
                        Duplicate reads will be removed, based only on read
                        id. The first occurrence is kept. (default: False)
  --removeDuplicatesUseMD5
                        MD5 sums will be stored instead of the full sequence
                        or read id when either --removeDuplicates or
                        removeDuplicatesById are given. One of
                        --removeDuplicates and --removeDuplicatesById must
                        also be given. Note that this makes duplicate removal
                        probabilistic. This option can be used to reduce the
                        amount of RAM consumed during duplicate removal.
                        (default: False)
  --randomSubset N      An integer giving the number of sequences that should
                        be kept. These will be selected at random. (default:
                        None)
  --trueLength N        The number of reads in the FASTA input. Only to be
                        used with randomSubset (default: None)
  --sampleFraction FRACTION
                        A [0.0, 1.0] C{float} indicating a fraction of the
                        reads that should be allowed to pass through the
                        filter. The sample size will only be approximately the
                        product of the sample fraction and the number of
                        reads. The sample is taken at random. (default: None)
  --sequenceNumbersFile FILENAME
                        A file of (1-based) sequence numbers to retain.
                        Numbers must be one per line. (default: None)
  --samfile [SAMFILE ...]
                        The SAM/BAM file to filter. (default: None)
  --referenceId ID      A reference sequence id whose alignments should be
                        kept (alignments against other references will be
                        dropped). If omitted, alignments against all
                        references will be kept. May be repeated. (default:
                        None)
  --dropUnmapped        If given, unmapped matches will not be output.
                        (default: False)
  --dropSecondary       If given, secondary matches will not be output.
                        (default: False)
  --dropSupplementary   If given, supplementary matches will not be output.
                        (default: False)
  --dropDuplicates      If given, matches flagged as optical or PCR duplicates
                        will not be output. (default: False)
  --keepQCFailures      If given, reads that are considered quality control
                        failures will be included in the output. (default:
                        False)
  --minScore FLOAT      If given, alignments with --scoreTag (default AS)
                        values less than this value will not be output. If
                        given, alignments that do not have a score will not be
                        output. (default: None)
  --maxScore FLOAT      If given, alignments with --scoreTag (default AS)
                        values greater than this value will not be output. If
                        given, alignments that do not have a score will not be
                        output. (default: None)
  --scoreTag TAG        The alignment tag to extract for --minScore and
                        --maxScore comparisons. (default: AS)
```

</details>

<details><summary><strong><code>extract-reference-from-fasta-and-bam.py</code></strong> — Extract a reference sequence from a FASTA file and its matches from a BAM file. <span title="Show --help">ℹ️</span></summary>

```
usage: extract-reference-from-fasta-and-bam.py [-h] --reference-id ID
                                               --out-dir DIR --bam BAM
                                               --fasta FASTA [--quiet]
                                               [--keep-description] [--force]
                                               [--index-fasta] [--keep-sam]

Extract a reference sequence from a FASTA file and its matches from a BAM
file.

options:
  -h, --help          show this help message and exit
  --reference-id ID   The ID of the reference sequence.
  --out-dir DIR       The directory into which to put the output files.
  --bam, --sam BAM    The BAM or SAM file.
  --fasta FASTA       The FASTA file.
  --quiet             Do not print a summary.
  --keep-description  Do not ignore the description (if any) of the sequence
                      in the FASTA file. Normally, anything after the first
                      space in the ID in the FASTA file will be ignored (for
                      finding the ID in the FASTA file, in the output FASTA,
                      and will not be included in the output BAM (where it is
                      illegal). If you use this option, the description will
                      not be ignored, so the --reference-id you give will be
                      used in full. This may result in not finding what you
                      are looking for (in the input FASTA).
  --force             Overwrite pre-existing output files.
  --index-fasta       Index the extracted FASTA (using samtools faidx).
  --keep-sam          Do not remove the intermediate SAM file with the
                      matching records (produced en route to making a BAM
                      file)
```

</details>

<details><summary><strong><code>reference-read-scores-to-JSON.py</code></strong> — Print a JSON object containing reference-to-read distances extracted
from a SAM file. <span title="Show --help">ℹ️</span></summary>

```
usage: reference-read-scores-to-JSON.py [-h] --samFile SAMFILE
                                        [--minMatchingReads MINMATCHINGREADS]
                                        [--scoreTag SCORETAG] [--verbose]

Print a JSON object containing reference to read distances extracted from a
SAM file.

options:
  -h, --help            show this help message and exit
  --samFile SAMFILE     The SAM file(s) to load. May be repeated.
  --minMatchingReads MINMATCHINGREADS
                        The minimum number of reads that must match a
                        reference for it to be included.
  --scoreTag SCORETAG   The score tag to use for the alignment score. If not
                        given, 1 will be used to indicate that a read matched
                        a reference (non-matches are not included). The
                        default is no score tag, which is not that useful. A
                        good choice is "AS", for the alignment score, but that
                        has to be present in the SAM file, which means that
                        the aligner (bowtie2, bwa, etc. has to have produced
                        such a tag.
  --verbose             Print extra information.
```

</details>

<details><summary><strong><code>plot-references-by-inter-read-distance.py</code></strong> — Plot references according to their pairwise distances, computed from
shared/non-shared read matches in a SAM file. <span title="Show --help">ℹ️</span></summary>

```
usage: plot-references-by-inter-read-distance.py [-h] [--samFile SAMFILE]
                                                 [--scoreFile SCOREFILE]
                                                 --htmlFile HTMLFILE
                                                 [--mdsFile MDSFILE]
                                                 [--categoryName CATEGORYNAME]
                                                 [--categorySortKey {HBV genotype}]
                                                 [--categoryFile CATEGORYFILE]
                                                 [--minMatchingReads MINMATCHINGREADS]
                                                 [--verbose] [--2d | --logZ]
                                                 [--metric {jaccard,soergel}]

Plot references according to their distances from one another, based on common
read matches (and non-matches).

options:
  -h, --help            show this help message and exit
  --samFile SAMFILE     The SAM file(s) to load. May be repeated.
  --scoreFile SCOREFILE
                        A (JSON) file to read (JSON) the reference/read score
                        matrix from. This should have been produced by an
                        earlier call with a --scoreFile argument or else by
                        saving the output of reference-read-scores-to-JSON.py
                        If the file does not exist, it will be created.
  --htmlFile HTMLFILE   The HTML plotly file to write.
  --mdsFile MDSFILE     A JSON to cache the result of the MDS optimization.
  --categoryName CATEGORYNAME
                        The name of the categories to which reference
                        sequences are assigned (e.g., "Genotype"). This name
                        will appear at the top of the legend.
  --categorySortKey {HBV genotype}
                        The name of the function to use to sort categories
                        (for the legend). If not given, categories will be
                        sorted by name.
  --categoryFile CATEGORYFILE
                        A file containing labels for reference sequences. The
                        format is lines containing a reference name, a TAB,
                        and a new name (to be shown in the plot). If no label
                        file is given, or if there is no new name for a
                        reference, the original reference names will be used.
  --minMatchingReads MINMATCHINGREADS
                        The minimum number of reads that must match a
                        reference for it to be included.
  --verbose             Print extra information.
  --2d, --2D            Make the figure in two dimensions, instead of three.
                        The third dimensions (reference matching read count)
                        is omitted.
  --logZ                Log the z (read count) axis (only valid if --2d is not
                        used).
  --metric {jaccard,soergel}
                        The distance metric to use.
```

</details>

## 7. BLAST/DIAMOND hit processing and alignment panels

Scripts built around dark matter's BLAST/DIAMOND JSON record format (an
older, simplified JSON encoding of BLAST/DIAMOND hits used throughout this
part of the codebase).

<details><summary><strong><code>convert-blast-xml-to-json.py</code></strong> — Convert a BLAST XML file to dark matter's JSON format. <span title="Show --help">ℹ️</span></summary>

```
usage: convert-blast-xml-to-json.py [-h] [--json JSON-output-file]
                                    [--xml BLAST-XML-file] [--bzip2]

Convert a BLAST XML file to JSON.

options:
  -h, --help            show this help message and exit
  --json JSON-output-file
                        the JSON filename to write the converted XML to. If
                        omitted, standard output will be used.
  --xml BLAST-XML-file  the BLAST XML output file to convert. If omitted,
                        standard input will be read.
  --bzip2               If True, compress the JSON output using bzip2.

Give a BLAST XML file and convert it to JSON. Optionally compress the JSON
output.
```

</details>

<details><summary><strong><code>convert-diamond-to-json.py</code></strong> — Convert DIAMOND tabular output to dark matter's JSON format. <span title="Show --help">ℹ️</span></summary>

```
usage: convert-diamond-to-json.py [-h] [--json JSON-output-file]
                                  [--diamond DIAMOND-file] [--bzip2]

Convert DIAMOND tabular output to JSON.

options:
  -h, --help            show this help message and exit
  --json JSON-output-file
                        The JSON filename to write the converted DIAMOND file
                        to. If omitted, standard output will be used.
  --diamond DIAMOND-file
                        The DIAMOND tabular output file to convert. If
                        omitted, standard input will be read.
  --bzip2               Compress output using bzip2.

Give a DIAMOND tabular file and convert it to JSON, optionally compressing the
output. You *must* invoke DIAMOND with the following output specification:
--outfmt 6 qtitle stitle bitscore evalue qframe qseq qstart qend sseq sstart
send slen btopNote that only each line of the output is JSON (the full output
is not valid JSON by itself).
```

</details>

**`cat-json-blast-records.py`** — Read dark matter's JSON BLAST output (one JSON object per line) from
stdin and pretty-print it to stdout.

**`check-fasta-json-blast-consistency.py`** — Verify that a FASTA file and the JSON BLAST file produced from it are
mutually consistent; exits non-zero with a stderr message if not.

**`read-blast-json.py`** — Read dark matter's JSON BLAST records and report the elapsed time (a
benchmarking/profiling script, not an end-user tool).

**`read-blast-xml.py`** — Read dark matter's "simplified XML" BLAST records and report elapsed
time (same purpose as above, for an XML-based variant).

<details><summary><strong><code>bit-score-to-e-value.py</code></strong> — Convert a BLAST bit score to an e-value (see also <code>e-value-to-bit-score.py</code>). <span title="Show --help">ℹ️</span></summary>

```
usage: bit-score-to-e-value.py [-h] --bitScore BITSCORE --dbSize DBSIZE
                               --dbSequenceCount DBSEQUENCECOUNT
                               --queryLength QUERYLENGTH
                               --lengthAdjustment LENGTHADJUSTMENT

Convert a bit score to an e-value.

options:
  -h, --help            show this help message and exit
  --bitScore BITSCORE   The bit score to convert
  --dbSize DBSIZE       The total number of bases in the sequence database.
  --dbSequenceCount DBSEQUENCECOUNT
                        The number of sequences in the database.
  --queryLength QUERYLENGTH
                        The length of the query sequence.
  --lengthAdjustment LENGTHADJUSTMENT
                        The length adjustment.
```

</details>

<details><summary><strong><code>e-value-to-bit-score.py</code></strong> — Convert a BLAST e-value to a bit score (see also <code>bit-score-to-e-value.py</code>). <span title="Show --help">ℹ️</span></summary>

```
usage: e-value-to-bit-score.py [-h] --eValue EVALUE --dbSize DBSIZE
                               --dbSequenceCount DBSEQUENCECOUNT
                               --queryLength QUERYLENGTH
                               --lengthAdjustment LENGTHADJUSTMENT

Convert a bit score to an e-value.

options:
  -h, --help            show this help message and exit
  --eValue EVALUE       The e-value to convert.
  --dbSize DBSIZE       The total number of bases in the sequence database.
  --dbSequenceCount DBSEQUENCECOUNT
                        The number of sequences in the database.
  --queryLength QUERYLENGTH
                        The length of the query sequence.
  --lengthAdjustment LENGTHADJUSTMENT
                        The length adjustment.
```

</details>

<details><summary><strong><code>filter-hits-to-fasta.py</code></strong> — Extract FASTA from matching BLAST hits. <span title="Show --help">ℹ️</span></summary>

```
usage: filter-hits-to-fasta.py [-h] --fasta FASTA-file [--minStart MINSTART]
                               [--maxStop MAXSTOP] [--oneAlignmentPerRead]
                               [--scoreCutoff SCORECUTOFF]
                               [--maxHspsPerHit MAXHSPSPERHIT]
                               [--whitelist WHITELIST [WHITELIST ...]]
                               [--blacklist BLACKLIST [BLACKLIST ...]]
                               [--titleRegex TITLEREGEX]
                               [--negativeTitleRegex NEGATIVETITLEREGEX]
                               [--truncateTitlesAfter TRUNCATETITLESAFTER]
                               [--minSequenceLen MINSEQUENCELEN]
                               [--maxSequenceLen MAXSEQUENCELEN]
                               [--taxonomy TAXONOMY]
                               [--minMatchingReads MINMATCHINGREADS]
                               [--maxMatchingReads MAXMATCHINGREADS]
                               [--minMedianScore MINMEDIANSCORE]
                               [--withScoreBetterThan WITHSCOREBETTERTHAN]
                               [--minNewReads MINNEWREADS]
                               BLAST-JSON-file [BLAST-JSON-file ...]

Extract FASTA from matching BLAST hits.

positional arguments:
  BLAST-JSON-file       the JSON file of BLAST output.

options:
  -h, --help            show this help message and exit
  --fasta FASTA-file    the FASTA file of sequences that were given to BLAST.
  --minStart MINSTART   Reads that start before this subject offset should not
                        be shown.
  --maxStop MAXSTOP     Reads that end after this subject offset should not be
                        shown.
  --oneAlignmentPerRead
                        If C{True}, only keep the best alignment for each
                        read.
  --scoreCutoff SCORECUTOFF
                        A float score. Matches with scores worse than this
                        will be ignored.
  --maxHspsPerHit MAXHSPSPERHIT
                        A numeric max number of HSPs to show for each hit on
                        hitId.
  --whitelist WHITELIST [WHITELIST ...]
                        sequence titles that should be whitelisted
  --blacklist BLACKLIST [BLACKLIST ...]
                        sequence titles that should be blacklisted
  --titleRegex TITLEREGEX
                        a regex that sequence titles must match.
  --negativeTitleRegex NEGATIVETITLEREGEX
                        a regex that sequence titles must not match.
  --truncateTitlesAfter TRUNCATETITLESAFTER
                        a string that titles will be truncated beyond. If a
                        truncated title has already been seen, that title will
                        be skipped.
  --minSequenceLen MINSEQUENCELEN
                        sequences of lesser length will be elided.
  --maxSequenceLen MAXSEQUENCELEN
                        sequences of greater length will be elided.
  --taxonomy TAXONOMY   a string of the taxonomic group on which should be
                        filtered. eg "Vira" will filter on viruses.
  --minMatchingReads MINMATCHINGREADS
                        sequences that are matched by fewer reads will be
                        elided.
  --maxMatchingReads MAXMATCHINGREADS
                        sequences that are matched by more reads will be
                        elided.
  --minMedianScore MINMEDIANSCORE
                        sequences that are matched with a median score that is
                        worse will be elided.
  --withScoreBetterThan WITHSCOREBETTERTHAN
                        sequences that are matched without at least one score
                        at least this good will be elided.
  --minNewReads MINNEWREADS
                        The fraction of its reads by which a new read set must
                        differ from all previously seen read sets in order to
                        be considered acceptably different.

Given JSON BLAST output files, a FASTA sequence file, and filtering criteria,
produce a FASTA file on stdout containing the reads that match the criteria.
```

</details>

<details><summary><strong><code>filter-reads-alignments.py</code></strong> — Given BLAST/DIAMOND JSON output and the FASTA/FASTQ files that were
searched, print the <em>input reads</em> that pass a set of filtering criteria
(the complement of <code>filter-hits-to-fasta.py</code>, which extracts matched
<em>subjects</em> rather than queries). <span title="Show --help">ℹ️</span></summary>

```
usage: filter-reads-alignments.py [-h] [--matcher {blast,diamond}]
                                  --json FILE.JSON[.BZ2] [FILE.JSON[.BZ2] ...]
                                  (--fasta FILE.FASTA [FILE.FASTA ...] |
                                  --fastq FILE.FASTQ [FILE.FASTQ ...])
                                  [--diamondDatabaseFastaFilename FILE.FASTA |
                                  --diamondSqliteDatabaseFilename FILE.SQL]
                                  [--diamondDatabaseFastaDirectory DIR]
                                  [--minStart OFFSET] [--maxStop OFFSET]
                                  [--oneAlignmentPerRead]
                                  [--maxAlignmentsPerRead N]
                                  [--scoreCutoff SCORE] [--maxHspsPerHit N]
                                  [--whitelist TITLE [TITLE ...]]
                                  [--blacklist TITLE [TITLE ...]]
                                  [--titleRegex TITLE-REGEX]
                                  [--negativeTitleRegex TITLE-REGEX]
                                  [--truncateTitlesAfter TRUNCATETITLESAFTER]
                                  [--minSequenceLen N] [--maxSequenceLen N]
                                  [--taxonomy NAME]

Filter reads alignments to select a subset of reads. Given BLAST or DIAMOND
JSON output files and the FASTA (or FASTQ) query sequence files given to BLAST
or DIAMOND, and some filtering criteria, print the input reads that pass the
filtering.

options:
  -h, --help            show this help message and exit
  --matcher {blast,diamond}
                        The matching algorithm that was used to produce the
                        JSON. (default: blast)
  --json FILE.JSON[.BZ2] [FILE.JSON[.BZ2] ...]
                        the JSON file(s) of BLAST or DIAMOND output.
  --fasta FILE.FASTA [FILE.FASTA ...]
                        the FASTA file(s) of sequences that were given to
                        BLAST or DIAMOND. (default: None)
  --fastq FILE.FASTQ [FILE.FASTQ ...]
                        the FASTQ file(s) of sequences that were given to
                        BLAST or DIAMOND. (default: None)
  --diamondDatabaseFastaFilename FILE.FASTA
                        The filename of the FASTA file used to make the
                        DIAMOND database. If --matcher diamond is used, either
                        this argument or --diamondSqliteDatabaseFilename must
                        be specified. (default: None)
  --diamondSqliteDatabaseFilename FILE.SQL
                        The filename of the sqlite3 database file of FASTA
                        metadata, made from the FASTA that was used to make
                        the DIAMOND database. If --matcher diamond is used,
                        either this argument or --diamondDatabaseFilename must
                        be specified. (default: None)
  --diamondDatabaseFastaDirectory DIR
                        The directory where the FASTA file used to make the
                        DIAMOND database can be found. This argument is only
                        useful when --diamondSqliteDatabaseFilename is
                        specified. (default: None)
  --minStart OFFSET     Reads that start before this subject offset should not
                        be shown. (default: None)
  --maxStop OFFSET      Reads that end after this subject offset should not be
                        shown. (default: None)
  --oneAlignmentPerRead
                        If True, only keep the best alignment for each read.
                        (default: False)
  --maxAlignmentsPerRead N
                        Reads with more than this many alignments will be
                        elided. Pass zero to only keep reads with no matches
                        (alignments). (default: None)
  --scoreCutoff SCORE   A float score. Matches with scores worse than this
                        will be ignored. (default: None)
  --maxHspsPerHit N     A numeric maximum number of HSPs to keep for each
                        subject match. (default: None)
  --whitelist TITLE [TITLE ...]
                        subject titles that should be whitelisted (default:
                        None)
  --blacklist TITLE [TITLE ...]
                        subject titles that should be blacklisted (default:
                        None)
  --titleRegex TITLE-REGEX
                        a regex that subject titles must match. (default:
                        None)
  --negativeTitleRegex TITLE-REGEX
                        a regex that subject titles must not match. (default:
                        None)
  --truncateTitlesAfter TRUNCATETITLESAFTER
                        a string that subject titles will be truncated beyond.
                        If the truncated version of a title has already been
                        seen, that title will be skipped. (default: None)
  --minSequenceLen N    subjects of lesser length will be elided. (default:
                        None)
  --maxSequenceLen N    subjects of greater length will be elided. (default:
                        None)
  --taxonomy NAME       the taxonomic group that subjects must match E.g.,
                        "Vira" will filter on viruses. (default: None)
```

</details>

**`alignments-per-read.py`** — Read in BLAST records and print, for each read, a count of how many
sequences it matches.

**`get-reads.py`** — Print the reads that fall in a specified X-axis (read offset) and
bit-score range of a BLAST alignment panel.

<details><summary><strong><code>alignment-panel-civ.py</code></strong> — Given BLAST/DIAMOND JSON output, the corresponding FASTA/FASTQ query
files, and filtering criteria, produce a summary of matched titles and
(optionally) an alignment panel image. The "civ" variants in this group
are tailored to a specific subject-title format
(<code>civ|GENBANK|&lt;protein-id&gt;|GENBANK|&lt;pathogen-accession&gt;|&lt;product&gt;
[&lt;pathogen name&gt;]</code>) used by protein databases built with
<code>make-protein-database.py</code>. <span title="Show --help">ℹ️</span></summary>

```
usage: alignment-panel-civ.py [-h]
                              --proteinGenomeDatabase PROTEINGENOMEDATABASE
                              [--earlyExit] [--matcher {blast,diamond}]
                              --json JSON-file [JSON-file ...]
                              (--fasta FASTA-file [FASTA-file ...] |
                              --fastq FASTQ-file [FASTQ-file ...])
                              [--minStart MINSTART] [--maxStop MAXSTOP]
                              [--oneAlignmentPerRead]
                              [--maxAlignmentsPerRead MAXALIGNMENTSPERREAD]
                              [--scoreCutoff SCORECUTOFF]
                              [--percentageIdenticalCutoff PERCENTAGEIDENTICALCUTOFF]
                              [--percentagePositiveCutoff PERCENTAGEPOSITIVECUTOFF]
                              [--maxHspsPerHit MAXHSPSPERHIT]
                              [--whitelist WHITELIST [WHITELIST ...]]
                              [--whitelistFile WHITELISTFILE]
                              [--blacklist BLACKLIST [BLACKLIST ...]]
                              [--blacklistFile BLACKLISTFILE]
                              [--titleRegex TITLEREGEX]
                              [--negativeTitleRegex NEGATIVETITLEREGEX]
                              [--truncateTitlesAfter TRUNCATETITLESAFTER]
                              [--minSequenceLen MINSEQUENCELEN]
                              [--maxSequenceLen MAXSEQUENCELEN]
                              [--taxonomy TAXONOMY]
                              [--minMatchingReads MINMATCHINGREADS]
                              [--maxMatchingReads MAXMATCHINGREADS]
                              [--minMedianScore MINMEDIANSCORE]
                              [--withScoreBetterThan WITHSCOREBETTERTHAN]
                              [--minNewReads MINNEWREADS]
                              [--maxTitles MAXTITLES]
                              [--minCoverage MINCOVERAGE]
                              [--sortOn {maxScore,medianScore,readCount,length,title}]
                              [--rankValues] --outputDir OUTPUTDIR
                              [--color COLOR] [--equalizeXAxes]
                              [--xRange {reads,subject}] [--logLinearXAxis]
                              [--logBase LOGBASE] [--showFeatures]
                              [--sortFilenames]
                              [--titlesJSONFile TITLESJSONFILE]
                              [--maxBluePlotReads MAXBLUEPLOTREADS]

Produce an alignment panel in HTML

options:
  -h, --help            show this help message and exit
  --proteinGenomeDatabase PROTEINGENOMEDATABASE
                        The filename of an Sqlite3 database holding protein
                        and genome information, as built by make-protein-
                        database.py
  --earlyExit           If True, just print the number of interesting matches,
                        but do not create the alignment panel. (default:
                        False)
  --matcher {blast,diamond}
                        The matching algorithm that was used to produce the
                        JSON. (default: blast)
  --json JSON-file [JSON-file ...]
                        the JSON file(s) of BLAST or DIAMOND output.
  --fasta FASTA-file [FASTA-file ...]
                        the FASTA file(s) of read sequences that were given to
                        BLAST or DIAMOND. (default: None)
  --fastq FASTQ-file [FASTQ-file ...]
                        the FASTQ file(s) of read sequences that were given to
                        BLAST or DIAMOND. (default: None)
  --minStart MINSTART   Reads that start before this subject offset should not
                        be shown. (default: None)
  --maxStop MAXSTOP     Reads that end after this subject offset should not be
                        shown. (default: None)
  --oneAlignmentPerRead
                        If True, only keep the best alignment for each read.
                        (default: False)
  --maxAlignmentsPerRead MAXALIGNMENTSPERREAD
                        Reads with more than this many alignments will be
                        elided. Pass zero to only keep reads with no matches
                        (alignments). (default: None)
  --scoreCutoff SCORECUTOFF
                        A float score. Matches with scores worse than this
                        will be ignored. (default: None)
  --percentageIdenticalCutoff PERCENTAGEIDENTICALCUTOFF
                        A float percentage identity (0.0 to 100.0). Matches
                        with percent identity scores less than this will be
                        ignored. (default: None)
  --percentagePositiveCutoff PERCENTAGEPOSITIVECUTOFF
                        A float percentage (0.0 to 100.0). Matches with
                        percent positive scores less than this will be
                        ignored. A single AA match is considered positive (by
                        DIAMOND) if its entry in the substitution (scoring)
                        matrix is positive. DIAMOND uses the BLOSUM62 matrix
                        by default. (default: None)
  --maxHspsPerHit MAXHSPSPERHIT
                        A numeric max number of HSPs to show for each hit on
                        hitId. (default: None)
  --whitelist WHITELIST [WHITELIST ...]
                        sequence titles that should be whitelisted (default:
                        None)
  --whitelistFile WHITELISTFILE
                        A file containing sequence titles (one per line) to be
                        whitelisted (default: None)
  --blacklist BLACKLIST [BLACKLIST ...]
                        sequence titles that should be blacklisted (default:
                        None)
  --blacklistFile BLACKLISTFILE
                        A file containing sequence titles (one per line) to be
                        blacklisted (default: None)
  --titleRegex TITLEREGEX
                        a regex that sequence titles must match. (default:
                        None)
  --negativeTitleRegex NEGATIVETITLEREGEX
                        a regex that sequence titles must not match. (default:
                        None)
  --truncateTitlesAfter TRUNCATETITLESAFTER
                        a string that titles will be truncated beyond. If the
                        truncated version of a title has already been seen,
                        that title will be skipped. (default: None)
  --minSequenceLen MINSEQUENCELEN
                        sequences of lesser length will be elided. (default:
                        None)
  --maxSequenceLen MAXSEQUENCELEN
                        sequences of greater length will be elided. (default:
                        None)
  --taxonomy TAXONOMY   a string of the taxonomic group on which should be
                        filtered. eg "Vira" will filter on viruses. (default:
                        None)
  --minMatchingReads MINMATCHINGREADS
                        sequences that are matched by fewer reads will be
                        elided. (default: None)
  --maxMatchingReads MAXMATCHINGREADS
                        sequences that are matched by more reads will be
                        elided. (default: None)
  --minMedianScore MINMEDIANSCORE
                        sequences that are matched with a median score that is
                        worse will be elided. (default: None)
  --withScoreBetterThan WITHSCOREBETTERTHAN
                        sequences that are matched without at least one score
                        at least this good will be elided. (default: None)
  --minNewReads MINNEWREADS
                        The fraction of its reads by which a new read set must
                        differ from all previously seen read sets in order to
                        be considered acceptably different. (default: None)
  --maxTitles MAXTITLES
                        The maximum number of titles to keep. If more titles
                        than this result from the filtering, titles will be
                        sorted (according to the --sortOn value) and only the
                        best will be retained. (default: None)
  --minCoverage MINCOVERAGE
                        The (0.0 to 1.0) minimum fraction of a subject
                        sequence that must be matched by at least one read.
                        (default: None)
  --sortOn {maxScore,medianScore,readCount,length,title}
                        The attribute to sort subplots on. (default: maxScore)
  --rankValues          If True, display reads with a Y axis coord that is the
                        rank of the score. (default: False)
  --outputDir OUTPUTDIR
                        Specifies a directory to write the HTML summary to.
                        Required.
  --color COLOR         a string which has a color as the first element and
                        read ids and/or FASTA file names as the following
                        element(s), separated by spaces. (default: None)
  --equalizeXAxes       If True, all alignment graphs will have their X axes
                        drawn with the same range. (default: False)
  --xRange {reads,subject}
                        Set the X axis range to show either the subject or the
                        extent of the reads that hit the subject. (default:
                        subject)
  --logLinearXAxis      If True, convert read offsets so that empty regions in
                        the alignment panel plots will only be as wide as
                        their logged actual values (default: False)
  --logBase LOGBASE     The base of the logarithm to use if logLinearXAxis is
                        True (default: 1.1)
  --showFeatures        If specified, look up features for the individual
                        images in the alignment panel. (default: False)
  --sortFilenames       If specified, the JSON and FASTA/Q file names will be
                        processed in sorted order. The sorting is based on
                        finding a numerical prefix in the filename. This can
                        be useful when processing output files produced by
                        systems like HTCondor, which makes files with names
                        like 1.out, 10.out, etc. that do not sort properly and
                        so cannot conveniently be given to this program along
                        with a single FASTA/Q file (because the order of the
                        results in the files from HTCondor does not match the
                        order of sequences in the FASTA/Q file. (default:
                        False)
  --titlesJSONFile TITLESJSONFILE
                        Give a file name for JSON holding information about
                        titles to be written to (default: None)
  --maxBluePlotReads MAXBLUEPLOTREADS
                        The maximum number of reads to display in a blue plot.
                        Note that because addQueryLines is unconditionally
                        True (in the calls to alignmentGraph in
                        alignmentPanelHTML in src/dark/graphics.py), the
                        number of lines drawn in the plots will be twice the
                        value specified by --maxBluePlotReads (one line for
                        the gray 'whiskers' drawn first for the non-matching
                        region, then one line for the famous blue-plot blue
                        line that shows the actually matched part of each
                        read). (default: None)

Given BLAST or DIAMOND JSON output files, the corresponding FASTA (or FASTQ)
sequence files, and filtering criteria, produce an alignment panel.
```

</details>

<details><summary><strong><code>proteins-to-pathogens-civ.py</code></strong> — Group protein match output (as produced by the alignment-panel scripts
above) by the pathogen each protein comes from, using the pathogen
accession number embedded in the "civ"-style subject title. <span title="Show --help">ℹ️</span></summary>

```
usage: proteins-to-pathogens-civ.py [-h]
                                    --proteinGenomeDatabase PROTEINGENOMEDATABASE
                                    [--sampleName SAMPLENAME |
                                    --sampleNameRegex SAMPLENAMEREGEX]
                                    [--pathogenPanelFilename [PATHOGENPANELFILENAME]]
                                    [--sampleIndexFilename SAMPLEINDEXFILENAME]
                                    [--html] [--format {fasta,fastq}]
                                    [--minProteinFraction MINPROTEINFRACTION]
                                    [--minProteinCount MINPROTEINCOUNT]
                                    [--pathogenType {bacterial,viral,generic}]
                                    [--showReadLengths] [--assetDir ASSETDIR]
                                    [--pathogenDataDir PATHOGENDATADIR]
                                    [--title TITLE] [--preamble PREAMBLE]
                                    [--titleRegex TITLEREGEX]
                                    [--negativeTitleRegex NEGATIVETITLEREGEX]
                                    [--omitVirusLinks]
                                    [--defaultReadCountColor DEFAULTREADCOUNTCOLOR]
                                    [--readCountColor READCOUNTCOLOR]
                                    [--bootstrapTreeviewDir BOOTSTRAPTREEVIEWDIR]
                                    [--taxonomyDatabase DATABASE-FILE]
                                    [filenames ...]

Group proteins by the pathogen they're from.

positional arguments:
  filenames             File names to read input from. The input will
                        typically be generated by alignment-panel-civ.py
                        (default: None)

options:
  -h, --help            show this help message and exit
  --proteinGenomeDatabase PROTEINGENOMEDATABASE
                        The filename of an Sqlite3 database holding protein
                        and genome information, as built by make-protein-
                        database.py
  --sampleName SAMPLENAME
                        An (optional) sample name. Use when all input files
                        are for a single sample. Cannot be used with
                        --sampleNameRegex. (default: None)
  --sampleNameRegex SAMPLENAMEREGEX
                        An (optional) regular expression that can be used to
                        extract a short sample name from full sample file
                        name. The regular expression must have a matching
                        group (delimited by parentheses) that captures the
                        part of the file name that should be used as the
                        sample name. (default: None)
  --pathogenPanelFilename [PATHOGENPANELFILENAME]
                        An (optional) filename to write a pathogen-sample
                        panel PNG image to. (default: None)
  --sampleIndexFilename SAMPLEINDEXFILENAME
                        An (optional) filename to write a sample index file
                        to. Lines in the file will have an integer index, a
                        space, and then the sample name. Only produced if
                        --html is used (because the pathogen-NNN-sample-
                        MMM.fastq are only written in that case). (default:
                        None)
  --html                If specified, output HTML instead of plain text.
                        (default: False)
  --format {fasta,fastq}
                        Give the format of the sequence files written by
                        alignment-panel-civ.py. (default: fasta)
  --minProteinFraction MINPROTEINFRACTION
                        The minimum fraction of proteins in a pathogen that
                        must be matched by a particular sample in order for
                        that pathogen to be displayed for that sample.
                        (default: 0.0)
  --minProteinCount MINPROTEINCOUNT
                        The minimum number of proteins in a pathogen that must
                        be matched by a particular sample in order for that
                        pathogen to be displayed for that sample. (default: 0)
  --pathogenType {bacterial,viral,generic}
                        Specify the pathogen type. This option only affects
                        the language used in HTML output. (default: viral)
  --showReadLengths     If specified, the HTML output (use --html to get this)
                        will contain the lengths of all reads that match
                        proteins for a pathogen. (default: False)
  --assetDir ASSETDIR   The output directory where noninteractive-alignment-
                        panel.py put its plots and FASTA/FASTQ files.
                        (default: out)
  --pathogenDataDir PATHOGENDATADIR
                        The directory where per-pathogen information (e.g.,
                        collected reads across all samples) should be written.
                        (default: pathogen-data)
  --title TITLE         The title to show at the top of the output. (default:
                        None)
  --preamble PREAMBLE   Optional preamble text to show after the title. The
                        argument value may also name a file, in which case the
                        file contents will be inserted into the output. May be
                        repeated. (default: None)
  --titleRegex TITLEREGEX
                        A regex that pathogen names must match. (default:
                        None)
  --negativeTitleRegex NEGATIVETITLEREGEX
                        a regex that pathogen names must not match. (default:
                        None)
  --omitVirusLinks      If specified, the HTML output (use --html to get this)
                        for viruses will not contain links to ICTV and
                        ViralZone. This should be used when working with
                        viruses that do not yet have names that can be looked
                        up. (default: False)
  --defaultReadCountColor DEFAULTREADCOUNTCOLOR
                        The font color for read counts. This will be used for
                        all read counts that do not otherwise have a color due
                        to use of --readCountColor. Only valid if --html is
                        used. (default: black)
  --readCountColor READCOUNTCOLOR
                        Specify read count coloring. This option must be given
                        as a space separated "value color" pair. The value is
                        an integer read count and the color is any color
                        specification that can be given to CSS. This argument
                        can be repeated. E.g., --readCountColor "0.9 red"
                        --readCountColor "0.75 rgb(23, 190, 207)"
                        --readCountColor "0.1 #CF3CF3". Read counts will be
                        colored using the color of the highest count threshold
                        they satisfy. The default is to color all read counts
                        with the --defaultReadCountColor color. Only valid if
                        --html is used. (default: None)
  --bootstrapTreeviewDir BOOTSTRAPTREEVIEWDIR
                        The directory where the bootstrap treeview JS and CSS
                        files can be found. This can be a relative path from
                        the location where the output HTML will be served from
                        or can be a URL directory. In both cases, HTML will be
                        emitted looking for files with names bootstrap-
                        treeview.min.css and bootstrap-treeview.min.js in the
                        directory (or URL directory). Only valid if --html is
                        used. (default: None)
  --taxonomyDatabase DATABASE-FILE
                        The file holding an sqlite3 taxonomy database. See
                        https://github.com/acorg/ncbi-taxonomy-database for
                        how to build one. If not specified, the value in the
                        'DARK_MATTER_TAXONOMY_DATABASE' environment variable
                        (if any) will be used. (default: None)
```

</details>

## 8. Protein/genome SQLite databases and GenBank data

<details><summary><strong><code>make-fasta-database.py</code></strong> — Create an sqlite3 database from FASTA sequences. <span title="Show --help">ℹ️</span></summary>

```
usage: make-fasta-database.py [-h] --out OUT [--force] [--quiet]
                              --fasta FASTA-file [FASTA-file ...]

Create an sqlite3 database from FASTA sequences.

options:
  -h, --help            show this help message and exit
  --out OUT             The output file. This file must not exist (use --force
                        to overwrite).
  --force               If True and the output file already exists, overwrite
                        it.
  --quiet               If True do not print indexing progress.
  --fasta FASTA-file [FASTA-file ...]
                        the FASTA file(s) to make the database from. These may
                        be uncompressed, or compressed with bgzip (from
                        samtools), with a .gz suffix.
```

</details>

<details><summary><strong><code>make-protein-database.py</code></strong> — Create a genome/protein sqlite3 database from GenBank and JSON files;
protein sequences for the database are printed to stdout for indexing by
other tools. (This is the database consumed by the <code>-civ</code> scripts in
category 7.) <span title="Show --help">ℹ️</span></summary>

```
usage: make-protein-database.py [-h] --databaseFile DATABASEFILE
                                [--databaseName DATABASENAME]
                                [--duplicationPolicy {error,ignore}]
                                [--progress] [--logFile LOGFILE]
                                [--noWarnings] [--dnaOnly | --rnaOnly]
                                [--maxGenomeLength MAXGENOMELENGTH]
                                [--minGenomeLength MINGENOMELENGTH]
                                [--excludeFungusOnlyViruses]
                                [--excludePlantOnlyViruses]
                                [--gb GenBank-file [GenBank-file ...]]
                                [--json JSON-file [JSON-file ...]]
                                [--allowedTaxonomicRanks name:rank [name:rank ...]]
                                [--proteinSource PROTEINSOURCE]
                                [--genomeSource GENOMESOURCE]
                                [--excludeExclusiveHost hostname [hostname ...]]
                                [--taxonomyDatabase DATABASE-FILE]

Create a genome/protein sqlite3 database from GenBank and JSON files. The
protein sequences for the database are printed to standard output for indexing
by other tools.

options:
  -h, --help            show this help message and exit
  --databaseFile DATABASEFILE
                        The output file. This file must not exist (use --force
                        to overwrite).
  --databaseName DATABASENAME
                        The database that the records in the (--gb) GenBank
                        files came from (e.g., "refseq" or "RVDB"). (default:
                        None)
  --duplicationPolicy {error,ignore}
                        What to do if a to-be-inserted accession number is
                        already present in the database. "error" results in a
                        ValueError being raised, "ignore" will ignore the
                        duplicate. It should also be possible to update (i.e.,
                        replace) but that is not supported yet. (default:
                        ignore)
  --progress            Print indexing progress. (default: False)
  --logFile LOGFILE     Write indexing details to a log file. (default: None)
  --noWarnings          Do not print warnings about unparseable GenBank or
                        JSON records. (default: False)
  --dnaOnly             If given, only include DNA viruses. (default: False)
  --rnaOnly             If given, only include RNA viruses. (default: False)
  --maxGenomeLength MAXGENOMELENGTH
                        Genomes longer than this will not be considered.
                        (default: None)
  --minGenomeLength MINGENOMELENGTH
                        Genomes shorter than this will not be considered.
                        (default: None)
  --excludeFungusOnlyViruses
                        If given, exclude fungus-only viruses (i.e., viruses
                        that only infect fungi host species). (default: False)
  --excludePlantOnlyViruses
                        If given, exclude plant-only viruses (i.e., viruses
                        that only infect plant host species). (default: False)
  --gb GenBank-file [GenBank-file ...]
                        The GenBank file(s) to make the database from. These
                        may be uncompressed, or compressed with bgzip (from
                        samtools), with a .gz suffix. (default: None)
  --json JSON-file [JSON-file ...]
                        The JSON file(s) to make the database from. These
                        contain genome and protein information for cases where
                        sequences that are not in GenBank should be added.
                        (default: None)
  --allowedTaxonomicRanks name:rank [name:rank ...]
                        Strings of the form name:rank to restrict included
                        viruses to be from any of a specific set of taxonomic
                        ranks. May be repeated. E.g., --allowedTaxonomicRanks
                        nidovirales:order --allowedTaxonomicRanks
                        retroviridae:family Viruses will be included only if
                        they match at least one of the name, rank pairs. The
                        strings are case insensitive. (default: None)
  --proteinSource PROTEINSOURCE
                        The source of the accession numbers for the proteins
                        found in the input files. This becomes part of the
                        sequence id printed in the protein FASTA output.
                        (default: GenBank)
  --genomeSource GENOMESOURCE
                        The source of the accession numbers for the genomes in
                        the input files. This becomes part of the sequence id
                        printed in the protein FASTA output. (default:
                        GenBank)
  --excludeExclusiveHost hostname [hostname ...]
                        A host type that should be excluded, but only if this
                        is the only host of the virus. (default: None)
  --taxonomyDatabase DATABASE-FILE
                        The file holding an sqlite3 taxonomy database. See
                        https://github.com/acorg/ncbi-taxonomy-database for
                        how to build one. If not specified, the value in the
                        'DARK_MATTER_TAXONOMY_DATABASE' environment variable
                        (if any) will be used. (default: None)
```

</details>

<details><summary><strong><code>describe-protein-database.py</code></strong> — Print a description of a protein/genome database. <span title="Show --help">ℹ️</span></summary>

```
usage: describe-protein-database.py [-h] databaseFile

Print a description of a protein/genome database.

positional arguments:
  databaseFile  The Sqlite3 database file.

options:
  -h, --help    show this help message and exit
```

</details>

<details><summary><strong><code>subset-protein-database.py</code></strong> — Query an sqlite3 protein/genome database and write the matching sequences
to stdout as FASTA. <span title="Show --help">ℹ️</span></summary>

```
usage: subset-protein-database.py [-h] --databaseFile DATABASEFILE
                                  [--query QUERY] [--queryFile QUERYFILE]
                                  [--keepDuplicateIds] [--idsOnly]

Query an sqlite3 protein/genome database and write the matching sequences to
standard output as FASTA.

options:
  -h, --help            show this help message and exit
  --databaseFile DATABASEFILE
                        The protein/genome database file. This must be an
                        sqlite3 database file created by make-protein-
                        database.py
  --query QUERY         An SQL database query. The query must result in rows
                        with two columns, the sequence id and the sequence.
                        These can be protein or nucleotide sequences, but it
                        is probably better not to mix them! (default: None)
  --queryFile QUERYFILE
                        The name of a file containing an SQL database query.
                        The query must result in rows with two columns, the
                        sequence id and the sequence. These can be protein or
                        nucleotide sequences, but it is probably better not to
                        mix them! (default: None)
  --keepDuplicateIds    If specified, do not remove duplicate ids. (default:
                        False)
  --idsOnly             If specified, only print matching sequence ids.
                        (default: False)
```

</details>

<details><summary><strong><code>parse-genbank-flat-file.py</code></strong> — Parse GenBank flat files; exits non-zero if parsing fails (also used as a
sanity check by <code>download-genbank.sh</code>). <span title="Show --help">ℹ️</span></summary>

```
usage: parse-genbank-flat-file.py [-h] [--expectedCount EXPECTEDCOUNT]
                                  [--quiet]
                                  files [files ...]

Parse GenBank flat files. Exit non-zero if parsing fails.

positional arguments:
  files                 The files to parse.

options:
  -h, --help            show this help message and exit
  --expectedCount EXPECTEDCOUNT
                        The expected number of records in each GenBank file.
                        If given and this number is not present in any file,
                        exit non-zero. (default: None)
  --quiet               If given, write no output. (default: False)
```

</details>

<details><summary><strong><code>genbank-grep.py</code></strong> — Read a GenBank flat file and print GenBank or FASTA for matching records. <span title="Show --help">ℹ️</span></summary>

```
usage: genbank-grep.py [-h] [--idRegex REGEX]
                       [--match {fullmatch,match,search}] [--matchDescription]
                       [--invert] [--ignoreCase] [--format {fasta,gb}]
                       [--expectedCount EXPECTEDCOUNT] [--verbose]
                       FILE.gb [FILE.gb ...]

Read a GenBank flat file and print GenBank or FASTA for matching records.

positional arguments:
  FILE.gb               The GenBank file(s) to examine.

options:
  -h, --help            show this help message and exit
  --idRegex REGEX       Only print records whose id (or description, if
                        --matchDescription is used) matches this regular
                        expression. (default: .)
  --match {fullmatch,match,search}
                        How to match the regular expression. See pydoc re for
                        an explanation. (default: search)
  --matchDescription    Match on the record description as well as the id.
                        (default: False)
  --invert, -v          Invert the match. Only print records that do not
                        match. (default: False)
  --ignoreCase, -i      Use case-insensitive matching. (default: False)
  --format {fasta,gb}   The output format. (default: gb)
  --expectedCount EXPECTEDCOUNT
                        The expected number of records to find. If this number
                        is not found, print a message and exit non-zero.
                        (default: None)
  --verbose             Print a summary of records read and printed. (default:
                        False)
```

</details>

**`download-genbank.sh`** — Download GenBank flat files for a list of accession numbers from NCBI, in
parallel batches, with retries. One of a documented set of tools (see
`doc/protein-database.md`) for building a protein database.

**`download-refseq-viral-fasta.sh`** — Download the current RefSeq viral FASTA release from NCBI's FTP site
(see also `download-refseq-viral-gbff.sh`).

**`download-refseq-viral-gbff.sh`** — Download the current RefSeq viral GenBank flat-file release from NCBI's
FTP site (see also `download-refseq-viral-fasta.sh`).

<details><summary><strong><code>ncbi-fetch-id.py</code></strong> — Fetch a sequence by id from NCBI and write it to stdout as FASTA. <span title="Show --help">ℹ️</span></summary>

```
usage: ncbi-fetch-id.py [-h] [--database {nucleotide,protein}]
                        [--format {fasta,gb}]
                        id

Fetch a sequence by id from NCBI and write it to stdout as FASTA.

positional arguments:
  id                    The id of the sequence to fetch.

options:
  -h, --help            show this help message and exit
  --database {nucleotide,protein}
                        The name of the NCBI database to query. (default:
                        nucleotide)
  --format {fasta,gb}   The format to return. (default: fasta)
```

</details>

## 9. Taxonomy

<details><summary><strong><code>get-taxonomy.py</code></strong> — Print taxonomy information (e.g., lineage) for given accession numbers or
taxonomy ids, using a local taxonomy database. <span title="Show --help">ℹ️</span></summary>

```
usage: get-taxonomy.py [-h] [--separator SEPARATOR] [--namesOnly] [--print]
                       [--taxonomyDatabase DATABASE-FILE]
                       [id ...]

Print taxonomy information for accession numbers or taxonomy ids.

positional arguments:
  id                    The ids (accession numbers, names, or taxonomy ids) to
                        print taxonomy information for. If not given, ids are
                        read from standard input, one per line. (default:
                        None)

options:
  -h, --help            show this help message and exit
  --separator SEPARATOR
                        The string to separate output columns with. Default is
                        to print aligned output padded with multiple spaces
                        (or a TAB if --namesOnly is used). (default: None)
  --namesOnly           If specified, only print the taxonomic names.
                        (default: False)
  --print               If specified, also print the id. (default: False)
  --taxonomyDatabase DATABASE-FILE
                        The file holding an sqlite3 taxonomy database. See
                        https://github.com/acorg/ncbi-taxonomy-database for
                        how to build one. If not specified, the value in the
                        'DARK_MATTER_TAXONOMY_DATABASE' environment variable
                        (if any) will be used. (default: None)
```

</details>

<details><summary><strong><code>get-hosts.py</code></strong> — Print known hosts for given accession numbers or taxonomy ids, using a
local taxonomy database. <span title="Show --help">ℹ️</span></summary>

```
usage: get-hosts.py [-h] --database DATABASE [--printId]
                    [--taxonomyDatabase DATABASE-FILE]
                    [ids ...]

Print hosts for accession numbers or taxonomy ids

positional arguments:
  ids                   The ids (accession numbers, names, or taxonomy ids) to
                        print host information for. If not given, ids are read
                        from standard input, one per line. (default: None)

options:
  -h, --help            show this help message and exit
  --database DATABASE   The file holding the sqlite3 taxonomy database. See
                        https://github.com/acorg/ncbi-taxonomy-database for
                        how to build one.
  --printId             If specified, also print the id. (default: False)
  --taxonomyDatabase DATABASE-FILE
                        The file holding an sqlite3 taxonomy database. See
                        https://github.com/acorg/ncbi-taxonomy-database for
                        how to build one. If not specified, the value in the
                        'DARK_MATTER_TAXONOMY_DATABASE' environment variable
                        (if any) will be used. (default: None)
```

</details>

## 10. Phylogenetic trees (Newick)

<details><summary><strong><code>curate-trees-dendropy.py</code></strong> — Curate phylogenetic trees (i.e., their specification/newick files) using <code>dendropy</code>. <span title="Show --help">ℹ️</span></summary>

```
usage: curate-trees-dendropy.py [-h] [--treeFiles TREEFILES [TREEFILES ...]]
                                [--outSuffix OUTSUFFIX] [--dryRun] [--force]
                                [--rootNode ROOTNODE] [--collapse COLLAPSE]
                                [--discardLengths]
                                [--ladderize {ascending,descending}]
                                [--format {fasta,phylip,newick,nexml,nexus}]
                                [--scale SCALE]

Curate phylogenetic trees (i.e., their specification files).

options:
  -h, --help            show this help message and exit
  --treeFiles TREEFILES [TREEFILES ...]
                        The tree files to curate. (default: None)
  --outSuffix OUTSUFFIX
                        The suffix to use on curated output files. (default:
                        .treefile-curated)
  --dryRun              Just print what would be done. (default: False)
  --force               Overwrite pre-existing result files. (default: False)
  --rootNode ROOTNODE   Re-root the tree on the incoming branch to the tip
                        node with this name. (default: None)
  --collapse COLLAPSE   A floating-point support value. Branches with support
                        less than this value will be collapsed (to form
                        polytomies). The value to use will depend on what the
                        program used to create the tree used to indicate the
                        support. This is frequently the number of times the
                        branch was supported in the bootstrap analysis. In
                        this case, if you did 100 bootstraps, you would see
                        integer support values from 0 to 100. If you did 1000
                        bootstraps, the values could go up to 1000. Other
                        programs might write the support value as a fraction
                        (of bootstraps in which the branch is supported),
                        therefore ranging from 0.0 to 1.0. In other words, you
                        need to know what's in your tree file to know for sure
                        what value to use for this option. More difficult is
                        to know what cut-off value actually makes sense in
                        your analysis. If you're using iqtree2 with its
                        ultrafast bootstrap (UFBoot), you should read the
                        related FAQ item at
                        http://www.iqtree.org/doc/Frequently-Asked-Questions
                        (default: None)
  --discardLengths      When collapsing edges due to low support, the length
                        of the edge being removed will be added to the length
                        of the edges leading to the nodes of the edge-to-be-
                        removed's child edges. This may be hard to follow, in
                        which case you should draw yourself a little picture
                        of part of a rooted tree and think about what happens
                        to the nodes that become detached when you remove an
                        edge (and the node that it leads to in the direction
                        away from the root). If the length of the edge you are
                        collapsing is important, because it came from a tree-
                        making algorithm that sets branch lengths according to
                        sequence change, the 'right' thing to do is to add the
                        length of the to-be-deleted edge to the lengths of its
                        descendant edges when you attach those nodes as
                        children of the node at the root-end of the edge you
                        are deleting. If for some reason you also want to
                        discard the length of the edge you are collapsing, you
                        can use this --discardLengths option. Make sure you
                        really understand the implications of doing so,
                        however, because discarding that length will cause the
                        patristic distance of these child nodes to be shorter
                        than it was estimated to be by the algorithm that made
                        the tree and set the original branch lengths. And why
                        would you want to second-guess the algorithm and mess
                        with its conclusions? (default: False)
  --ladderize {ascending,descending}
                        How to (optionally) ladderize the tree after re-
                        rooting. (default: None)
  --format {fasta,phylip,newick,nexml,nexus}
                        The output format. (default: newick)
  --scale SCALE         Factor to scale the length of the branches leading to
                        the descendants of the new root (if --rootNode is
                        used). This can be used to reduce the length of very
                        long branches that result from rerooting on a distant
                        outgroup. If you use this option for an image in a
                        publication, you of course must point out that you
                        have done so! (default: 1.0)
```

</details>

<details><summary><strong><code>curate-tree-ete.py</code></strong> — The same kind of tree curation, but using the <code>ete4</code> library instead of
<code>dendropy</code> (these look like two alternative implementations — worth
deciding whether both are needed going forward, especially given the
recent move to <code>ete4</code> noted in the CHANGELOG). <span title="Show --help">ℹ️</span></summary>

```
usage: curate-tree-ete.py [-h] [--parser PARSER] [--rootNode ROOTNODE]
                          [--scale SCALE] [--detach] [--collapse]
                          [--supportThreshold SUPPORTTHRESHOLD]
                          [--alrtThreshold ALRTTHRESHOLD] [--removeSupport]
                          [--removeAlrt] [--verbose]
                          treeFile

Curate a phylogenetic tree using the ete4 library and print the result.

positional arguments:
  treeFile              The tree file to curate.

options:
  -h, --help            show this help message and exit
  --parser, --format PARSER
                        The Newick tree flavor. See https://etetoolkit.github.
                        io/ete/reference/reference_parsers.html for
                        information on the Newick variations ete4 understands.
                        The output tree will be in the same format as the
                        input. (default: 0)
  --rootNode ROOTNODE   Re-root the tree on the branch entering the node that
                        is the most-recent common ancestor of the tips with
                        these names (may be repeated). (default: None)
  --scale SCALE         Factor to scale (i.e., divide) the length of the
                        branches leading to the descendants of the new root
                        (if --rootNode is used). This can be used to reduce
                        the length of very long branches that result from
                        rerooting on a distant outgroup. If you use this
                        option for an image in a publication, you of course
                        must point out that you have done so! (default: 1.0)
  --detach              Detach the named root node from the tree after rooting
                        on the edge coming into it. The newly added (by
                        rooting) parent of the named node will be the
                        (temporary) new root node. This option will first
                        remove the named node as a child of that root, leaving
                        it with just one child. The remaining child will then
                        be made the new root. This essentially like rooting on
                        an outgroup taxon then completely removing the
                        outgroup from the tree. (default: False)
  --collapse            Collapse branches that are insufficiently supported.
                        (default: False)
  --supportThreshold SUPPORTTHRESHOLD
                        The support level below which edges will be collapsed
                        (only valid if --collapse is used). (default: None)
  --alrtThreshold ALRTTHRESHOLD
                        The alrt level below which edges will be collapsed
                        (only valid if --collapse is used). Alrt values are
                        added by iqtree2 if you use its -alrt option and can
                        be used to check for well-supported branches.
                        (default: None)
  --removeSupport       Remove support values from internal nodes. (default:
                        False)
  --removeAlrt          Remove alrt values from internal nodes. (default:
                        False)
  --verbose             Print information about processing to standard error.
                        (default: False)
```

</details>

<details><summary><strong><code>newick-to-ascii.py</code></strong> — Print a Newick phylogenetic tree in ASCII. <span title="Show --help">ℹ️</span></summary>

```
usage: newick-to-ascii.py [-h] [--treeFile FILENAME]
                          [--outgroupRegex TAXON-REGEX] [--parser N]
                          [--verbose]

Print a Newick phylogenetic tree in ASCII.

options:
  -h, --help            show this help message and exit
  --treeFile FILENAME   The Newick tree file. (default: <_io.TextIOWrapper
                        name='<stdin>' mode='r' encoding='utf-8'>)
  --outgroupRegex TAXON-REGEX
                        A regular expression (for taxa names). The lowest
                        common parent of any matching taxa will be used as an
                        outgroup. (default: None)
  --parser, --format N  The ete4 Newick parser to use. (default: None)
  --verbose             Print information about the outgroup (if any) taxa to
                        standard error. (default: False)
```

</details>

<details><summary><strong><code>relabel-newick-tree.py</code></strong> — Relabel a Newick tree according to a file of from/to names. <span title="Show --help">ℹ️</span></summary>

```
usage: relabel-newick-tree.py [-h] [--nameFile NAMEFILE]
                              [--inputFormat INPUTFORMAT]
                              [--outputFormat OUTPUTFORMAT] [--quotedNames]
                              [--verbose]
                              newickFile

Relabel a Newick tree according to a file of from/to names.

positional arguments:
  newickFile            The Newick tree file.

options:
  -h, --help            show this help message and exit
  --nameFile NAMEFILE   The file of renamings. Each line should contain an
                        existing (in the tree) sequence id, a TAB, and then a
                        new sequence id. If not given, the renaming is read
                        from standard input. (default: None)
  --inputFormat INPUTFORMAT
                        The parser argument to pass to the ete Tree class
                        initializer to read the input Newick. (default: 0)
  --outputFormat OUTPUTFORMAT
                        The parser argument to pass to the ete Tree writer for
                        the resulting tree. (default: 0)
  --quotedNames         If given, pass quoted_node_names=True to the ete Tree
                        class initializer (default: False)
  --verbose             Print details of the renaming(s). (default: False)
```

</details>

<details><summary><strong><code>create-newick-relabeling-output.py</code></strong> — Find sequences with a specific property (via a user-supplied "relabel"
function) and print a TAB-separated list of renamings, suitable as input
to <code>relabel-newick-tree.py</code>. <span title="Show --help">ℹ️</span></summary>

```
usage: create-newick-relabeling-output.py [-h]
                                          --relabelFunctionFile RELABELFUNCTIONFILE
                                          [--fastaFile FILENAME]
                                          [--readClass CLASSNAME] [--fasta |
                                          --fastq | --fasta-ss]
                                          [FILENAME]

Find sequences that have a specific property (as detected by the "relabel"
function in the file specified by --relabelFunctionFile) and print a TAB-
separated list of renamings to standard output.

positional arguments:
  FILENAME              The name of the FASTA/FASTQ input file. Standard input
                        will be read if '-' is used or if no file name is
                        given. (default: None)

options:
  -h, --help            show this help message and exit
  --relabelFunctionFile RELABELFUNCTIONFILE
                        A file containing a Python function named "relabel"
                        that takes a list of reads and returns a dictionary of
                        renamings.
  --fastaFile FILENAME  The name of the FASTA/FASTQ input file. Kept for
                        backward compatibility - giving the file name as a
                        plain argument (with no option) does the same thing.
                        (default: None)
  --readClass CLASSNAME
                        If specified, give the type of the reads in the input.
                        Possible choices: AARead, AAReadORF, AAReadWithX,
                        DNARead, RNARead, Read, SSAARead, SSAAReadWithX,
                        TranslatedRead. (default: DNARead)
  --fasta               If specified, input will be treated as FASTA. If none
                        of --fasta, --fastq, or --fasta-ss is given, the input
                        format is detected automatically by looking at its
                        first non-whitespace character ('>' for FASTA, '@' for
                        FASTQ). (default: False)
  --fastq               If specified, input will be treated as FASTQ.
                        (default: False)
  --fasta-ss            If specified, input will be treated as PDB FASTA
                        (i.e., regular FASTA with each sequence followed by
                        its structure). This cannot be auto-detected, so must
                        be given explicitly. (default: False)
```

</details>

<details><summary><strong><code>remove-alrt-from-tree.py</code></strong> — Read a Newick tree produced by <code>iqtree2</code> with both aLRT and bootstrap
support values on its nodes, drop the aLRT values, optionally adjust the
support values against given thresholds, and print the resulting tree. <span title="Show --help">ℹ️</span></summary>

```
usage: remove-alrt-from-tree.py [-h] [--supportThreshold SUPPORTTHRESHOLD]
                                [--alrtThreshold ALRTTHRESHOLD]
                                [--adjustSupportValues]
                                [--replacementSupportValue REPLACEMENTSUPPORTVALUE]
                                [--expectedUnlabeledNodeCount EXPECTEDUNLABELEDNODECOUNT]
                                [--verbose]
                                treefile

Read a Newick tree file produced by iqtree2 with both alrt and support values
on nodes, drop the alrt values and optionally adjust the support values (based
on the passed bootstrap support threshold and alrt threshold). Print the new
tree as a (Newick) string.

positional arguments:
  treefile              The Newick tree file to read.

options:
  -h, --help            show this help message and exit
  --supportThreshold SUPPORTTHRESHOLD
                        The integer minimum UF boot support value. (default:
                        95)
  --alrtThreshold ALRTTHRESHOLD
                        The float minimum alrt value. (default: 80.0)
  --adjustSupportValues
                        The support on branches that have support values at
                        least --supportThreshold but alrt values below
                        --alrtThreshold will be set to the
                        --replacementSupportValue value. (default: False)
  --replacementSupportValue REPLACEMENTSUPPORTVALUE
                        The string value to use as a label on branches whose
                        support values are adjusted. (default: None)
  --expectedUnlabeledNodeCount EXPECTEDUNLABELEDNODECOUNT
                        The number of unlabeled nodes (in the input file) to
                        allow. (default: None)
  --verbose             Print information about how many branches have good
                        support but too-low alrt values. (default: False)
```

</details>

<details><summary><strong><code>tree-info.py</code></strong> — Print information about phylogenetic tree nodes. <span title="Show --help">ℹ️</span></summary>

```
usage: tree-info.py [-h] [--format {newick,nexml,nexus,phylip}]
                    [--indent INDENT] [--lengths] [--rootNode ROOTNODE]
                    [--noTaxon NOTAXON] [--noLabel NOLABEL]
                    treeFile

Print information phylogenetic tree nodes.

positional arguments:
  treeFile              The tree file to examine.

options:
  -h, --help            show this help message and exit
  --format {newick,nexml,nexus,phylip}
                        The input tree format. (default: newick)
  --indent INDENT       The number of spaces to increase node indentation by
                        when descending to children. (default: 2)
  --lengths             Also show edge lengths. (default: False)
  --rootNode ROOTNODE   Re-root the tree on the edge leading into the tip node
                        with this name. (default: None)
  --noTaxon NOTAXON     What to print if a tip node has no taxon name.
                        (default: NO-TAXON!)
  --noLabel NOLABEL     What to print if an internal node has no label.
                        (default: NO-LABEL!)
```

</details>

## 11. Small, general-purpose utilities

<details><summary><strong><code>ids.py</code></strong> — Print a generated sequence of ids (e.g., <code>id1</code>, <code>id2</code>, ... or a custom
prefix and numeric range), with options for separators and zero-padding.
General-purpose scripting/test-data helper; not specific to biological
sequences. <span title="Show --help">ℹ️</span></summary>

```
usage: ids.py [-h] [--prefix PREFIX] [--sep SEP] [-0] [-n MAXRESULTS]
              ID [ID ...]

Print a generated sequence of ids (e.g., id1, id2, ... or a custom prefix and
numeric range).

positional arguments:
  ID                    One or two arguments to specify the start and
                        (optionally) the end range of ids to create.

options:
  -h, --help            show this help message and exit
  --prefix PREFIX       The id prefix. If this is not given and two arguments
                        are supplied, their longest common prefix will be
                        used.
  --sep SEP             The string to print between successive ids.
  -0, --zeroes          Print numeric suffixes with leading zeroes.
  -n, --count MAXRESULTS
                        The number of items to generate.
```

</details>

<details><summary><strong><code>dark-matter-version.py</code></strong> — Print the installed <code>dark-matter</code> package version. <span title="Show --help">ℹ️</span></summary>

```
usage: dark-matter-version.py [-h]

Print the installed dark-matter package version.

options:
  -h, --help  show this help message and exit
```

</details>

# Deprecated

<details><summary><strong><code>compare-consensuses.py</code></strong> — Compare different consensus-making methods (presumably running several of
the above and diffing/scoring the results). <span title="Show --help">ℹ️</span></summary>

```
usage: compare-consensuses.py [-h] --referenceFile FILENAME
                              --alignmentFile FILENAME [--verbose N] [--force]
                              [--outputDir DIRNAME]

Compare different consensus making methods.

options:
  -h, --help            show this help message and exit
  --referenceFile FILENAME
                        The name of the FASTA file containing the reference
                        sequence.
  --alignmentFile FILENAME
                        The name of the SAM or BAM file containing an
                        alignment to the reference.
  --verbose N           The integer verbosity level (0 = no output, 1 = some
                        output, 2 = maximal output). (default: 0)
  --force               If given, overwrite pre-existing files. (default:
                        False)
  --outputDir DIRNAME   The directory to save result files to. (default: None)
```

</details>

**`convert-sam-to-fastq.sh`** — Read SAM (see http://samtools.github.io/hts-specs/SAMv1.pdf) from stdin
and write FASTQ to stdout. The `+` line of each FASTQ record also repeats
the sequence id, since some FASTQ-consuming tools insist on it being
present even though it's technically optional. (Not previously listed in
this document — added when generating it from doc/scripts.toml.)

<details><summary><strong><code>filter-fasta-by-taxonomy.py</code></strong> — Depends on NCBI's GI numbering scheme, which NCBI retired in 2016, so
this will silently do nothing useful on modern NCBI data unless you still
have old-style FASTA/hit files lying around. (Also listed in category 1,
where it still appears as a FASTA filter.) <span title="Show --help">ℹ️</span></summary>

```
usage: filter-fasta-by-taxonomy.py [-h] --taxonomy TAXONOMY [--invert]
                                   [--detailsFile FILE]

Filter FASTA based on taxonomy

options:
  -h, --help           show this help message and exit
  --taxonomy TAXONOMY  The regex to match the taxonomy on. Case is ignored.
  --invert             If True, only write sequences whose taxonomy does not
                       match.
  --detailsFile FILE   The name of a file to save taxonomy details to

Read DNA FASTA from stdin and print FASTA to stdout, only including sequences
that match the requested taxonomy regular expression at any level.
```

</details>

<details><summary><strong><code>noninteractive-alignment-panel.py</code></strong> — The same alignment-panel functionality as <code>alignment-panel-civ.py</code>,
generated non-interactively (these two scripts look very close in
structure/intent — worth checking whether one has simply superseded the
other). <span title="Show --help">ℹ️</span></summary>

```
usage: noninteractive-alignment-panel.py [-h] [--earlyExit]
                                         [--matcher {blast,diamond}]
                                         --json JSON-file [JSON-file ...]
                                         (--fasta FASTA-file [FASTA-file ...] |
                                         --fastq FASTQ-file [FASTQ-file ...])
                                         [--databaseFastaFilename DATABASEFASTAFILENAME]
                                         [--sqliteDatabaseFilename SQLITEDATABASEFILENAME]
                                         [--databaseFastaDirectory DATABASEFASTADIRECTORY]
                                         [--minStart MINSTART]
                                         [--maxStop MAXSTOP]
                                         [--oneAlignmentPerRead]
                                         [--maxAlignmentsPerRead MAXALIGNMENTSPERREAD]
                                         [--scoreCutoff SCORECUTOFF]
                                         [--percentageIdenticalCutoff PERCENTAGEIDENTICALCUTOFF]
                                         [--percentagePositiveCutoff PERCENTAGEPOSITIVECUTOFF]
                                         [--maxHspsPerHit MAXHSPSPERHIT]
                                         [--whitelist WHITELIST [WHITELIST ...]]
                                         [--whitelistFile WHITELISTFILE]
                                         [--blacklist BLACKLIST [BLACKLIST ...]]
                                         [--blacklistFile BLACKLISTFILE]
                                         [--titleRegex TITLEREGEX]
                                         [--negativeTitleRegex NEGATIVETITLEREGEX]
                                         [--truncateTitlesAfter TRUNCATETITLESAFTER]
                                         [--minSequenceLen MINSEQUENCELEN]
                                         [--maxSequenceLen MAXSEQUENCELEN]
                                         [--taxonomy TAXONOMY]
                                         [--minMatchingReads MINMATCHINGREADS]
                                         [--maxMatchingReads MAXMATCHINGREADS]
                                         [--minMedianScore MINMEDIANSCORE]
                                         [--withScoreBetterThan WITHSCOREBETTERTHAN]
                                         [--minNewReads MINNEWREADS]
                                         [--maxTitles MAXTITLES]
                                         [--minCoverage MINCOVERAGE]
                                         [--sortOn {maxScore,medianScore,readCount,length,title}]
                                         [--rankValues RANKVALUES]
                                         --outputDir OUTPUTDIR [--color COLOR]
                                         [--equalizeXAxes]
                                         [--xRange {reads,subject}]
                                         [--logLinearXAxis]
                                         [--logBase LOGBASE] [--showFeatures]
                                         [--showOrfs] [--sortFilenames]
                                         [--titlesJSONFile TITLESJSONFILE]

Non-interactively generate an alignment panel

options:
  -h, --help            show this help message and exit
  --earlyExit           If True, just print the number of interesting matches,
                        but do not create the alignment panel. (default:
                        False)
  --matcher {blast,diamond}
                        The matching algorithm that was used to produce the
                        JSON. (default: blast)
  --json JSON-file [JSON-file ...]
                        the JSON file(s) of BLAST or DIAMOND output.
  --fasta FASTA-file [FASTA-file ...]
                        the FASTA file(s) of sequences that were given to
                        BLAST or DIAMOND. (default: None)
  --fastq FASTQ-file [FASTQ-file ...]
                        the FASTQ file(s) of sequences that were given to
                        BLAST or DIAMOND. (default: None)
  --databaseFastaFilename DATABASEFASTAFILENAME
                        Only used when --showOrfs is also given. The filename
                        of the FASTA file used to make the BLAST or DIAMOND
                        database. If --matcher diamond is used, either this
                        argument or --sqliteDatabaseFilename must be
                        specified. If --matcher blast is used these options
                        can be omitted, in which case the code will fall back
                        to using blastdbcmd, which can be unreliable. See also
                        --sqliteDatabaseFilename for a way to enable fast
                        subject lookup for either matcher. (default: None)
  --sqliteDatabaseFilename SQLITEDATABASEFILENAME
                        Only used when --showOrfs is also given. The filename
                        of the sqlite3 database file of FASTA metadata, made
                        from the FASTA that was used to make the BLAST or
                        DIAMOND database. If --matcher diamond is used, either
                        this argument or --databaseFilename must be specified.
                        (default: None)
  --databaseFastaDirectory DATABASEFASTADIRECTORY
                        Only used when --showOrfs is also given. The directory
                        where the FASTA file used to make the BLAST or DIAMOND
                        database can be found. This argument is only useful
                        when --sqliteDatabaseFilename is specified. (default:
                        None)
  --minStart MINSTART   Reads that start before this subject offset should not
                        be shown. (default: None)
  --maxStop MAXSTOP     Reads that end after this subject offset should not be
                        shown. (default: None)
  --oneAlignmentPerRead
                        If True, only keep the best alignment for each read.
                        (default: False)
  --maxAlignmentsPerRead MAXALIGNMENTSPERREAD
                        Reads with more than this many alignments will be
                        elided. Pass zero to only keep reads with no matches
                        (alignments). (default: None)
  --scoreCutoff SCORECUTOFF
                        A float score. Matches with scores worse than this
                        will be ignored. (default: None)
  --percentageIdenticalCutoff PERCENTAGEIDENTICALCUTOFF
                        A float percentage identity (0.0 to 100.0). Matches
                        with percent identity scores less than this will be
                        ignored. (default: None)
  --percentagePositiveCutoff PERCENTAGEPOSITIVECUTOFF
                        A float percentage (0.0 to 100.0). Matches with
                        percent positive scores less than this will be
                        ignored. A single AA match is considered positive (by
                        DIAMOND) if its entry in the substitution (scoring)
                        matrix is positive. DIAMOND uses the BLOSUM62 matrix
                        by default. (default: None)
  --maxHspsPerHit MAXHSPSPERHIT
                        A numeric max number of HSPs to show for each hit on
                        hitId. (default: None)
  --whitelist WHITELIST [WHITELIST ...]
                        sequence titles that should be whitelisted (default:
                        None)
  --whitelistFile WHITELISTFILE
                        A file containing sequence titles (one per line) to be
                        whitelisted (default: None)
  --blacklist BLACKLIST [BLACKLIST ...]
                        sequence titles that should be blacklisted (default:
                        None)
  --blacklistFile BLACKLISTFILE
                        A file containing sequence titles (one per line) to be
                        blacklisted (default: None)
  --titleRegex TITLEREGEX
                        a regex that sequence titles must match. (default:
                        None)
  --negativeTitleRegex NEGATIVETITLEREGEX
                        a regex that sequence titles must not match. (default:
                        None)
  --truncateTitlesAfter TRUNCATETITLESAFTER
                        a string that titles will be truncated beyond. If the
                        truncated version of a title has already been seen,
                        that title will be skipped. (default: None)
  --minSequenceLen MINSEQUENCELEN
                        sequences of lesser length will be elided. (default:
                        None)
  --maxSequenceLen MAXSEQUENCELEN
                        sequences of greater length will be elided. (default:
                        None)
  --taxonomy TAXONOMY   a string of the taxonomic group on which should be
                        filtered. eg "Vira" will filter on viruses. (default:
                        None)
  --minMatchingReads MINMATCHINGREADS
                        sequences that are matched by fewer reads will be
                        elided. (default: None)
  --maxMatchingReads MAXMATCHINGREADS
                        sequences that are matched by more reads will be
                        elided. (default: None)
  --minMedianScore MINMEDIANSCORE
                        sequences that are matched with a median score that is
                        worse will be elided. (default: None)
  --withScoreBetterThan WITHSCOREBETTERTHAN
                        sequences that are matched without at least one score
                        at least this good will be elided. (default: None)
  --minNewReads MINNEWREADS
                        The fraction of its reads by which a new read set must
                        differ from all previously seen read sets in order to
                        be considered acceptably different. (default: None)
  --maxTitles MAXTITLES
                        The maximum number of titles to keep. If more titles
                        than this result from the filtering, titles will be
                        sorted (according to the --sortOn value) and only the
                        best will be retained. (default: None)
  --minCoverage MINCOVERAGE
                        The (0.0 to 1.0) minimum fraction of a subject
                        sequence that must be matched by at least one read.
                        (default: None)
  --sortOn {maxScore,medianScore,readCount,length,title}
                        The attribute to sort subplots on. (default: maxScore)
  --rankValues RANKVALUES
                        If True, display reads with a Y axis coord that is the
                        rank of the score. (default: False)
  --outputDir OUTPUTDIR
                        Specifies a directory to write the HTML summary to.
                        Required.
  --color COLOR         a string which has a color as the first element and
                        read ids and/or FASTA file names as the following
                        element(s), separated by spaces. (default: None)
  --equalizeXAxes       If True, all alignment graphs will have their X axes
                        drawn with the same range. (default: False)
  --xRange {reads,subject}
                        Set the X axis range to show either the subject or the
                        extent of the reads that hit the subject. (default:
                        subject)
  --logLinearXAxis      If True, convert read offsets so that empty regions in
                        the alignment panel plots will only be as wide as
                        their logged actual values (default: False)
  --logBase LOGBASE     The base of the logarithm to use if logLinearXAxis is
                        True (default: 1.1)
  --showFeatures        If specified, look up features for the individual
                        images in the alignment panel. (default: False)
  --showOrfs            If specified, show subject ORFs in the individual
                        panel plots. This option requires that you also
                        provide information about the subject database, e.g.,
                        via --databaseFastaFilename. (default: False)
  --sortFilenames       If specified, the JSON and FASTA/Q file names will be
                        processed in sorted order. The sorting is based on
                        finding a numerical prefix in the filename. This can
                        be useful when processing output files produced by
                        systems like HTCondor, which makes files with names
                        like 1.out, 10.out, etc. that do not sort properly and
                        so cannot conveniently be given to this program along
                        with a single FASTA/Q file (because the order of the
                        results in the files from HTCondor does not match the
                        order of sequences in the FASTA/Q file. (default:
                        False)
  --titlesJSONFile TITLESJSONFILE
                        Give a file name for JSON holding information about
                        titles to be written to (default: None)

Given BLAST or DIAMOND JSON output files, the corresponding FASTA (or FASTQ)
sequence files, and filtering criteria, produce an alignment panel.
```

</details>

<details><summary><strong><code>proteins-to-pathogens.py</code></strong> — Group protein match output (as produced by the alignment-panel scripts
above) by the pathogen each protein comes from, using the pathogen
accession number embedded in the "civ"-style subject title. <span title="Show --help">ℹ️</span></summary>

```
usage: proteins-to-pathogens.py [-h] [--sampleName SAMPLENAME]
                                [--sampleNameRegex SAMPLENAMEREGEX]
                                [--pathogenPanelFilename PATHOGENPANELFILENAME]
                                [--sampleIndexFilename SAMPLEINDEXFILENAME]
                                [--pathogenIndexFilename PATHOGENINDEXFILENAME]
                                [--html] [--format {fasta,fastq}]
                                [--proteinFastaFilename PROTEINFASTAFILENAME [PROTEINFASTAFILENAME ...]]
                                [--minProteinFraction MINPROTEINFRACTION]
                                [--minProteinCount MINPROTEINCOUNT]
                                [--pathogenType {bacterial,viral}]
                                [--showReadLengths] [--assetDir ASSETDIR]
                                [--pathogenDataDir PATHOGENDATADIR]
                                [--title TITLE] [--preamble PREAMBLE]
                                [--titleRegex TITLEREGEX]
                                [--negativeTitleRegex NEGATIVETITLEREGEX]
                                [--omitVirusLinks] [--omitSampleProteinCount]
                                [filenames ...]

Group proteins by the pathogen they're from.

positional arguments:
  filenames             Sample file names to read input from. (default: None)

options:
  -h, --help            show this help message and exit
  --sampleName SAMPLENAME
                        An (optional) sample name. This is only used in
                        producing HTML output. Should be used when all input
                        files are for a single sample. Cannot be used with
                        --sampleNameRegex. (default: None)
  --sampleNameRegex SAMPLENAMEREGEX
                        An (optional) regular expression that can be used to
                        extract a short sample name from full sample file
                        name. The regular expression must have a matching
                        group (delimited by parentheses) that captures the
                        part of the file name that should be used as the
                        sample name. (default: None)
  --pathogenPanelFilename PATHOGENPANELFILENAME
                        An (optional) filename to write a pathogen-sample
                        panel PNG image to. (default: None)
  --sampleIndexFilename SAMPLEINDEXFILENAME
                        An (optional) filename to write a sample index file
                        to. Lines in the file will have an integer index, a
                        space, and then the sample name. Only produced if
                        --html is used (because the pathogen-NNN-sample-
                        MMM.fastq are only written in that case). (default:
                        None)
  --pathogenIndexFilename PATHOGENINDEXFILENAME
                        An (optional) filename to write a pathogen index file
                        to. Lines in the file will have an integer index, a
                        space, and then the pathogen name. Only produced if
                        --html is used (because the pathogen-NNN-sample-
                        MMM.fastq are only written in that case). (default:
                        None)
  --html                If specified, output HTML instead of plain text.
                        (default: False)
  --format {fasta,fastq}
                        Give the format of the sequence files written by
                        noninteractive-alignment-panel.py when it created the
                        summary-proteins files given on output. (default:
                        fasta)
  --proteinFastaFilename, --pff PROTEINFASTAFILENAME [PROTEINFASTAFILENAME ...]
                        Optional filename(s) giving the name of the FASTA
                        file(s) with the protein AA sequences with their
                        associated pathogens in square brackets. This is the
                        format used by NCBI for bacterial and viral reference
                        sequence protein files. If given, the contents of this
                        file will be used to determine how many proteins each
                        matched pathogen has. This makes it much easier to
                        spot significant matches (as opposed to those where,
                        say, just one protein from a pathogen is matched).
                        (default: None)
  --minProteinFraction MINPROTEINFRACTION
                        The minimum fraction of proteins in a pathogen that
                        must be matched by a particular sample in order for
                        that pathogen to be displayed for that sample.
                        (default: 0.0)
  --minProteinCount MINPROTEINCOUNT
                        The minimum number of proteins in a pathogen that must
                        be matched by a particular sample in order for that
                        pathogen to be displayed for that sample. (default: 0)
  --pathogenType {bacterial,viral}
                        Specify the pathogen type. This option only affects
                        the language used in HTML output. (default: viral)
  --showReadLengths     If specified, the HTML output (use --html to get this)
                        will contain the lengths of all reads that match
                        proteins for a pathogen. (default: False)
  --assetDir ASSETDIR   The output directory where noninteractive-alignment-
                        panel.py puts its HTML, plots and FASTA or FASTQ
                        files, needed for using --html. (default: out)
  --pathogenDataDir PATHOGENDATADIR
                        The directory where per-pathogen information (e.g.,
                        collected reads across all samples) should be written.
                        (default: pathogen-data)
  --title TITLE         The title to show at the top of the output. (default:
                        Summary of pathogens)
  --preamble PREAMBLE   Optional preamble text to show after the title.
                        (default: None)
  --titleRegex TITLEREGEX
                        A regex that pathogen names must match. (default:
                        None)
  --negativeTitleRegex NEGATIVETITLEREGEX
                        a regex that pathogen names must not match. (default:
                        None)
  --omitVirusLinks      If specified, the HTML output (use --html to get this)
                        for viruses will not contain links to ICTV and
                        ViralZone. This should be used when working with
                        viruses that do not yet have names that can be looked
                        up. (default: False)
  --omitSampleProteinCount
                        If specified, the HTML output (use --html to get this)
                        for viruses will not contain counts of the number of
                        proteins matched by each sample for a given pathogen.
                        This should be used when working with RVDB where there
                        are many sequences for some proteins and a sample
                        matches many of them, leading to incorrect reporting
                        of the number of proteins of a pathogen that are
                        matched by samples. (default: False)
```

</details>

<details><summary><strong><code>sam-to-fasta-alignment.py</code></strong> — Produce aligned FASTA queries from a SAM/BAM file. <span title="Show --help">ℹ️</span></summary>

```
usage: sam-to-fasta-alignment.py [-h] [--rcSuffix RCSUFFIX] [--rcNeeded]
                                 [--listReferenceInsertions] --samfile SAMFILE
                                 [--referenceId ID] [--dropUnmapped]
                                 [--dropSecondary] [--dropSupplementary]
                                 [--dropDuplicates] [--keepQCFailures]
                                 [--minScore FLOAT] [--maxScore FLOAT]
                                 [--scoreTag TAG] [--minLength N]
                                 [--maxLength N] [--maxNFraction N]
                                 [--whitelist SEQUENCE-ID]
                                 [--blacklist SEQUENCE-ID]
                                 [--whitelistFile SEQUENCE-ID-FILE]
                                 [--blacklistFile SEQUENCE-ID-FILE]
                                 [--titleRegex REGEX]
                                 [--negativeTitleRegex REGEX]
                                 [--sequenceWhitelist SEQUENCE]
                                 [--sequenceBlacklist SEQUENCE]
                                 [--sequenceWhitelistFile SEQUENCE-FILE]
                                 [--sequenceBlacklistFile SEQUENCE-FILE]
                                 [--sequenceRegex REGEX]
                                 [--sequenceNegativeRegex REGEX]
                                 [--keepSequences NUMBER,RANGE,... |
                                 --removeSequences NUMBER,RANGE,...]
                                 [--head N] [--removeDuplicates |
                                 --removeDuplicatesById]
                                 [--removeDuplicatesUseMD5] [--randomSubset N]
                                 [--trueLength N] [--sampleFraction FRACTION]
                                 [--sequenceNumbersFile FILENAME]

Produce aligned FASTA queries from a SAM/BAM file.

options:
  -h, --help            show this help message and exit
  --rcSuffix RCSUFFIX   A string to add to the end of query names that are
                        reverse complemented. This is added before the /1, /2,
                        etc., that are added for duplicated ids (if there are
                        duplicates and --allowDuplicateIds is not used)
                        (default: )
  --rcNeeded            If given, queries that are flagged as matching when
                        reverse complemented will be reverse complemented in
                        the output. This must be used if the program that
                        created the SAM/BAM input flags reversed matches but
                        does not also store the reverse complemented query.
                        The bwa program (mem and aln followed by samse) stores
                        the queries reversed complemented if the match was, so
                        this option is not needed for bwa. If in doubt, test
                        the output of your matching program as this is very
                        important! (default: False)
  --listReferenceInsertions
                        If given, information about reference sequence
                        insertions will be printed to standard error. These
                        correspond to "I" CIGAR operations that for the match
                        would require inserting query bases into the
                        reference. Because we cannot change the reference (in
                        fact we typically do not have the reference in the
                        SAM/BAM file), we cut the inserted bases out of the
                        aligned query and save the information about what
                        would have been inserted and where. That information
                        is printed by this option. The output gives the
                        0-based offset where the inserted base would be
                        placed, followed by a list of the nucleotides that
                        were suggested as being inserted and the number of
                        times each nucleotide was suggested. So for example
                        the output might contain "27: T:3, G:10" which
                        indicates that 13 query (3 with T and 10 with G)
                        matches would insert a nucleotide into the reference
                        at offset 27. (default: False)
  --samfile SAMFILE     The SAM/BAM file to filter.
  --referenceId ID      A reference sequence id whose alignments should be
                        kept (alignments against other references will be
                        dropped). If omitted, alignments against all
                        references will be kept. May be repeated. (default:
                        None)
  --dropUnmapped        If given, unmapped matches will not be output.
                        (default: False)
  --dropSecondary       If given, secondary matches will not be output.
                        (default: False)
  --dropSupplementary   If given, supplementary matches will not be output.
                        (default: False)
  --dropDuplicates      If given, matches flagged as optical or PCR duplicates
                        will not be output. (default: False)
  --keepQCFailures      If given, reads that are considered quality control
                        failures will be included in the output. (default:
                        False)
  --minScore FLOAT      If given, alignments with --scoreTag (default AS)
                        values less than this value will not be output. If
                        given, alignments that do not have a score will not be
                        output. (default: None)
  --maxScore FLOAT      If given, alignments with --scoreTag (default AS)
                        values greater than this value will not be output. If
                        given, alignments that do not have a score will not be
                        output. (default: None)
  --scoreTag TAG        The alignment tag to extract for --minScore and
                        --maxScore comparisons. (default: AS)
  --minLength N         The minimum sequence length. (default: None)
  --maxLength N         The maximum sequence length. (default: None)
  --maxNFraction N      The maximum fraction of Ns that can be present in the
                        sequence. (default: None)
  --whitelist SEQUENCE-ID
                        Sequence titles (ids) that should be whitelisted. May
                        be repeated. (default: None)
  --blacklist SEQUENCE-ID
                        Sequence titles (ids) that should be blacklisted. May
                        be repeated. (default: None)
  --whitelistFile SEQUENCE-ID-FILE
                        The name of a file that contains sequence titles (ids)
                        that should be whitelisted, one per line (default:
                        None)
  --blacklistFile SEQUENCE-ID-FILE
                        The name of a file that contains sequence titles (ids)
                        that should be blacklisted, one per line. (default:
                        None)
  --titleRegex REGEX    A regex that sequence titles (ids) must match.
                        (default: None)
  --negativeTitleRegex REGEX
                        A regex that sequence titles (ids) must not match.
                        (default: None)
  --sequenceWhitelist SEQUENCE
                        Sequences that should be whitelisted. May be repeated.
                        (default: None)
  --sequenceBlacklist SEQUENCE
                        Sequences that should be blacklisted. May be repeated.
                        (default: None)
  --sequenceWhitelistFile SEQUENCE-FILE
                        The name of a file that contains sequences that should
                        be whitelisted, one per line. (default: None)
  --sequenceBlacklistFile SEQUENCE-FILE
                        The name of a file that contains sequences that should
                        be blacklisted, one per line. (default: None)
  --sequenceRegex REGEX
                        A regex that sequences must match. (default: None)
  --sequenceNegativeRegex REGEX
                        A regex that sequences must not match. (default: None)
  --keepSequences NUMBER,RANGE,...
                        Specify (1-based) ranges of sequence numbers that
                        should be kept. E.g., --keepSequences 1-3,5 will
                        output just the 1st, 2nd, 3rd, and 5th sequences. All
                        others will be omitted. (default: None)
  --removeSequences NUMBER,RANGE,...
                        Specify (1-based) ranges of sequence numbers that
                        should be removed. E.g., --removeSequences 1-3,5 will
                        output all but the 1st, 2nd, 3rd, and 5th sequences.
                        All others will be ouput. (default: None)
  --head N              Only the first N sequences will be printed. (default:
                        None)
  --removeDuplicates    Duplicate reads will be removed, based only on
                        sequence identity. The first occurrence is kept.
                        (default: False)
  --removeDuplicatesById
                        Duplicate reads will be removed, based only on read
                        id. The first occurrence is kept. (default: False)
  --removeDuplicatesUseMD5
                        MD5 sums will be stored instead of the full sequence
                        or read id when either --removeDuplicates or
                        removeDuplicatesById are given. One of
                        --removeDuplicates and --removeDuplicatesById must
                        also be given. Note that this makes duplicate removal
                        probabilistic. This option can be used to reduce the
                        amount of RAM consumed during duplicate removal.
                        (default: False)
  --randomSubset N      An integer giving the number of sequences that should
                        be kept. These will be selected at random. (default:
                        None)
  --trueLength N        The number of reads in the FASTA input. Only to be
                        used with randomSubset (default: None)
  --sampleFraction FRACTION
                        A [0.0, 1.0] C{float} indicating a fraction of the
                        reads that should be allowed to pass through the
                        filter. The sample size will only be approximately the
                        product of the sample fraction and the number of
                        reads. The sample is taken at random. (default: None)
  --sequenceNumbersFile FILENAME
                        A file of (1-based) sequence numbers to retain.
                        Numbers must be one per line. (default: None)
```

</details>

# Incomplete

<details><summary><strong><code>simple-consensus.py</code></strong> — Make a consensus sequence directly from a BAM file using dark matter's
own (pure-Python) consensus code, rather than external tools. <span title="Show --help">ℹ️</span></summary>

```
usage: simple-consensus.py [-h] --bamFilename BAMFILENAME [--version]
                           [--bamId ID] [--referenceFasta FILE] [--fastaId ID]
                           [--consensusId ID | --idLambda "lambda id: ..."]
                           [--minCoverage N] [--lowCoverage N]
                           [--noCoverage NOCOVERAGE] [--threshold N]
                           [--ignoreQuality] [--includeSoftClipped]
                           [--compareWithPileupFile FILE]
                           [--deletionSymbol CHAR] [--strategy {fetch}]
                           [--deletionThreshold N]
                           [--insertionCountThreshold N] [--progress]
                           [--quiet]

Make a consensus sequence.

options:
  -h, --help            show this help message and exit
  --bamFilename BAMFILENAME
                        The BAM file from which the consensus should be
                        called.
  --version             Print the version number and exit (default: False)
  --bamId ID            The BAM file reference name indicating which aligned
                        reads to make a consensus from. If not given, will be
                        inferred from the BAM file header. (default: None)
  --referenceFasta FILE
                        The reference FASTA file. (default: None)
  --fastaId ID          The id of the sequence in --referenceFasta to use as a
                        reference. Only considered if --referenceFasta is
                        used. If not given and --referenceFasta is, the
                        reference id will be inferred from reference names in
                        the BAM header, or will be taken as the id of the
                        first sequence in --referenceFasta. (default: None)
  --consensusId ID      The id to use in the consensus sequence in the output
                        FASTA. If not given, the reference sequence id will be
                        used. (default: None)
  --idLambda "lambda id: ..."
                        A one-argument function taking and returning a read
                        id. This can be used to set the id of the consensus
                        sequence based on the id of the reference sequence.
                        The function will be called with the id of the BAM
                        reference sequence. E.g., --idLambda "lambda id:
                        id.split('_')[0]" or --idLambda "lambda id: id[:10] +
                        '-consensus'". (default: None)
  --minCoverage N       The minimum number of reads that must cover a site for
                        a consensus base to be called. If zero reads cover a
                        site, the --noCoverage value is used or if the number
                        is greater than zero but less than minCoverage, the
                        --lowCoverage value is used. (default: 0)
  --lowCoverage N       What to do when some reads cover a site, but fewer
                        than the --minCoverage value. Either set this to
                        "reference" (to use the base from the reference
                        sequence)or a single character (e.g., "N"). (default:
                        n)
  --noCoverage NOCOVERAGE
                        What to do when no reads cover a site, Either set this
                        to "reference" (to use the base from the reference
                        sequence), or a single character (e.g., "N").
                        (default: N)
  --threshold N         The frequency threshold when calling the consensus. If
                        the frequency of the most-common nucleotide at a site
                        meets this threshold, that nucleotide will be called.
                        Otherwise, an ambiguous nucleotide code will be
                        produced, based on the smallest set of most-frequent
                        nucleotides whose summed frequencies meet the
                        threshold. If the frequency of the nucleotide that
                        causes the threshold to be reached is the same as that
                        of other nucleotides, all such nucleotides will be
                        included in the ambiguous code. (default: 0.6)
  --ignoreQuality       Ignore FASTQ quality scores. All bases of all reads
                        will be considered equally reliable. (default: False)
  --includeSoftClipped  Include information from read bases that were marked
                        as soft-clipped by the algorithm that made the BAM
                        file. See https://samtools.github.io/hts-
                        specs/SAMv1.pdf if this makes no sense to you.
                        (default: False)
  --compareWithPileupFile FILE
                        Make a consensus using the pysam pileup function and
                        compare it to the one made by usign pysam fetch. Write
                        the summary to this file. If the file is "-", write to
                        standard error. (default: None)
  --deletionSymbol CHAR
                        When a deletion is detected, put this symbol into the
                        consensus sequence. Use the empty string to have
                        deletions omitted. (default: -)
  --strategy {fetch}    The consensus-making strategy to use. (default: fetch)
  --deletionThreshold N
                        If some reads have a deletion at a site and some do
                        not, call the site as a deletion if the fraction of
                        reads with the deletion is at least this value.
                        (default: 0.5)
  --insertionCountThreshold N
                        The number of reads that must have an insertion at an
                        offset in order for the insertion to be called in the
                        consensus. (default: 5)
  --progress            Show a progress bar (unless standard error has been
                        redirected). (default: False)
  --quiet               Suppress diagnostic output. Note that this will
                        silence warnings about differing reference names.
                        (default: False)
```

</details>

<details><summary><strong><code>convert-diamond-to-sam.py</code></strong> — Convert DIAMOND tabular output (produced with a specific required
<code>--outfmt 6</code> field order, printable via <code>--printFields</code>) to SAM format. <span title="Show --help">ℹ️</span></summary>

```
usage: convert-diamond-to-sam.py [-h] [--printFields]
                                 [--mappingQuality MAPPINGQUALITY] [--ram]
                                 [--keepDescriptions]
                                 --genomesProteinsDatabaseFilename GENOMESPROTEINSDATABASEFILENAME
                                 [--outputDirectory OUTPUTDIRECTORY]
                                 [--perReferenceOutput]
                                 [--filePointerCacheMaxSize FILEPOINTERCACHEMAXSIZE]

Convert DIAMOND tabular format to SAM. The DIAMOND invocation *must* include
--outfmt 6 bitscore btop qframe qend full_qqual qlen full_qseq qseqid qstart
slen sstart stitle

options:
  -h, --help            show this help message and exit
  --printFields         Print the field names in the order that they must be
                        given to diamond --outfmt 6 to produce correct input
                        for this script, then exit. (default: False)
  --mappingQuality MAPPINGQUALITY
                        The mapping quality to use for MAPQ (field 5). The
                        default (255) indicates that mapping quality
                        information is not available. (default: 255)
  --ram                 Do not use a temporary file to hold the non-header SAM
                        output. This will run faster but use more memory since
                        all non-header SAM output will be stored in RAM and
                        only written out when the full header can be
                        determined. (default: False)
  --keepDescriptions    Do not discard text after the first space in query or
                        subject sequence ids. Note that this violates the SAM
                        specification, but since SAM files are TAB-separated
                        it may not cause problems downstream, depending on
                        what tools you use on the SAM. (default: False)
  --genomesProteinsDatabaseFilename GENOMESPROTEINSDATABASEFILENAME
                        The name of the file containing an sqlite3 database
                        with information on genomes and proteins, as made by
                        make-protein-database.py
  --outputDirectory OUTPUTDIRECTORY
                        The name of the directory to create per-reference SAM
                        files in. This is only used when --perReferenceOutput
                        is given, otherwise the output is written to stdout as
                        a single SAM file. (default: .)
  --perReferenceOutput  If given, a separate SAM file is written for matches
                        against each reference sequence. The files will have
                        names corresponding to the accession number of the
                        reference, with a .sam suffix and will be put into the
                        directory specified by --outputDirectory. (default:
                        False)
  --filePointerCacheMaxSize FILEPOINTERCACHEMAXSIZE
                        The maximum size for the file pointer cache. SAM
                        output may be being written to an arbitrary number of
                        output files so we need to make sure the OS limit on
                        open files is not exceeded. The default value is size
                        half the maximum number of open files will be used.
                        (default: 524288)
```

</details>
