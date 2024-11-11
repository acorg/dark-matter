import multiprocessing
from tempfile import mkdtemp
import os
from os.path import join
import requests
from typing import Optional, TextIO

from dark.process import Executor
from dark.sam import samfile as openSamfile


class Bowtie2:
    """
    Run Bowtie2.
    """

    def __init__(
        self,
        executor: Optional[Executor] = None,
        threads: Optional[int] = None,
        verboseFp: Optional[TextIO] = None,
        dryRun: bool = False,
        reference: Optional[str] = None,
        tempdir: Optional[str] = None,
        tmpChmod: Optional[str] = None,
    ) -> None:
        self._executor = executor or Executor(dryRun)
        if dryRun:
            self.tempdir = tempdir or "/tmp/xxx"
        else:
            self.tempdir = mkdtemp(prefix="bt2-", dir=tempdir)
            if tmpChmod:
                self._executor.execute(f"chmod {tmpChmod} {self.tempdir}")
        self._samFile = join(self.tempdir, "result.sam")
        self._bamFile = join(self.tempdir, "result.bam")
        self._indexFile = join(self.tempdir, "index")
        self._verboseFp = verboseFp
        self._indexCalled = self._samExists = self._bamExists = False
        self.nThreads = threads or multiprocessing.cpu_count()
        self._reference = reference

    def buildIndex(self, index: str) -> None:
        """
        Find or make a Bowtie2 index.
        """
        if os.path.exists(index):
            # Check if this is a pre-existing bowtie2 index. Look for and
            # remove a bowtie2 suffix, if any (otherwise bowtie2 will
            # complain). We do things this way to allow a user to use TAB
            # completion on the command line to give them the full path to
            # any bowtie index file.
            for suffix in ("1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2").split():
                suffix = "." + suffix
                if index.endswith(suffix):
                    self._indexFile = index[: -len(suffix)]
                    self._report(
                        "Using pre-existing Bowtie2 index %r." % self._indexFile
                    )
                    break
            else:
                # Assume a FASTA file and make an index.
                self._indexFile = self._makeIndexFromFastaFile(index)
        else:
            # Not a filename. So either the start of the path to a bowtie2
            # index or else an accession number.
            if os.path.exists(index + ".1.bt2"):
                self._report("Using pre-existing Bowtie2 index %r." % index)
                self._indexFile = index
            else:
                # Assume an accession number.
                self._indexFile = self._makeIndexFromAccession(index)

        self._indexCalled = True

    def align(
        self,
        bowtie2Args: str = "--no-unal",
        fastq1: Optional[str] = None,
        fastq2: Optional[str] = None,
        threads: Optional[int] = None,
        discardSAM: bool = False,
        readGroup: str = "orig",
        sampleName: str = "orig",
    ) -> None:
        """
        Run Bowtie2.
        """
        if not self._indexCalled:
            raise ValueError("buildIndex() has not yet been called.")

        self._report("Aligning with Bowtie2.")

        nThreads = threads or self.nThreads
        samFile = "/dev/null" if discardSAM else self._samFile

        if fastq1 and fastq2:
            self._executor.execute(
                "bowtie2 %s --threads %d --rg-id '%s' --rg 'SM:%s' -x '%s' "
                "-1 '%s' -2 '%s' > '%s'"
                % (
                    bowtie2Args,
                    nThreads,
                    readGroup,
                    sampleName,
                    self._indexFile,
                    fastq1,
                    fastq2,
                    samFile,
                )
            )
        elif fastq1:
            self._executor.execute(
                "bowtie2 %s --threads %d --rg-id '%s' --rg 'SM:%s' -x '%s' "
                "-U '%s' > '%s'"
                % (
                    bowtie2Args,
                    nThreads,
                    readGroup,
                    sampleName,
                    self._indexFile,
                    fastq1,
                    samFile,
                )
            )
        else:
            raise ValueError("At least fastq1 must be passed.")

        self._samExists = True

    def _report(self, mesg: str) -> None:
        """
        Print a progress message.
        """
        if self._verboseFp:
            print(mesg, file=self._verboseFp)

    def makeBAM(self, samtoolsViewArgs: str = "") -> None:
        """
        Convert the SAM to BAM.
        """
        # The following will raise if there is no SAM file.
        self._SAMorBAM()
        self._report("Converting SAM to BAM.")
        self._executor.execute(
            "samtools view -b %s '%s' > '%s'"
            % (samtoolsViewArgs, self._samFile, self._bamFile)
        )
        self._bamExists = True

    def _SAMorBAM(self) -> str:
        """
        Do we have a SAM or BAM file available?
        """
        if self._bamExists:
            return "BAM"
        elif self._samExists:
            return "SAM"
        else:
            raise ValueError("bowtie2() has not yet been called.")

    def outputFile(self) -> str:
        """
        The name of the output file.
        """
        if self._bamExists:
            return self._bamFile
        elif self._samExists:
            return self._samFile
        else:
            raise ValueError("bowtie2() has not yet been called.")

    def indexBAM(self) -> None:
        """
        Index the BAM file.
        """
        which = self._SAMorBAM()

        if which != "BAM":
            raise ValueError("makeBAM() has not yet been called.")

        self._report("Indexing BAM.")
        self._executor.execute("samtools index '%s'" % self._bamFile)

    def isEmpty(self) -> bool:
        """
        Are any reads present in the SAM file?
        """
        # Reading the SAM file manually doesn't require BAM that is sorted and indexed
        # (which is what we would need to do if we wanted to use the mapped and unmapped
        # attributes of an AlignmentFile instance).
        with open(self._samFile) as fp:
            for line in fp:
                if not line.startswith("@"):
                    return False
        return True

    def mappedAndUnmappedCounts(self) -> tuple[int, int]:
        """
        Return the number of mapped and unmapped reads.

        @return: A 2-C{tuple} of C{int}s, giving the number of mapped and unmapped
            reads.
        """
        which = self._SAMorBAM()

        if which != "BAM":
            raise ValueError("makeBAM() has not yet been called.")

        with openSamfile(self._bamFile) as samfile:
            # This will fail if the BAM file has not been sorted and indexed.
            return samfile.mapped, samfile.unmapped

    def sort(self, byName=False) -> None:
        """
        Sort the BAM or SAM.
        """
        which = self._SAMorBAM()
        self._report("Sorting %s (by %s)." % (which, "name" if byName else "coord"))
        inFile = self._bamFile if which == "BAM" else self._samFile
        sortedFile = join(self.tempdir, "result-sorted." + which.lower())
        self._executor.execute(
            "samtools sort --output-fmt %s %s'%s' > '%s'"
            % (which, "-n " if byName else "", inFile, sortedFile)
        )
        self._executor.execute("mv '%s' '%s'" % (sortedFile, inFile))

    def removePrimers(self, bedFile: str) -> None:
        """
        Removes primers specified in the bed file
        """
        which = self._SAMorBAM()

        if which != "BAM":
            raise ValueError("makeBAM() has not yet been called.")

        self._report("removing primers specified in %s" % bedFile)
        tempTrimmedBamPrefix = "%s.trimmed" % self._bamFile
        self._executor.execute(
            "ivar trim -b '%s' -p '%s' -i '%s' -q 20 -m 30 -s 4 -e"
            % (bedFile, tempTrimmedBamPrefix, self._bamFile)
        )
        self._executor.execute(
            "mv '%s'.bam '%s'" % (tempTrimmedBamPrefix, self._bamFile)
        )

    def markDuplicatesPicard(self, picardFile: str) -> None:
        """
        Use Picard to mark duplicates.
        """
        which = self._SAMorBAM()
        self._report("Marking duplicates with Picard.")

        inFile = self._bamFile if which == "BAM" else self._samFile
        tempFile = join(self.tempdir, "picard-duplicates." + which.lower())
        tempErrFile = join(self.tempdir, "picard.errs")

        self._executor.execute(
            "java -Xmn2g -Xms2g -Xmx2g -jar %s "
            "MarkDuplicates I='%s' O='%s' M=/dev/null >'%s' 2>&1"
            % (picardFile, inFile, tempFile, tempErrFile)
        )

        self._executor.execute("mv '%s' '%s'" % (tempFile, inFile))

    def markDuplicatesGATK(self, threads: Optional[int] = None) -> None:
        """
        Use GATK to mark duplicates.
        """
        nThreads = threads or self.nThreads
        which = self._SAMorBAM()
        self._report("Marking duplicates with GATK.")

        inFile = self._bamFile if which == "BAM" else self._samFile
        tempFile = join(self.tempdir, "gatk-duplicates." + which.lower())

        self._executor.execute(
            "gatk MarkDuplicatesSpark "
            "-I '%s' -O '%s' --conf spark.executor.cores=%d"
            % (inFile, tempFile, nThreads)
        )

        self._executor.execute("mv '%s' '%s'" % (tempFile, inFile))

    def callHaplotypesGATK(
        self,
        picardJar: str,
        vcfFile: Optional[str] = None,
        referenceFasta: Optional[str] = None,
    ) -> None:
        """
        Use GATK to call haplotypes.
        """
        which = self._SAMorBAM()
        self._report("Calling haplotypes with GATK.")

        inFile = self._bamFile if which == "BAM" else self._samFile
        vcfFile = vcfFile or join(self.tempdir, "output.vcf.gz")

        if referenceFasta is None:
            if self._reference:
                self._report("Using %s as a reference." % self._reference)
                referenceFasta = self._reference
            else:
                raise ValueError(
                    "No reference was passed, given to the "
                    "Bowtie2 __init__, or used in buildIndex."
                )

        indexFile = referenceFasta + ".fai"
        if os.path.exists(indexFile):
            removeIndex = False
        else:
            removeIndex = True
            self._executor.execute("samtools faidx '%s'" % referenceFasta)

        if referenceFasta.lower().endswith(".fasta"):
            dictFile = referenceFasta[: -len(".fasta")] + ".dict"
        else:
            dictFile = referenceFasta + ".dict"

        if os.path.exists(dictFile):
            removeDict = False
        else:
            removeDict = True
            self._executor.execute(
                "java -jar '%s' CreateSequenceDictionary R='%s' O='%s'"
                % (picardJar, referenceFasta, dictFile)
            )

        self._executor.execute(
            "gatk --java-options -Xmx4g HaplotypeCaller "
            "--reference '%s' "
            "--input '%s' "
            "--output '%s' "
            "--sample-ploidy 1 "
            "--dont-use-soft-clipped-bases true "
            "-ERC GVCF" % (referenceFasta, inFile, vcfFile)
        )

        if removeIndex:
            self._executor.execute("rm '%s'" % indexFile)

        if removeDict:
            self._executor.execute("rm '%s'" % dictFile)

    def callHaplotypesBcftools(
        self, vcfFile: Optional[str] = None, referenceFasta: Optional[str] = None
    ) -> None:
        """
        Use bcftools call to call haplotypes.
        """
        which = self._SAMorBAM()
        self._report("Calling haplotypes with bcftools call.")

        inFile = self._bamFile if which == "BAM" else self._samFile
        vcfFile = vcfFile or join(self.tempdir, "output.vcf.gz")

        if referenceFasta is None:
            if self._reference:
                self._report("Using %s as a reference." % self._reference)
                referenceFasta = self._reference
            else:
                raise ValueError(
                    "No reference was passed, given to the "
                    "Bowtie2 __init__, or used in buildIndex."
                )

        self._executor.execute(
            'bcftools mpileup --max-depth 5000 -Ou -f "%s" "%s" | '
            'bcftools call --ploidy 1 -mv -Oz -o "%s"'
            % (referenceFasta, inFile, vcfFile)
        )

        self._executor.execute(f"bcftools index --force {vcfFile!r}")

    def removeDuplicates(self) -> None:
        """
        Use samtools to remove marked duplicates.
        """
        which = self._SAMorBAM()
        self._report("Removing marked duplicates.")

        inFile = self._bamFile if which == "BAM" else self._samFile
        tempFile = join(self.tempdir, "non-duplicates." + which.lower())

        # See the comment on the testRemoveDuplicates test in
        # ../test/test_bowtie2.py regarding there being two spaces in the
        # samtools command below if we are not producing BAM.  If you
        # change spacing here you will also need to change that test.
        bArg = "-b" if which == "BAM" else ""
        self._executor.execute(
            "samtools view %s -F 1024 '%s' > '%s'" % (bArg, inFile, tempFile)
        )
        self._executor.execute("mv '%s' '%s'" % (tempFile, inFile))

    def _makeIndexFromFastaFile(self, fastaFilename: str) -> str:
        index = join(self.tempdir, "index")

        self._report("Building Bowtie2 index from %s." % fastaFilename)
        self._executor.execute(
            "bowtie2-build --quiet '%s' '%s'" % (fastaFilename, index)
        )

        if self._reference is None:
            self._report("Will use this FASTA file as the reference (if needed).")
            self._reference = fastaFilename

        return index

    def _makeIndexFromAccession(self, accessionId: str) -> str:
        fastaFilename = join(self.tempdir, accessionId + ".fasta")
        index = join(self.tempdir, "index")

        self._report("Downloading FASTA for accession %s from NCBI." % accessionId)

        URL = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
            "db=nucleotide&id=%s&rettype=fasta&retmode=text" % accessionId
        )

        if not self._executor.dryRun:
            with open(fastaFilename, "w") as fp:
                print(requests.get(URL).text.rstrip("\n"), file=fp)

        self._report("Building bowtie2 index for %s." % accessionId)
        self._executor.execute(
            "bowtie2-build --quiet '%s' '%s'" % (fastaFilename, index)
        )

        return index

    def close(self) -> None:
        """
        Clean up.
        """
        self._executor.execute("rm -r '%s'" % self.tempdir)
