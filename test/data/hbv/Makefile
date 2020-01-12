REFERENCE := KJ586809.1
FASTA := $(REFERENCE).fasta

FASTQ1 := query1.fastq
FASTQ2 := query2.fastq
FASTQ3 := query3.fastq

SAM1 := query1.sam
SAM2 := query2.sam
SAM3 := query3.sam

BAM1 := query1.bam
BAM2 := query2.bam
BAM3 := query3.bam

SITES1 := 1-200
SITES2 := 1401-1500
SITES3 := 1001-1450

LEN1 := 200
LEN2 := 100
LEN3 := 450


all: hbv-proteins.db $(REFERENCE).1.bt2 $(BAM1) $(BAM2) $(BAM3)

hbv-proteins.db: ids_aaaaa.gb
	env PYTHONPATH=../../.. make-protein-database.py --gb ids_aaaaa.gb --databaseFile $@ --databaseName NCBI --logFile hbv-proteins.log > hbv-proteins.fasta

$(REFERENCE).1.bt2: $(FASTA)
	bowtie2-build $(FASTA) $(REFERENCE) >/dev/null

ids_aaaaa.gb:
	echo $(REFERENCE) | download-genbank.sh

$(FASTA):
	ncbi-fetch-id.py $(REFERENCE) > $@

$(FASTQ1): $(FASTA)
	echo '@query1' > $@
	env PYTHONPATH=../../.. filter-fasta.py --quiet --keepSites $(SITES1) < $(FASTA) | tail -n 1 >> $@
	echo '+' >> $@
	python -c "print('x' * $(LEN1))" >> $@

$(FASTQ2): $(FASTA)
	echo '@query2' > $@
	env PYTHONPATH=../../.. filter-fasta.py --quiet --keepSites $(SITES2) < $(FASTA) | tail -n 1 >> $@
	echo '+' >> $@
	python -c "print('x' * $(LEN2))" >> $@

$(FASTQ3): $(FASTA)
	echo '@query3' > $@
	env PYTHONPATH=../../.. filter-fasta.py --quiet --keepSites $(SITES3) < $(FASTA) | tail -n 1 >> $@
	echo '+' >> $@
	python -c "print('x' * $(LEN3))" >> $@

$(SAM1): $(FASTQ1) $(REFERENCE).1.bt2
	bowtie2 --quiet -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 --end-to-end -x $(REFERENCE) $(FASTQ1) > $@

$(SAM2): $(FASTQ2) $(REFERENCE).1.bt2
	bowtie2 --quiet -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 --end-to-end -x $(REFERENCE) $(FASTQ2) > $@

$(SAM3): $(FASTQ3) $(REFERENCE).1.bt2
	bowtie2 --quiet -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 --end-to-end -x $(REFERENCE) $(FASTQ3) > $@

$(BAM1): $(SAM1)

$(BAM2): $(SAM2)

$(BAM3): $(SAM3)

%.bam: %.sam
	samtools view -O BAM $< | samtools sort 2>/dev/null > $@
	samtools index $@

clean:
	rm -f ids_aaaaa.gb $(REFERENCE).* \
            hbv-proteins.fasta hbv-proteins.log \
            query[123].fastq query[123].sam *~

clobber: clean
	rm -f hbv-proteins.db query[123].bam query[123].bam.bai
