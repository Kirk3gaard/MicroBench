# Read artefacts in basecalled nanopore data

Code behind <https://kirk3gaard.github.io/MicroBench/artefacts.html>. It measures artefacts that are still
present in dorado-basecalled reads (after dorado's default adapter trimming and read splitting):
chimeras between and within genomes, foldback reads, internal and end adapters, low-complexity
junk, skip-prone reads, and how error rates change through a run.

Requirements: Python 3 with `pysam`, `mappy` and `edlib` (`pip install pysam mappy edlib`).

```bash
# 1. sample complete BAM records from the public BAM (HTTP range requests, ~300 MB download)
python3 sample_records.py \
  https://ftp.sra.ebi.ac.uk/vol1/run/ERR154/ERR15459724/PAS01578.dorado1.0.0.bm5.2.0_sup.sim.bam \
  44000000000 300 hmw_sup520.bam

# 2. reference: ZymoBIOMICS D6322 genomes, concatenated, plus a contig -> species table
curl -O https://s3.amazonaws.com/zymo-files/BioPool/D6322.refseq.zip && unzip D6322.refseq.zip
cat D6322.refseq/Genomes/*.fasta > zymo_hmw.fa
for f in D6322.refseq/Genomes/*.fasta; do sp=$(basename $f | sed 's/_complete_genome.fasta//;s/_draft_genome.fasta//')
  grep ">" $f | sed 's/^>//;s/ .*//' | awk -v s=$sp '{print $1"\t"s}'; done > zymo_hmw.fa.species

# 3. align and classify every read
python3 analyse.py hmw_sup520.bam zymo_hmw.fa hmw          # -> hmw.reads.jsonl

# 4. synchrony test for deletions across channels (min length, window in seconds)
for m in 10 20 50; do python3 sync.py hmw.reads.jsonl $m 2; done > sync.jsonl

# 5. page data
python3 build_page_data.py hmw.reads.jsonl hmw_sup520.bam ../../docs/data/artefacts_zymohmw.json sync.jsonl
```
