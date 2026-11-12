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

## Whole-run scan (all 12.6 M reads)

The page now uses every read of the run. On a 12-core Mac this takes ~45 min after downloading the 44 GB BAM
(`aria2c -x 16`, MD5 checked against ENA). Extra requirements: `pandas numpy scipy`.

```bash
python3 fullrun.py sup520.bam zymo_hmw.fa full 11            # per-read table, indel events, artefact layouts
python3 fullstats.py full full_stats.json zymo_hmw.fa        # summaries, per-second synchrony test, channels
python3 build_full_page.py full sup520.bam full_stats.json artefacts_zymohmw_full.json   # examples
python3 scan.py          # detects the periodic MinKNOW scans and bins reads by time since the last scan -> scan_stats.json
python3 saw.py           # 5-min error/junk series -> saw.json
python3 junk_time.py     # timing and position of low-complexity stretches -> junk_time.json
python3 scan_examples.py # example reads that start with junk right after a scan -> scan_examples.json
python3 occ.py           # artefact rates vs active channels and per channel (printed)
```
The four JSON files are merged into `docs/data/artefacts_zymohmw.json` (keys `saw`, `scan`, `examples.post_scan_junk`).

## Signal gaps, molecules read twice, second run

`fullrun.py` also writes collinear split alignments and chimera junctions as timed events.

```bash
python3 gaps.py full hmw_gaps.json [ontrim_report.html]   # flow-cell-wide signal gaps (+ comparison with ontrim's list)
python3 duplex.py full hmw_duplex.json                    # complement strand read right after the template in the same pore
```

ZymoGut (ERR15459768, P2 Solo) was run with the same scripts using the D6331 references
(`D6331.refseq.zip` from Zymo; E. coli strains collapsed to one species) and a `gut_` prefix for the scan scripts.
`docs/data/artefacts_index.json` lists the runs shown on the page.

## Kit inference

`scripts/kit_detect.py` samples reads from one BAM per sample and looks for native barcode constructs,
rapid barcode constructs, the rapid adapter primer and ligation adapter remnants (sequences from dorado's kit
definitions, `scripts/ont_barcodes.json`). Results are cached in `scripts/.cache/kits.json` and shown on the
datasets page as "Kit from reads", with a warning where barcodes are still in the reads.

## Headline numbers

`count_numbers.py <prefix> out.json [page.json]` recomputes every total shown at the top of the page from the
per-read table, counting each read once (unions, not sums of overlapping classes), and checks the page JSON's
per-class numbers against it. The result is stored in the page JSON under `numbers`.

## More runs, recurrent breaks, streaming

- `post.sh DIR PREFIX SEQS REF TAG DATASET_JSON [ONTRIM_HTML]` runs every step after `fullrun.py`
  (refine_chimeras, gaps, duplex, fullstats, scan, saw, junk_time, scan_examples, build_full_page, count_numbers);
  `merge_page.py` then assembles `docs/data/artefacts_<run>.json` and updates `artefacts_index.json`.
- `refine_chimeras.py` moves reads whose chimera junctions all sit at recurrent genome positions
  (strain differences, repeats, ends of linear reference contigs) out of the chimera classes into
  `recurrent_junction`, so that only junctions at random positions count as chimeras.
- `rangecat.py URL SIZE | fullrun.py - ...` streams a BAM too large to store through the scan
  (parallel HTTP range requests via curl, MD5 checked on the fly); `POOL_FRAC=0.02` keeps the full
  sequence of 2% of artefact reads for example plots.
- ZymoOral (D6332 references) runs exactly like ZymoGut. ZymoFecal has no reference genomes: its
  reference is assembled from its own reads with metaMDBG (complete circular contigs ≥1 Mb named,
  the rest pooled as "unbinned", and joins involving unbinned contigs are not called chimeras).
