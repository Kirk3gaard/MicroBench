README
================
Rasmus Kirkegaard
02 October, 2024

<img align="right" src="images/microbench.jpg" width="200">

# MicroBench

Overview of our data for microbial genomic benchmarking. We have
sequenced a bunch of mono cultures, and metagenomic samples.

Monocultures:

- [Anabaena variabilis PCC 7120 DSM
  107007](#anabaena-variabilis-pcc-7120-dsm-107007)
- [Bacillus cereus NRS 248 ATCC
  10987](bacillus-cereus-nrs-248-atcc-10987)
- [Burkholderia cenocepacia J2315 DSM 16553]()
- [Cellulophaga lytica Cy l20, LIM -21 DSMZ 7489]()
- [Desulfobacca acetoxidans ASRB2 DSM 11109]()
- [Helicobacter pylori J99 ATCC 700824]()
- [Helicobacter pylori 26695 ATCC 700392]()
- [Methanocaldococcus jannaschii JAL-1 DSM 2661]()
- [Neisseria gonorrhoeae FA1090 ATCC 700825]()
- [Pelobacter carbinolicus Bd1 DSM 2380]()
- [Sagittula stellata E-37 DSM 11524]()
- [Salmonella bongori 1224.72 DSM 13772]()
- [Sphaerobacter thermophilus S 6022 DSM 20745]()
- [Sphaerotilus natans 6 DSM 6575]()
- [Thermanaerovibrio acidaminovorans Su883 DSM 6589]()
- [Zymomonas mobilis subsp. pomaceae Barker I DSM 22645]()

Mock Metagenomes:

- [ZymoHMW](#zymohmw)
- [ZymoGut](#zymogut)
- [ZymoOral](#zymooral)

“Real” metagenomes:

- [ZymoFecal](#zymofecal)
- [AAU-Anaerobic Digester standard]()
- [AAU-WTTP standard]()
- [AAU-Soil standard]()

## To do:

- Sequence all standards
- Get AAU to support public links or upload data to NCBI/ENA

## Mono cultures

### Anabaena variabilis PCC 7120 DSM 107007

You find a description of this sample here:
<https://www.dsmz.de/collection/catalogue/details/culture/DSM-107007>

| Sample | date of sequencing | sequencing kit | flowcell type | flowcell ID | sample rate | date of basecalling | dorado version | basecall model version | fast | hac | hacdup | sup | supdup | mods |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| Anabaena variabilis PCC 7120 DSM 107007 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |

### Bacillus cereus NRS 248 ATCC 10987

You find a description of this sample here:
<https://www.atcc.org/products/10987d-5>

| Sample | date of sequencing | sequencing kit | flowcell type | flowcell ID | sample rate | date of basecalling | dorado version | basecall model version | fast | hac | hacdup | sup | supdup | mods |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| Bacillus cereus NRS 248 ATCC 10987 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |

## Mock metagenomes

### ZymoHMW

You find a description of this sample here:
<https://zymoresearch.eu/products/zymobiomics-hmw-dna-standard?srsltid=AfmBOooNJLjOVzGjtNm0BJY2WoMS0vHsUDi7mR59p4HDqyBX7FlJe7jy>

Some brief analysis of this data can be found here:
[analysis/zymohmw](https://github.com/Kirk3gaard/MicroBench/tree/main/analysis/zymohmw)

| Sample | date of sequencing | sequencing kit | flowcell type | flowcell ID | sample rate | date of basecalling | dorado version | basecall model version | fast | hac | hacdup | sup | supdup | mods |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| ZymoHMW | 2023-08-08 | LSK114 | FLO-PRO114HD | PAS01578 | 5 khz | 2024-09-13 | 0.7.3 | 5.0.0 | /<projects/MicroBench/data/PAS01578.dorado0.7.3.bmdna_r10.4.1_e8.2_400bps@5.0.0>\_fast.sim.fastq.gz | /<projects/MicroBench/data/PAS01578.dorado0.7.3.bmdna_r10.4.1_e8.2_400bps@5.0.0>\_hac.sim.fastq.gz | /<projects/MicroBench/data/PAS01578.dorado0.7.3.bmdna_r10.4.1_e8.2_400bps@5.0.0>\_hac.dup.fastq.gz | /<projects/MicroBench/data/PAS01578.dorado0.7.3.bmdna_r10.4.1_e8.2_400bps@5.0.0>\_sup.sim.fastq.gz | /<projects/MicroBench/data/PAS01578.dorado0.7.3.bmdna_r10.4.1_e8.2_400bps@5.0.0>\_sup.dup.fastq.gz | /<projects/MicroBench/data/PAS01578.dorado0.7.3.bmdna_r10.4.1_e8.2_400bps@5.0.0>\_sup.sim.mod4mC_5mC_6mA.bam |
| ZymoHMW | 2023-08-08 | LSK114 | FLO-PRO114HD | PAS01578 | 5 khz | 2023-12-16 | 0.3.4 | 4.2.0 | NA | NA | NA | <ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR119/ERR11901474/PAS01578.dna_r10.4.1_e8.2_400bps_sup@v4.2.0.fastq.gz> | <ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR119/ERR11901474/PAS01578.dna_r10.4.1_e8.2_400bps_sup@v4.2.0.fastq.gz> | NA |

### ZymoGut

You find a description of this sample here:
<https://zymoresearch.eu/collections/zymobiomics-microbial-community-standards/products/zymobiomics-gut-microbiome-standard>

| Sample | date of sequencing | sequencing kit | flowcell type | flowcell ID | sample rate | date of basecalling | dorado version | basecall model version | fast | hac | hacdup | sup | supdup | mods |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| ZymoGut | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |

### ZymoOral

You find a description of this sample here:
<https://zymoresearch.eu/collections/zymobiomics-microbial-community-standards/products/zymobiomics-oral-microbiome-standard>

| Sample | date of sequencing | sequencing kit | flowcell type | flowcell ID | sample rate | date of basecalling | dorado version | basecall model version | fast | hac | hacdup | sup | supdup | mods |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| ZymoOral | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |

## “Real” metagenomes

### ZymoFecal

You find a description of this sample here:
<https://zymoresearch.eu/products/zymobiomics-fecal-reference-with-trumatrix-technology?srsltid=AfmBOorYH_EvSWZvzNTVTAb7NYQk_zlNcRO4NXzI8WCiUPgG13fQBJhd>

| Sample | date of sequencing | sequencing kit | flowcell type | flowcell ID | sample rate | date of basecalling | dorado version | basecall model version | fast | hac | hacdup | sup | supdup | mods |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| ZymoFecal | 2024-10-01 | LSK114 | FLO-PRO114M | PAW77640 | 5 khz | NA | NA | NA | NA | NA | NA | NA | NA | NA |
