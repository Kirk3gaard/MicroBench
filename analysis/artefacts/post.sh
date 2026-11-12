#!/bin/zsh
# post-processing of a whole-run scan: post.sh DIR PREFIX SEQS REF TAG DATASET_JSON [ONTRIM_HTML]
# SEQS = the BAM (or the .pool.fastq.gz when the BAM was streamed)
set -u
B=~/microbench-push/bin; PY=~/micromamba/envs/mbart/bin/python3
cd "$1"; P=$2; S=$3; REF=$4; T=$5; DS=$6; ONT=${7:-}
log() { echo "$(date '+%F %T') $*" >> post.log; }
log start
$PY $B/refine_chimeras.py $P $REF > refine.log 2>&1 || log "refine failed"
$PY $B/gaps.py $P ${T}gaps.json $ONT > gaps.log 2>&1 || log "gaps failed"
$PY $B/duplex.py $P ${T}duplex.json > dup.log 2>&1 || log "duplex failed"
$PY $B/fullstats.py $P ${T}stats.json $REF > stats.log 2>&1 || log "stats failed"
$PY $B/scan.py $P $S $T > scan.log 2>&1 || log "scan failed"
$PY $B/saw.py $P $S $T > saw.log 2>&1 || log "saw failed"
$PY $B/junk_time.py $P $S $T > junk.log 2>&1 || log "junk failed"
$PY $B/scan_examples.py $P $S $T > scanex.log 2>&1 || log "scan examples failed"
$PY $B/build_full_page.py $P $S ${T}stats.json ${T}page.json "$DS" > page.log 2>&1 || log "page failed"
$PY $B/count_numbers.py $P ${T}numbers.json > numbers.log 2>&1 || log "numbers failed"
log done
