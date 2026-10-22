#!/bin/bash
# Restore a case folder to the clean starting state (checkpoint at step 1680000 + SEM eddies)
# and set the last time step, so that runs can be compared from the same state.
#
# usage: ./restore_case.sh <case_dir> <nstep>      e.g.  ./restore_case.sh prof300 1680300
#
# Every run overwrites data/restart_dir and sem_restart_1.txt at the end, hence this script.

BACKUP=/leonardo_scratch/large/userexternal/aroccon0/checkpoints/single_nuc_site_1680000

CASE=$1
NSTEP=$2
if [ -z "$CASE" ] || [ -z "$NSTEP" ]; then
  echo "usage: $0 <case_dir> <nstep>"; exit 1
fi
if [ ! -d "$BACKUP" ]; then
  echo "backup not found: $BACKUP"; exit 1
fi
if [ ! -f "$CASE/dns.in" ] || [ ! -f "$CASE/go_leo.sh" ]; then
  echo "$CASE is not a complete case folder (dns.in or go_leo.sh missing): nothing changed"; exit 1
fi

rm -rf "$CASE/data/restart_dir"
mkdir -p "$CASE/data"
cp -r "$BACKUP" "$CASE/data/restart_dir"
rm -f "$CASE/data/restart_dir/sem_restart_1.txt"
cp "$BACKUP/sem_restart_1.txt" "$CASE/"
sed -i "11s/^[0-9]*/$NSTEP/" "$CASE/dns.in"

echo "checkpoint : $(cat $CASE/data/restart_dir/scalar.out)"
echo "dns.in l.11: $(sed -n 11p $CASE/dns.in)"
echo "dns.in l.32: $(sed -n 32p $CASE/dns.in)"
echo "inputs     : $(ls $CASE/*.in | wc -l) .in files (expected 8)"
