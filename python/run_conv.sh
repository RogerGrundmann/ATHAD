#!/bin/bash
# ATM_CONV_SCHEME: adjust (default since 2026-10-10) against tiedtke. 20 iterations from scratch, 8 threads (config_chk40.xml with nm = 20, checkpoint = 10);
# byte check that ATM_CONV_SCHEME=tiedtke reproduces the binary before the change (cli/had_pre), 2 iterations, 1 thread.
cd /home/roger/SynologyDrive/Cloudstation_Notebook/ATHAD/python
( cd run_ca_adj  && env OMP_NUM_THREADS=8 ../../cli/had config.xml > run.log 2>&1; echo $? > exit ) &
( cd run_ca_tied && env OMP_NUM_THREADS=8 ATM_CONV_SCHEME=tiedtke ../../cli/had config.xml > run.log 2>&1; echo $? > exit ) &
( cd run_ca_vpre && OMP_NUM_THREADS=1 ../../cli/had_pre config.xml > run.log 2>&1; echo $? > exit ) &
( cd run_ca_vt   && OMP_NUM_THREADS=1 ATM_CONV_SCHEME=tiedtke ../../cli/had config.xml > run.log 2>&1; echo $? > exit ) &
wait
cmpdir(){ n=0; s=0; for f in $(cd $1 && find . -type f); do n=$((n+1)); cmp -s $1/$f $2/$f && s=$((s+1)); done; echo "ATHAD: $1 = $2 : $s/$n files identical"; }
cmpdir run_ca_vpre/out run_ca_vt/out
cmpdir run_ca_tied/out run_ca_adj/out
echo "ATHAD exits: $(cat run_ca_*/exit | tr '\n' ' ')"
