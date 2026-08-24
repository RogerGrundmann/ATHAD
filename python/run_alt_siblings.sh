#!/bin/bash
# ATM_CELL_ALTERNATE A/B in ATHAD and ATHAD_COND. 4 runs x 40 iterations, 24 threads.
# ATHAD ~4 min/arm (im=41), ATHAD_COND ~7 min/arm (im=61)  ->  ~22 min total.
export OMP_NUM_THREADS=24
log=/tmp/claude-1000/-home-roger-SynologyDrive-Cloudstation-Notebook-ATHAD/1e34b8f4-a0f1-49fb-bfdd-40702fab034e/scratchpad
B=/home/roger/SynologyDrive/Cloudstation_Notebook
cd $B/ATHAD/python      || exit 1
echo "### ATHAD off $(date +%T)"; ATM_CELL_ALTERNATE=0 ../cli/had  config_alt_off.xml > $log/ALT_ATHAD_off.log 2>&1; echo " rc=$?"
echo "### ATHAD on  $(date +%T)";                      ../cli/had  config_alt_on.xml  > $log/ALT_ATHAD_on.log  2>&1; echo " rc=$?"
cd $B/ATHAD_COND/python || exit 1
echo "### COND off  $(date +%T)"; ATM_CELL_ALTERNATE=0 ../cli/cond config_alt_off.xml > $log/ALT_COND_off.log 2>&1; echo " rc=$?"
echo "### COND on   $(date +%T)";                      ../cli/cond config_alt_on.xml  > $log/ALT_COND_on.log  2>&1; echo " rc=$?"
echo "### ALT SIBLINGS DONE $(date +%T)"
