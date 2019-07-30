#!/bin/bash

source ~/work/tools/my_env.sh

script="../experiment_figures.py"
datadir="../../../data"

tannic_acid_msd_outdir="$datadir/msd/TannicAcid"
olaparib_msd_outdir="$datadir/msd/Olaparib"

log_message "Tannic Acid + R5020"
python3 $script "Tannic-Acid+R5020" $tannic_acid_msd_outdir\
  $datadir/0_Tannic_Acid_6h/1_Tannic_R5020_6h/Sample1 \
  $datadir/0_Tannic_Acid_6h/1_Tannic_R5020_6h/Sample2 \

log_message "Tannic Acid Control"
python3 $script "Tannic-Acid-Control" $tannic_acid_msd_outdir\
  $datadir/0_Tannic_Acid_6h/2_Control_R5020_6h

log_message "Olaparib"
python3 $script "Olaparib+R5020" $olaparib_msd_outdir\
  $datadir/1_Olaparib_R5020

log_message "Olaparib Control"
python3 $script "Olaparib-R5020-Control" $olaparib_msd_outdir\
  $datadir/2_DMSO_R5020_Control

log_message "EtOH Control"
python3 $script "EtOH-NoHormone-Control" $olaparib_msd_outdir\
  $datadir/3_EtOH_Nohormone_Control

log_message "Olaparib Subset: with Olaparib"
python3 $script "Olaparib-Subset+R5020" $olaparib_msd_outdir\
  $datadir/Subset_Trajectories_Olaparib/1_Olaparib_R5020

log_message "Olaparib Subset: R5020 Control"
python3 $script "Olaparib-Subset-R5020-Control" $olaparib_msd_outdir\
  $datadir/Subset_Trajectories_Olaparib/2_DMSO_R5020_Control

log_message "Olaparib Subset: EtOH Control"
python3 $script "Olaparib-Subset-EtOH-Control" $olaparib_msd_outdir\
  $datadir/Subset_Trajectories_Olaparib/3_EtOH_Nohormone_Control

log_message "Tannic Acid Subset: with Tannic Acid"
python3 $script "Tannic-Acid-Subset+R5020" $tannic_acid_msd_outdir \
  $datadir/Subset_Trajectories_TannicAcid/1_Tannic_R5020_6h/Sample1 \
  $datadir/Subset_Trajectories_TannicAcid/1_Tannic_R5020_6h/Sample2

log_message "Tannic Acid Subset: Control"
python3 $script "Tannic-Acid-Subset-R5020-Control" $tannic_acid_msd_outdir\
  $datadir/Subset_Trajectories_TannicAcid/2_Control_R5020_6h
