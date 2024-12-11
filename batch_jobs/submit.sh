#!/bin/bash

cd /vols/dune/nk3717/MaCh3_EL9/MaCh3_DUNE/
source setup.sh
#read input1
#read input2
#echo inputfiles1/${ProcId}.txt
#cat inputfiles1/$2.txt inputfiles2/$2.txt &> output_test.txt

#echo "/vols/dune/nk3717/MaCh3_EL9/MaCh3_DUNE/configs/SamplePDFConfigs_fidrad160_0_Bfield0_5_withecalcontainment/SamplePDFDuneNDGAr_FHC_CCnumuselec_$2.yaml"

Selections /vols/dune/nk3717/MaCh3_EL9/MaCh3_DUNE/configs/SamplePDFConfigs_fidrad160_0_Bfield0_5_withecalcontainment/Selections_NDGAr_acceptancecorrection_2dhist_$2.yaml /vols/dune/nk3717/MaCh3_EL9/MaCh3_DUNE/configs/SamplePDFConfigs_fidrad160_0_Bfield0_5_withecalcontainment/SamplePDFDuneNDGAr_FHC_CCnumuselec_$2.yaml
