for j in {170..230..20}
do
  for i in {50..100..10}
  do
    for k in {50..100..10}
    do
      output="Selections_CC_AcceptanceCorrection_2dhist_AcceptedEvents_TrueRad_${j}_pi0eff_${i}_gammaeff_${k}.root"
      output2="Selections_CC_AcceptanceCorrection_2dhist_AllEvents_TrueRad_${j}_pi0eff_${i}_gammaeff_${k}.root"
      line="$(grep -n 'FileName:' configs/Selections_NDGAr_acceptancecorrection_2dhist_AcceptedEvents.yaml | cut -d ':' -f1)"
      echo $line
      sed -i "${line}s/.*/    FileName: outputs\/AcceptanceCorrection_momresthreshold_5percent_pixelgrid\/"${output}"/" configs/Selections_NDGAr_acceptancecorrection_2dhist_AcceptedEvents.yaml
      line2="$(grep -n 'IsAccepted' configs/SamplePDFDuneNDGAr_FHC_CCnumuselec.yaml | cut -d ':' -f1)"
      line3=$(expr $line2 + 1)
      sed -i "${line3}s/.*/      Bounds: \[1.0, 2.0\]/" configs/SamplePDFDuneNDGAr_FHC_CCnumuselec.yaml
      line6="$(grep -n 'TrueRad' configs/SamplePDFDuneNDGAr_FHC_CCnumuselec.yaml | cut -d ':' -f1)"
      line7=$(expr $line6 + 1)
      sed -i "${line7}s/.*/      Bounds: \[0.0, ${j}.0\]/" configs/SamplePDFDuneNDGAr_FHC_CCnumuselec.yaml
      line4="$(grep -n 'pi0_reco_efficiency' configs/SamplePDFDuneNDGAr_FHC_CCnumuselec.yaml | cut -d ':' -f1)"
      sed -i "${line4}s/.*/  pi0_reco_efficiency: 0.${i}/" configs/SamplePDFDuneNDGAr_FHC_CCnumuselec.yaml
      if [$i -eq 100]
       then
         sed -i "${line4}s/.*/  pi0_reco_efficiency: 1.0/" configs/SamplePDFDuneNDGAr_FHC_CCnumuselec.yaml
      fi
      line5="$(grep -n 'gamma_reco_efficiency' configs/SamplePDFDuneNDGAr_FHC_CCnumuselec.yaml | cut -d ':' -f1)"
      sed -i "${line5}s/.*/  gamma_reco_efficiency: 0.${k}/" configs/SamplePDFDuneNDGAr_FHC_CCnumuselec.yaml
      if [$k -eq 100]
       then
         sed -i "${line5}s/.*/  gamma_reco_efficiency: 1.0/" configs/SamplePDFDuneNDGAr_FHC_CCnumuselec.yaml
      fi
      linetwod="$(grep -n 'FileName:' configs/Selections_NDGAr_acceptancecorrection_2dhist_AllEvents.yaml | cut -d ':' -f1)"
      echo $linetwod
      sed -i "${linetwod}s/.*/    FileName: outputs\/AcceptanceCorrection_momresthreshold_5percent_pixelgrid\/"${output2}"/" configs/Selections_NDGAr_acceptancecorrection_2dhist_AllEvents.yaml

      Selections configs/Selections_NDGAr_acceptancecorrection_2dhist_AllEvents.yaml &> acceptancecorr_allevents.txt &
      wait
      line8="$(grep -n 'IsAccepted' configs/SamplePDFDuneNDGAr_FHC_CCnumuselec.yaml | cut -d ':' -f1)"
      line9=$(expr $line8 + 1)
      sed -i "${line9}s/.*/      Bounds: \[0.0, 2.0\]/" configs/SamplePDFDuneNDGAr_FHC_CCnumuselec.yaml

      Selections configs/Selections_NDGAr_acceptancecorrection_2dhist_AcceptedEvents.yaml &> acceptancecorr_acceptedevents.txt &
      wait
    done
  done
done


