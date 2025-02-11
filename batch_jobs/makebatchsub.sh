#!/bin/bash
#file="folders_to_submit.txt"
#INPUT FILE IS A TXT FILE WITH A LIST OF ALL PATHS TO YAML CONFIGS WHERE YAML FILES ARE SAVED AS "Selections_NDGAr_acceptancecorrection_2dhist_<number>.yaml" and "SamplePDFDuneNDGAr_FHC_CCnumuselec_<number>.yaml"
# run e.g "ls -d /vols/dune/nk3717/MaCh3_EL9/MaCh3_DUNE/configs/*_Bfield0_5_instrumentedrad249_45_withecalcontainment &> $file.txt" to pipe filepaths to a txt file
file="folders_to_submit_all249_45_Bfield0_2.txt"
#ls -d /vols/dune/nk3717/MaCh3_EL9/MaCh3_DUNE/configs/*_Bfield0_5_instrumentedrad249_45_withecalcontainment &> $file
linenum=1
outputnum=1
output="/vols/dune/nk3717/outlogs_1/out/out.$(CLUSTER).$(PROCESS)"
error="/vols/dune/nk3717/outlogs_1/err/err.$(CLUSTER).$(PROCESS)"
log="/vols/dune/nk3717/outlogs_1/log/log.$(CLUSTER)"
arguments="$(CLUSTER) $(PROCESS)"
while IFS= read -r line; do
  # Process each line here
  echo "Line: $line"
  script_name="submit_$linenum.sh"
  sub_name="submit_$linenum.sub"
  if (( $linenum % 10 == 0 )); then
    outputnum=$(($outputnum+1))
    sleep 4000
  fi
  #new_string="${line//configs/outputs}"
  #new_string1="${new_string//SamplePDFConfigs_/}"
  cat << EOF > "$script_name"
#!/bin/bash
cd /vols/dune/nk3717/MaCh3_EL9/MaCh3_DUNE/
source setup.sh
Selections $line/Selections_NDGAr_acceptancecorrection_2dhist_\$2.yaml $line/SamplePDFDuneNDGAr_FHC_CCnumuselec_\$2.yaml
  
EOF
if [ ! -d "/vols/dune/nk3717/outlogs_$outputnum" ]
then
  mkdir -p /vols/dune/nk3717/outlogs_$outputnum/out
  mkdir -p /vols/dune/nk3717/outlogs_$outputnum/err
  mkdir -p /vols/dune/nk3717/outlogs_$outputnum/log
fi

  cat << EOF > "$sub_name"
executable = $script_name
output = /vols/dune/nk3717/outlogs_$outputnum/out/out.\$(CLUSTER).\$(PROCESS)
error = /vols/dune/nk3717/outlogs_$outputnum/err/err.\$(CLUSTER).\$(PROCESS)
log = /vols/dune/nk3717/outlogs_$outputnum/log/log.\$(CLUSTER)
arguments = \$(CLUSTER) \$(PROCESS)
+MaxRuntime = 5400
queue 24
EOF
#NOTE CHANGE MAX RUNTIME TO 4500 for 249.45 instrumented regions and 2700 for 199.45 radius
  chmod +x $script_name
  condor_submit $sub_name
  sleep 30
linenum=$(($linenum+1))
done < "$file"
