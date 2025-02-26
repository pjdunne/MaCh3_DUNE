#!/bin/bash

# Input YAML files

EVENTRATES_CONFIG="/scratch/abipeake/MaCh3_DUNE_FD/MaCh3_DUNE/configs/EventRates_BeamVD.yaml"

INPUT_FILE="/scratch/abipeake/MaCh3_DUNE_FD/MaCh3_DUNE/configs/Samples/FD_VD/FD_VD_FHC.yaml"
OUTPUT_PREFIX="modified_"
EVENTRATES_CONFIG="/scratch/abipeake/MaCh3_DUNE_FD/MaCh3_DUNE/configs/EventRates_BeamVD.yaml"

ROOT_FILE="efficiency_purity_results_0.01steps_new_2nue.root"
ROOT_SCRIPT="store_results.C"
CSV_FILE="efficiency_purity_resultsnew_0.01steps_new_2nue.csv"

# Initialize the ROOT file
echo "void store_results() { TFile f(\"$ROOT_FILE\", \"RECREATE\"); TTree t(\"EfficiencyPurity\", \"Efficiency and Purity Data\"); }" > "$ROOT_SCRIPT"
root -l -b -q "$ROOT_SCRIPT"

# Loop through values from 0.4 to 0.7 in increments of 0.1
for nue_cut in $(seq 0.35 0.1 1.0); do
    OUTPUT_FILE="configs/Samples/test/${OUTPUT_PREFIX}${nue_cut}_nue.yaml"

    cp $INPUT_FILE $OUTPUT_FILE

    # Modify YAML file with new nue_cut value
    sed -i "s|^nue_cut.*|nue_cut: $nue_cut|g" $OUTPUT_FILE
    #sed "s/^ nue_cut: .*/nue_cut: $nue_cut/" "$INPUT_FILE" > "$OUTPUT_FILE"
    echo "Created: $OUTPUT_FILE with nue_cut = $nue_cut"

    # Update EventRates config to use the modified YAML file
    #sed "s|DUNESamples: \[.*\]|DUNESamples: [\"$OUTPUT_FILE\"]|" "$EVENTRATES_CONFIG" > "updated_eventrates.yaml"
    #sed -i "s|DUNESamples.*|DUNESamples: ["$OUTPUT_FILE"]|g"  "$EVENTRATES_CONFIG"
    sed -i "s|DUNESamples.*|DUNESamples: [\"$OUTPUT_FILE\"]|g" "$EVENTRATES_CONFIG"
    echo "Updated EventRates config with $OUTPUT_FILE"

    grep "DUNESamples" "$EVENTRATES_CONFIG"


    # Debugging: Check nue_cut value in modified YAML
    grep "nue_cut" "$OUTPUT_FILE"
    echo "DEBUG: EVENTRATES_CONFIG='$EVENTRATES_CONFIG'"

    # Run EventRates and check for errors
    #OUTPUT=$(./build/src/EventRates updated_eventrates.yaml 2>&1)
    OUTPUT=$(./build/src/EventRates "$EVENTRATES_CONFIG" 2>&1)

    if [[ $? -ne 0 ]]; then
        echo "ERROR: EventRates execution failed for nue_cut=$nue_cut"
        continue
    fi

    echo "DEBUG: Full EventRates Output"
    echo "$OUTPUT"

    # Extract relevant values
    eventsthatpasscut=$(echo "$OUTPUT" | grep "eventsthatpassedcut =" | awk -F '= ' '{print $2}' | tail -n1)
    total_true_ccnue=$(echo "$OUTPUT" | grep "total_true_ccnue  =" | awk -F '= ' '{print $2}' | tail -n1)
    total_events_incut=$(echo "$OUTPUT" | grep "total_events_incut =" | awk -F '= ' '{print $2}' | tail -n1)

    echo "Extracted values: eventsthatpasscut=$eventsthatpasscut, total_true_ccnue=$total_true_ccnue, total_events_incut=$total_events_incut"

    # Compute efficiency and purity (handle division by zero)
    efficiency=0
    purity=0
    if [[ -n "$total_true_ccnue" && $(echo "$total_true_ccnue > 0" | bc -l) -eq 1 ]]; then
        efficiency=$(echo "scale=5; $eventsthatpasscut / $total_true_ccnue" | bc)
    fi
    if [[ -n "$total_events_incut" && $(echo "$total_events_incut > 0" | bc -l) -eq 1 ]]; then
        purity=$(echo "scale=5; $eventsthatpasscut / $total_events_incut" | bc)
    fi

    # Store the results in CSV
    echo "$nue_cut,$efficiency,$purity" >> "$CSV_FILE"
    echo "Saved efficiency ($efficiency) and purity ($purity) for nue_cut = $nue_cut"

    # Append results to ROOT script
    echo "
    void store_results() {
        TFile f(\"$ROOT_FILE\", \"UPDATE\");
        TTree *t = (TTree*) f.Get(\"EfficiencyPurity\");
        if (!t) { t = new TTree(\"EfficiencyPurity\", \"Efficiency and Purity Data\"); }
        double nue_cut = $nue_cut, efficiency = $efficiency, purity = $purity;
        t->Branch(\"nue_cut\", &nue_cut, \"nue_cut/D\");
        t->Branch(\"efficiency\", &efficiency, \"efficiency/D\");
        t->Branch(\"purity\", &purity, \"purity/D\");
        t->Fill();
        t->Write(\"\", TObject::kOverwrite);
        f.Close();
    }" > "$ROOT_SCRIPT"

    # Run the ROOT script to store data
    root -l -b -q "$ROOT_SCRIPT"

    echo "Results saved for nue_cut = $nue_cut"
    #rm -f updated_eventrates.yaml
done

echo "All calculations complete. Results stored in $ROOT_FILE."