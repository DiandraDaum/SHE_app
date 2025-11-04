#!/bin/bash

# Assign arguments to variables
COLOR_MAPPING_CSV=/Users/diandra/tree/complete_tree/color_mapping_all.csv
INPUT_TXT=/Users/diandra/tree/complete_tree/itol_annotation.txt
OUTPUT_TXT=/Users/diandra/tree/complete_tree/iTOL_evolrate_all.txt

# Run the Python script
python3 evol_rates.py "$COLOR_MAPPING_CSV" "$INPUT_TXT" "$OUTPUT_TXT"