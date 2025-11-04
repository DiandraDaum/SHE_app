import pandas as pd
import sys

# Check for correct number of arguments
if len(sys.argv) != 4:
    print("Usage: python update_colors.py <color_mapping_csv> <input_txt> <output_txt>")
    sys.exit(1)

# Load the CSV file containing Branch and Color
csv_file_path = sys.argv[1]
color_mapping = pd.read_csv(csv_file_path)

# Create a dictionary for quick lookup of colors by branch
color_dict = dict(zip(color_mapping['Branch'], color_mapping['Color']))

# Read the original text file
txt_file_path = sys.argv[2]
with open(txt_file_path, 'r') as file:
    lines = file.readlines()

# Update colors in the text file
updated_lines = []
for line in lines:
    #if line.startswith("unnamed_branch"):  # Check if the line contains a branch
    parts = line.split('\t')  # Split the line by tab
    branch_name = parts[0]  # Get the branch name
    if branch_name in color_dict:  # Check if branch is in the color mapping
        parts[2] = color_dict[branch_name]  # Update the color
        parts[1] = "branch	color"
        parts[3] = "width	5\n"
        updated_lines.append('\t'.join(parts))  # Rejoin the parts
        
     

# Write the updated lines to a new text file
new_txt_file_path = sys.argv[3]
with open(new_txt_file_path, 'w') as file:
    file.writelines(updated_lines)

print(f"Updated data has been written to {new_txt_file_path}")