import pandas as pd

# --- Configuration ---
# 1. Replace these with your actual filenames
protein_csv_file = '.../.../protein_families_238.csv' # Your main data file
id_text_file = 'core_proteome_95percent.txt'       # Your text file with IDs like 'fam_0', 'fam_7'
output_csv_file = '.../.../filtered_proteins.csv'     # The name for the new output file

# --- Main Script ---
try:
    # 2. Read the list of family IDs from the text file
    # We use .strip() to remove any extra spaces or newline characters
    with open(id_text_file, 'r') as f:
        ids_to_keep = [line.strip() for line in f]

    # 3. Read the main protein CSV data into a pandas DataFrame
    df = pd.read_csv(protein_csv_file)

    # 4. Filter the DataFrame to keep only rows where 'family_id' is in our list
    # The .isin() method is a very efficient way to perform this check.
    filtered_df = df[df['family_id'].isin(ids_to_keep)]

    # 5. Save the filtered results to a new CSV file
    # `index=False` prevents pandas from writing the row numbers into the file
    filtered_df.to_csv(output_csv_file, index=False)

    print(f"✅ Filtering complete!")
    print(f"Found {len(filtered_df)} matching rows.")
    print(f"Results saved to '{output_csv_file}'")

except FileNotFoundError as e:
    print(f"❌ Error: The file '{e.filename}' was not found. Please make sure it's in the same directory as the script or provide the full path.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")