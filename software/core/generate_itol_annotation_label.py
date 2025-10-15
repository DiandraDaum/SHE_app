import argparse
import sys
import pandas as pd
from ete3 import Tree

def get_species_id(protein_name):
    """Robustly parses a protein name to find a UniProt proteome ID."""
    if isinstance(protein_name, str) and protein_name.startswith("UP"):
        return protein_name.split('-')[0].split('_')[0]
    return None

def create_itol_file(families_file, names_file, newick_tree_file, family_to_visualize, output_file):
    """
    Generates an iTOL annotation file to color the tree leaf labels based on protein family presence.
    """
    # 1. Load the names file and create the name map
    print(f"--- Step 1: Loading Name Map from '{names_file}' ---")
    try:
        names_df = pd.read_csv(names_file)
    except FileNotFoundError:
        print(f"FATAL ERROR: The names file '{names_file}' was not found.")
        sys.exit(1)

    if 'ID' not in names_df.columns or 'Organism_formatted' not in names_df.columns:
        print(f"FATAL ERROR: The required columns 'ID' and/or 'Organism_formatted' were not found in '{names_file}'.")
        sys.exit(1)

    mapping_df = names_df[['ID', 'Organism_formatted']].drop_duplicates().set_index('ID')
    id_to_name_map = mapping_df['Organism_formatted'].to_dict()
    print(f"Successfully created a name map for {len(id_to_name_map)} unique IDs.")

    # 2. Load and process protein family data
    print(f"--- Step 2: Preparing data for '{family_to_visualize}' from '{families_file}' ---")
    try:
        family_df = pd.read_csv(families_file)
    except FileNotFoundError:
        print(f"FATAL ERROR: The families file '{families_file}' was not found.")
        sys.exit(1)

    family_df['species_id'] = family_df['protein'].apply(get_species_id)
    family_df.dropna(subset=['species_id'], inplace=True)

    if len(family_df) == 0:
        print("FATAL ERROR: No proteins with 'UP...' IDs were found in the data file.")
        sys.exit(1)

    presence_absence = pd.crosstab(family_df['family_id'], family_df['species_id']) > 0
    if family_to_visualize not in presence_absence.index:
        print(f"FATAL ERROR: Family '{family_to_visualize}' not found in the filtered dataset.")
        sys.exit(1)
    family_pattern = presence_absence.loc[family_to_visualize]

    # 3. Load the tree to get all leaf names
    print(f"--- Step 3: Loading Tree from '{newick_tree_file}' ---")
    try:
        t = Tree(newick_tree_file)
    except FileNotFoundError:
        print(f"FATAL ERROR: The tree file '{newick_tree_file}' was not found.")
        sys.exit(1)
    all_leaf_ids = t.get_leaf_names()
    print(f"Found {len(all_leaf_ids)} leaves in the tree.")

    # 4. Write the iTOL annotation file for coloring leaf labels
    print(f"--- Step 4: Writing iTOL label style file to '{output_file}' ---")

    present_color = "#ff0000"  # Red
    absent_color = "#0000ff"   # Blue
    
    with open(output_file, 'w') as f:
        f.write("DATASET_STYLE\n")
        f.write("SEPARATOR COMMA\n")
        
        # Dataset label
        f.write(f"DATASET_LABEL,{family_to_visualize}_presence_labels\n")
        
        # Legend definition
        f.write(f"LEGEND_TITLE,Presence of {family_to_visualize}\n")
        f.write("LEGEND_SHAPES,1,1\n") # 1 = square shape
        f.write(f"LEGEND_COLORS,{present_color},{absent_color}\n")
        f.write("LEGEND_LABELS,Present,Absent\n")
        f.write("\n")
        
        # Data section
        f.write("DATA\n")

        # Write data for each leaf
        for leaf_id in all_leaf_ids:
            output_name = id_to_name_map.get(leaf_id, leaf_id)
            species_id = leaf_id.split('_')[0].strip()
            is_present = family_pattern.get(species_id, False)
            color = present_color if is_present else absent_color

            # --- CHANGE: Modified 'branch' to 'label' to color the leaf text ---
            # This tells iTOL to modify the 'color' style of the 'label' (the text)
            f.write(f"{output_name},label,color,{color}\n")

    # --- CHANGE: Updated the success message ---
    print(f"✅ Success! iTOL annotation file for label styling saved to '{output_file}'.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate an iTOL annotation file to color tree leaf labels based on protein family presence.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "--families-file",
        required=True,
        help="Path to the CSV file containing protein families (e.g., protein_families.csv)."
    )
    parser.add_argument(
        "--names-file",
        required=True,
        help="Path to the CSV file for name mapping (e.g., names.csv)."
    )
    parser.add_argument(
        "--newick-tree-file",
        required=True,
        help="Path to your phylogenetic tree in Newick format."
    )
    parser.add_argument(
        "--family-to-visualize",
        required=True,
        help="The family ID you want to visualize (e.g., 'family_123')."
    )
    parser.add_argument(
        "--output-file",
        required=True,
        # --- CHANGE: Updated the example output filename ---
        help="Path to save the output iTOL annotation file (e.g., 'itol_label_styles.txt')."
    )

    args = parser.parse_args()

    create_itol_file(
        args.families_file,
        args.names_file,
        args.newick_tree_file,
        args.family_to_visualize,
        args.output_file
    )