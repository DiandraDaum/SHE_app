import csv
from ete3 import Tree
import argparse

def extract_subtree(csv_file, newick_file, taxonomic_group, output_file):
    """
    Extracts a subtree from a larger phylogenetic tree based on a taxonomic group.

    Args:
        csv_file (str): Path to the CSV file containing organism names and their
                        taxonomic lineage.
        newick_file (str): Path to the Newick file of the main phylogenetic tree.
        taxonomic_group (str): The name of the taxonomic group to extract (e.g., 'Primates').
                               This is case-sensitive.
        output_file (str): Path to save the resulting subtree in Newick format.
    """
    # Step 1: Find all organisms belonging to the specified taxonomic group from the CSV.
    organisms_to_keep = set()
    print(f"Searching for organisms in the '{taxonomic_group}' group...")
    try:
        with open(csv_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            # Check for required columns in the CSV header
            if 'Organism_formatted' not in reader.fieldnames or 'Taxonomic lineage' not in reader.fieldnames:
                print("Error: CSV file must contain 'Organism_formatted' and 'Taxonomic lineage' columns.")
                return

            for row in reader:
                # This checks if the taxonomic group string is part of the full lineage string.
                if taxonomic_group in row['Taxonomic lineage']:
                    organisms_to_keep.add(row['Organism_formatted'])
    except FileNotFoundError:
        print(f"Error: CSV file not found at '{csv_file}'")
        return
    except Exception as e:
        print(f"An error occurred while reading the CSV file: {e}")
        return

    if not organisms_to_keep:
        print(f"No organisms found for the taxonomic group '{taxonomic_group}'. No tree will be generated.")
        return

    print(f"Found {len(organisms_to_keep)} organisms in the '{taxonomic_group}' group.")

    # Step 2: Load the main tree and prune it to keep only the desired organisms.
    print(f"Loading tree from '{newick_file}'...")
    try:
        main_tree = Tree(newick_file, format=1)
    except FileNotFoundError:
        print(f"Error: Newick file not found at '{newick_file}'")
        return
    except Exception as e: # Catches potential ete3 parsing errors
        print(f"An error occurred while loading or parsing the Newick tree: {e}")
        return
        
    # The 'prune' method from ete3 efficiently keeps only the specified leaves.
    # We must match the organism names from the CSV to the full leaf names in the tree,
    # which might contain extra info (e.g., 'Organism_Name:ProteinID').
    
    # Create a map from the simple organism name to the full leaf name in the tree.
    all_leaf_names = main_tree.get_leaf_names()
    org_name_map = {name.split(':')[0] if ':' in name else name: name for name in all_leaf_names}
    
    leaves_in_tree_to_keep = []
    for org in organisms_to_keep:
        if org in org_name_map:
            leaves_in_tree_to_keep.append(org_name_map[org])

    if not leaves_in_tree_to_keep:
        print(f"Warning: None of the {len(organisms_to_keep)} organisms from the CSV were found as leaves in the tree file.")
        print("Please check for inconsistencies in organism names between your files.")
        return

    print(f"Pruning tree to keep {len(leaves_in_tree_to_keep)} matching leaves.")
    main_tree.prune(leaves_in_tree_to_keep, preserve_branch_length=True)

    # Step 3: Save the new, smaller subtree.
    if len(main_tree) > 0:
        main_tree.write(outfile=output_file, format=1)
        print(f"\nSuccess! Extracted subtree saved to '{output_file}'")
    else:
        print("The resulting subtree has no leaves. No file was written.")

def main():
    """Main function to extract a subtree from a phylogenetic tree based on a taxonomic group."""
    parser = argparse.ArgumentParser(
        description="Extract a subtree from a phylogenetic tree based on a taxonomic group.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "-c", "--csv_file", 
        required=True, 
        help="Path to the input CSV file.\nMust contain 'Organism_formatted' and 'Taxonomic lineage' columns."
    )
    parser.add_argument(
        "-n", "--newick_file", 
        required=True, 
        help="Path to the input Newick tree file."
    )
    parser.add_argument(
        "-t", "--taxonomic_group", 
        required=True, 
        help="Name of the taxonomic group to extract (e.g., 'Primates', 'Fungi').\nThis is case-sensitive."
    )
    parser.add_argument(
        "-o", "--output_file", 
        required=True, 
        help="Path for the output Newick subtree file (e.g., 'primates_tree.newick')."
    )

    args = parser.parse_args()
    extract_subtree(args.csv_file, args.newick_file, args.taxonomic_group, args.output_file)

if __name__ == "__main__":
    main()
