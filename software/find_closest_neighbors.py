#!/usr/bin/env python3
import argparse
from ete3 import Tree

def extract_closest_species(newick_file, target_species, num_neighbors, output_file):
    """
    Extracts a subtree containing a target species and its N closest relatives
    based on phylogenetic distance.

    Args:
        newick_file (str): Path to the input Newick tree file.
        target_species (str): The name of the target species to find neighbors for.
                              (e.g., 'Homo_sapiens'). This is case-sensitive.
        num_neighbors (int): The number of closest species to extract.
        output_file (str): Path to save the resulting subtree in Newick format.
    """
    # Step 1: Load the main tree.
    print(f"Loading tree from '{newick_file}'...")
    try:
        main_tree = Tree(newick_file, format=1)
    except FileNotFoundError:
        print(f"Error: Newick file not found at '{newick_file}'")
        return
    except Exception as e:
        print(f"An error occurred while loading or parsing the Newick tree: {e}")
        return

    # Step 2: Find the target species leaf node in the tree.
    # The leaf name might contain extra info (e.g., 'Homo_sapiens:proteinID'),
    # so we search for a leaf name that starts with the target species name.
    target_node = None
    all_leaf_names = main_tree.get_leaf_names()
    
    # Find the full leaf name that corresponds to the target species
    matching_full_names = [name for name in all_leaf_names if name.startswith(target_species)]
    
    if not matching_full_names:
        print(f"Error: Target species '{target_species}' could not be found in the tree.")
        print("Please check for typos or name mismatches.")
        return
        
    # In case of multiple partial matches, we take the first one.
    # This assumes names like 'Homo_sapiens' and 'Homo_sapiens_neanderthalensis' are distinct enough.
    full_target_name = matching_full_names[0]
    target_node = main_tree.get_leaves_by_name(full_target_name)[0]
    print(f"Found target node: '{full_target_name}'")

    # Step 3: Calculate the distance from the target to all other leaves.
    distances = []
    print("Calculating distances to all other species...")
    for leaf in main_tree.iter_leaves():
        # No need to calculate distance to self, but doing so is harmless.
        distance = target_node.get_distance(leaf)
        distances.append((distance, leaf.name))

    # Step 4: Sort species by distance and select the closest N neighbors.
    distances.sort(key=lambda x: x[0])

    # The first item will be the target itself (distance 0). We need it plus N neighbors.
    total_to_keep = num_neighbors + 1
    
    if len(distances) < total_to_keep:
        print(f"Warning: The tree has only {len(distances)-1} other species, which is less than the requested {num_neighbors}.")
        print("The output tree will contain all species.")
        leaves_to_keep = [name for dist, name in distances]
    else:
        # Get the names of the target and its N closest neighbors.
        closest_species = distances[:total_to_keep]
        leaves_to_keep = [name for dist, name in closest_species]

    print(f"\nTop {num_neighbors} closest species to '{full_target_name}':")
    for dist, name in distances[1:total_to_keep]:
        print(f"- {name} (distance: {dist:.4f})")

    # Step 5: Prune the tree and save the result.
    main_tree.prune(leaves_to_keep, preserve_branch_length=True)
    main_tree.write(outfile=output_file, format=1)
    print(f"\nSuccess! Extracted subtree with {len(leaves_to_keep)} species saved to '{output_file}'")

def main():
    """Main function to parse command-line arguments and run the script."""
    parser = argparse.ArgumentParser(
        description="Extracts a subtree containing a target species and its N closest phylogenetic neighbors.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "-n", "--newick_file",
        required=True,
        help="Path to the input Newick tree file."
    )
    parser.add_argument(
        "-t", "--target_species",
        required=True,
        help="Name of the target species to find neighbors for (e.g., 'Homo_sapiens').\nThis is case-sensitive and should match the start of the leaf name in the tree."
    )
    parser.add_argument(
        "-k", "--num_neighbors",
        type=int,
        default=10,
        help="The number of closest neighbors to find (default: 10)."
    )
    parser.add_argument(
        "-o", "--output_file",
        required=True,
        help="Path for the output Newick subtree file (e.g., 'human_and_neighbors.newick')."
    )

    args = parser.parse_args()
    extract_closest_species(args.newick_file, args.target_species, args.num_neighbors, args.output_file)

if __name__ == "__main__":
    main()