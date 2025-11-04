import pandas as pd
import argparse

def get_species_id(protein_name):
    """
    Parses a protein name to find a UniProt proteome ID.
    Returns the ID (e.g., 'UP000005447_10141') if found, otherwise returns None.
    """
    if isinstance(protein_name, str) and protein_name.startswith("UP"):
        # We need the full ID including taxon for uniqueness across species
        return protein_name.split('-')[0]
    return None

def find_core_families(families_file, presence_threshold=1.0):
    """
    Identifies protein families present in a certain percentage of species.
    
    Args:
        families_file (str): Path to the protein_families.csv file.
        presence_threshold (float): The minimum fraction of species a family must be in 
                                    to be considered 'core' (e.g., 1.0 for 100%, 0.99 for 99%).
    """
    print(f"--- Loading '{families_file}' ---")
    try:
        family_df = pd.read_csv(families_file)
    except FileNotFoundError:
        print(f"FATAL ERROR: The file '{families_file}' was not found.")
        return

    print("--- Parsing species IDs and filtering data ---")
    family_df['species_id'] = family_df['protein'].apply(get_species_id)
    family_df.dropna(subset=['species_id'], inplace=True)

    if family_df.empty:
        print("FATAL ERROR: No valid 'UP...' protein IDs found in the file.")
        return
        
    # Get the total number of unique species in the dataset
    total_species = family_df['species_id'].nunique()
    print(f"Found a total of {total_species} unique species in the dataset.")
    
    # For each family, count the number of unique species it appears in
    print("--- Counting species per family ---")
    species_per_family = family_df.groupby('family_id')['species_id'].nunique()
    
    # Calculate the presence fraction for each family
    presence_fraction = species_per_family / total_species
    
    # Filter for families that meet the presence threshold
    threshold_value = int(presence_threshold * 100)
    print(f"--- Identifying core families present in >{threshold_value}% of species ---")
    
    core_families = presence_fraction[presence_fraction >= presence_threshold]
    
    if core_families.empty:
        print(f"No families were found present in at least {threshold_value}% of species.")
        print("You might want to try a lower --threshold, for example: 0.95")
        return

    # Sort the core families by their prevalence
    core_families = core_families.sort_values(ascending=False)
    
    output_file = f"core_proteome_{threshold_value}percent.txt"
    core_families.index.to_series().to_csv(output_file, index=False, header=False)
    
    print(f"✅ Success! Found {len(core_families)} core protein families.")
    print(f"The list of core family IDs has been saved to '{output_file}'")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find the core set of protein families across a dataset.")
    parser.add_argument("--families-file", required=True, help="Path to 'protein_families.csv'.")
    parser.add_argument(
        "--threshold", 
        type=float, 
        default=0.99, 
        help="Fraction of species a family must be present in to be considered 'core' (default: 0.99 for 99%%)."
    )
    
    args = parser.parse_args()
    find_core_families(args.families_file, args.threshold)
