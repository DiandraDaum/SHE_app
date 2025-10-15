import argparse
import pandas as pd
import networkx as nx
from tqdm import tqdm

def create_protein_families(raw_hits_file, output_file):
    print("Loading raw homologous pairs...")
    hits_df = pd.read_csv(raw_hits_file)
    required_columns = {"protein1", "protein2", "bit_score"}
    missing = required_columns - set(hits_df.columns)
    if missing:
        raise ValueError(
            f"Input file is missing required columns: {', '.join(sorted(missing))}. "
            "Expected columns: protein1, protein2, bit_score"
        )
    print("Building protein graph...")
    G = nx.from_pandas_edgelist(hits_df, 'protein1', 'protein2', ['bit_score'])
    print("Finding connected components (protein families)...")
    components = nx.connected_components(G)
    protein_to_family = {}
    for i, component in enumerate(tqdm(components, desc="Assigning families")):
        family_id = f"fam_{i}"
        for protein in component:
            protein_to_family[protein] = family_id
    family_df = pd.DataFrame(protein_to_family.items(), columns=['protein', 'family_id'])
    family_df.to_csv(output_file, index=False)
    print(f"Protein family assignments saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Create protein families from homologous pairs.")
    parser.add_argument("--raw-hits", required=True, dest="raw_hits",
                        help="Path to CSV with columns: protein1, protein2, bit_score")
    parser.add_argument("--output", required=True,
                        help="Path to write protein_families.csv")
    args = parser.parse_args()
    create_protein_families(args.raw_hits, args.output)

if __name__ == "__main__":
    main()