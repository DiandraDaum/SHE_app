import os
import pandas as pd
import subprocess
import shutil
import glob
import argparse
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

parser = argparse.ArgumentParser(description="Builds full protein IDs from a 'one-archive-per-organism' directory structure.")
parser.add_argument("--root-dir", required=True, help="Root directory containing the organism .tar.zst archives.")
parser.add_argument("--output-file", required=True, help="Output CSV file for all reconstructed homologous pairs.")
parser.add_argument("--bit-score-threshold", type=float, default=238.0, help="Minimum bit score for a hit.")
parser.add_argument("--e-value-threshold", type=float, default=0, help="Maximum e-value for a hit.")
parser.add_argument("--zstd-path", default="/proj/berzelius-2023-267/users/x_patbr/software/zstd/zstd", help="Path to zstd program.")
parser.add_argument("--temp-dir", default=None, help="Directory for temporary files. Defaults to the current working directory.")

def find_tar_command():
    for cmd in ['gtar', 'tar']:
        if shutil.which(cmd):
            return cmd
    raise FileNotFoundError("Neither 'gtar' nor 'tar' found in PATH.")

def process_archive(tar_file_path, writer, args):
    tar_cmd = find_tar_command()

    # --- LOGIC MODIFIED TO USE --temp-dir ---
    unique_dir_name = f'temp_extract_{os.getpid()}'
    if args.temp_dir:
        # Ensure the user-specified base temp directory exists
        os.makedirs(args.temp_dir, exist_ok=True)
        temp_dir = os.path.join(args.temp_dir, unique_dir_name)
    else:
        # Fallback to the original behavior (current working directory)
        temp_dir = unique_dir_name
    
    # This part remains the same, but 'temp_dir' is now correctly located
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    os.makedirs(temp_dir)

    try:
        archive_basename = os.path.basename(tar_file_path)
        organism_A_id = archive_basename.replace('.tar.zst', '')

        logging.info(f"Processing archive for organism: {organism_A_id}")
        
        subprocess.run(
            [tar_cmd, '--use-compress-program', args.zstd_path, '-xvf', tar_file_path, '-C', temp_dir],
            check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

        csv_files = glob.glob(os.path.join(temp_dir, '**', '*.csv'), recursive=True)
        for csv_file in csv_files:
            try:
                csv_basename = os.path.basename(csv_file)
                
                base_name_no_ext = csv_basename.replace('_top_hits.csv', '')
                organism_B_id = base_name_no_ext.split('vs')[-1]

                for chunk in pd.read_csv(csv_file, header=None, usecols=[0, 1, 3, 4], chunksize=10000):
                    chunk.columns = ['proteinA', 'proteinB', 'e_value', 'bit_score']
                    
                    filtered = chunk[
                        (chunk['bit_score'] >= args.bit_score_threshold) &
                        (chunk['e_value'] <= args.e_value_threshold)
                    ].copy()
                    
                    if not filtered.empty:
                        filtered['full_protein_A'] = organism_A_id + '-' + filtered['proteinA'].astype(str)
                        filtered['full_protein_B'] = organism_B_id + '-' + filtered['proteinB'].astype(str)
                        
                        output_chunk = filtered[['full_protein_A', 'full_protein_B', 'bit_score']]
                        output_chunk.columns = ['protein1', 'protein2', 'bit_score']
                        output_chunk.to_csv(writer, mode='a', header=False, index=False)
            except Exception as e:
                logging.error(f"Error processing CSV {csv_file}: {e}")
    
    except Exception as e:
        logging.error(f"Error processing archive {tar_file_path}: {e}")
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

def main():
    args = parser.parse_args()
    
    with open(args.output_file, 'w') as writer:
        writer.write('protein1,protein2,bit_score\n')

    archive_files = glob.glob(os.path.join(args.root_dir, '*.tar.zst'))
    logging.info(f"Found {len(archive_files)} archives to process in '{args.root_dir}'")
    
    for tar_file_path in archive_files:
        with open(args.output_file, 'a') as writer:
            process_archive(tar_file_path, writer, args)

    logging.info(f"Extraction complete. All data saved to {args.output_file}")

if __name__ == "__main__":
    main()