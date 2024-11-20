from pathlib import Path
import pandas as pd

def get_TF_fasta(id_path: str = "data/id_mapping.csv", 
                 output_path: str = "data/protein_sequences.fasta", 
                 name: str = "Entry"
    ) -> None:
    """
    Reads protein sequences from a CSV file and writes them to a FASTA file format.

    Args:
        id_path (str): Path to the CSV file containing protein IDs and sequences. 
                       Default is "data/id_mapping.csv".
        output_path (str): Path to save the output FASTA file. Default is "data/protein_sequences.fasta".
    """
    # Define the project root directory
    project_root = Path(__file__).resolve().parent.parent
    abs_id_path = project_root / id_path
    abs_output_path = project_root / output_path

    # Read data and process
    df_id = pd.read_csv(abs_id_path)
    protein_sequences = df_id[[name, "Sequence"]]
    get_fasta = lambda x: ">" + x.iloc[0] + "\n" + x.iloc[1] + "\n"

    with open(abs_output_path, "w") as f:
        f.write("\n".join(protein_sequences.apply(get_fasta, axis=1)))
