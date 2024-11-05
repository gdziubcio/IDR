import pandas as pd

def get_scores(path: str = "data/fldpnn2.txt") -> pd.DataFrame:
    """
    Reads a txt file from fldpnn2 containing disorder sequences, 
    processes the data, and returns a DataFrame.

    Args:
        path (str, optional): The file path to the disorder sequences data.
                              Defaults to "data/fldpnn2.txt".

    Returns:
        pd.DataFrame: A DataFrame containing two columns:
                      - "ID": The cleaned sequence identifiers.
                      - "fldpnn2_score": The corresponding disorder scores.
    """
    # Initialize an empty dictionary to store sequence IDs and their scores
    fldpnn2_dict: dict = {}

    with open(path, mode="r") as file:

        lines: list = file.readlines()
        
        # header occupies always first 8 lines
        cleaned_lines = [line.strip() for line in lines[8:]]
        
        # Extract sequence IDs and disorder scores:
        # - Sequence IDs are assumed to be every 5th line starting from the first
        # - Disorder scores are assumed to be every 5th line starting from the fifth
        sequence_id = cleaned_lines[0::5]
        idr_ranges = cleaned_lines[1::5]
        aa_sequence = cleaned_lines[2::5]
        is_disordered = cleaned_lines[3::5]
        disorder_scores = cleaned_lines[4::5]
        
        # Create a dictionary mapping each sequence ID to its disorder score
        fldpnn2_dict = {
            name: disorder_seq 
            for name, disorder_seq in zip(sequence_id, disorder_scores)
        }


    df_fldpnn2 = pd.DataFrame(
        fldpnn2_dict.items(), 
        columns=["ID", "fldpnn2_score"]
    )
    
    # clean the sequence identifiers
    df_fldpnn2["ID"] = df_fldpnn2["ID"].str.replace(">", "", regex=False)
    
    return df_fldpnn2
