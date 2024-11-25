import pandas as pd
import re

def merge_ranges(text, gap=10):
    # Extract all ranges as tuples of integers
    ranges = [tuple(map(int, s.split('-'))) for s in re.findall(r'\d+-\d+', text)]
    # Merge ranges based on the gap
    merged = []
    for start, end in ranges:
        if merged and start - merged[-1][1] <= gap:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def get_scores(path: str = "data/fldpnn2/fldpnn2_master.txt") -> pd.DataFrame:
    """
    Reads a txt file from fldpnn2 containing disorder sequences, 
    processes the data, and returns a DataFrame.

    Args:
        path (str, optional): The file path to the disorder sequences data.
                              Defaults to "data/fldpnn2.txt".

    Returns:
        pd.DataFrame: A DataFrame containing two columns:
                      - "ID": The cleaned sequence identifiers.
                      - "IDR": The corresponding disordered regions.
                      - "AA": The corresponding amino acid sequences.
                      - "Disordered": The corresponding disordered regions.
                      - "fldpnn2_score": The corresponding disorder scores.
    """

    with open(path, mode="r") as file:

        lines: list = file.readlines()
        
        # header occupies always first 8 lines
        cleaned_lines = [line.strip() for line in lines[8:]]
        
        # Extract sequence information from apropiate lines

        sequence_id = cleaned_lines[0::5]
        idr_ranges = cleaned_lines[1::5]
        aa_sequence = cleaned_lines[2::5]
        is_disordered = cleaned_lines[3::5]
        disorder_scores = cleaned_lines[4::5]
    
    # create a DataFrame from the extracted data

    df_fldpnn2 = pd.DataFrame(
        {
            "ID": sequence_id,
            "idr_ranges": idr_ranges,
            "AA": aa_sequence,
            "Disordered": is_disordered,
            "fldpnn2_score": disorder_scores,
        }
    )
    
    # convert the sequence string data to apropriate type 
    df_fldpnn2["ID"] = df_fldpnn2["ID"].str.replace(">", "", regex=False)
    df_fldpnn2['merged_ranges'] = df_fldpnn2['idr_ranges'].apply(merge_ranges)
    df_fldpnn2["AA"] = df_fldpnn2["AA"].apply(lambda x: [aa for aa in x.split(",")])
    df_fldpnn2["Disordered"] = df_fldpnn2["Disordered"].apply(lambda x: [int(d) for d in x.split(",")]) 
    df_fldpnn2["fldpnn2_score"] = df_fldpnn2["fldpnn2_score"].apply(lambda x: [float(s) for s in x.split(",")])   
    
    return df_fldpnn2