import re 
import pandas as pd

def parse_protein_data(filepath:str='data/psp.txt') -> pd.DataFrame:
    proteins = []
    protein_data = None

    with open(filepath, 'r') as file:
        for line in file:
            line = line.strip()
            
            if line.startswith('>'):
                # Save the previous protein data if it exists
                if protein_data:
                    proteins.append(protein_data)
                
                # Initialize a new protein data entry
                protein_data = {
                    "name": line[1:],  # Protein name
                    "AA": [],
                    "DRegion": []
                }
                
            elif re.match(r'^\d+\s+\w+\s+[-.\d]+\s+\d$', line):
                # Parse the data line into position, amino acid, probability, and region
                pos, aa, prob, dregion = line.split()
                
                # Append data to the current protein's lists
                protein_data["AA"].append(aa)
                protein_data["DRegion"].append(bool(int(dregion)))
                
    # Append the last protein data after loop ends
    if protein_data:
        proteins.append(protein_data)
    
    return pd.DataFrame(proteins)