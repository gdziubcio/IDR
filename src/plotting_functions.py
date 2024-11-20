
import numpy as np 
import pandas as pd 
from matplotlib import pyplot as plt


def ps_plot(axes,
            data,
            thresh: float = 0.25,
            probabilities: None|dict = None,
            legends_pos: int = 2
        ) -> None:
    
    def find_true_regions(bool_array):
        """
        Finds the start and end indices of contiguous regions where bool_array is True.
        """
        bool_array = np.asarray(bool_array, dtype=bool)
        edges = np.diff(np.concatenate(([0], bool_array.view(np.int8), [0])))
        starts = np.where(edges == 1)[0]
        ends = np.where(edges == -1)[0]
        return list(zip(starts, ends - 1))
    
    if probabilities is None: 
        probabilities = {
            'B7WN96': 0.5495,
            'P22980': 0.6207,
            'Q10655': 0.77396,
            'Q17381': 0.82436,
            'Q18612': 0.40498,
            'P28515': 0.7854
        }
    
    axes = axes.flatten()

    for ax, (i, row) in zip(axes, data.iterrows()):
        uniprot_id = row['ID']
        name = row['gene_name']
        prob = probabilities.get(uniprot_id, None)  # Get the probability for this protein

        x = np.arange(len(row['fldpnn2_score']))
        y = np.array(row['fldpnn2_score'])
        ax.plot(x, y, color='#57cc99', linewidth=2)
        ax.axhline(thresh, color='red', linestyle='--', linewidth=1.5, label="Disorder Threshold")
        ax.fill_between(x, y, thresh, where=(y > thresh), color='grey', alpha=0.5, label="Potential IDRs")
    
        # Update the title to include the probability
        if probabilities is not None:
            title_text = f"{name}\nPhase Separation Probability: {prob:.2f}"
        else:
            title_text = name
        ax.set_title(title_text, fontsize=15)
    
        ax.set_yticks([0, thresh, 1])
        ax.set_ylim(0, 1)
        x_ticks = np.arange(0, len(x) + 1, 100)
        ax.set_xticks(x_ticks)
        if i >= 3:
            ax.set_xlabel('Position', fontsize=12)
        if i == 0 or i == 3:
            ax.set_ylabel('Disorder Propensity', fontsize=12)

        # Add vertical colored regions based on the boolean list
        bool_array = np.array(row['DRegion'], dtype=bool)  # Replace 'bool_list' with your column name
        regions = find_true_regions(bool_array)
        for start, end in regions:
            ax.axvspan(start, end, color='#147c89', alpha=0.3, label="Phase Sep. Key Residues")
    
        if i == legends_pos:
            ax.legend()
        
        ax.spines[['top', 'right']].set_visible(False)

    plt.tight_layout()
    plt.show()

###--------------------------

def residue_plot(axes,
            data,
            thresh: float = 0.25,
            legends_pos: None|int = None
        ) -> None:

    # Define properties of amino acids
    AROMATIC = {'F', 'Y', 'W'}
    CHARGE = {
        'D': -1, 'E': -1,
        'K': +1, 'R': +1,
        'H': +0.5  # Histidine has partial positive charge at physiological pH
    }
    
    axes = axes.flatten()
    for ax, (i, row) in zip(axes, data.iterrows()):
        name = row['gene_name']
        x = np.arange(len(row['fldpnn2_score']))
        y = np.array(row['fldpnn2_score'])
        ax.plot(x, y, color='#57cc99', linewidth=2, label='Disorder Propensity')
        ax.axhline(thresh, color='red', linestyle='--', linewidth=1.5, label='Disorder Threshold')
        ax.fill_between(x, y, thresh, where=(y > thresh), color='grey', alpha=0.5, label='Potential IDRs')
        ax.set_title(name, fontsize=20)
        ax.set_yticks([0, thresh, 1])
        ax.set_ylim(0, 1)
        x_ticks = np.arange(0, len(x) + 1, 100)
        ax.set_xticks(x_ticks)
        if i >= 3:
            ax.set_xlabel('Position', fontsize=19)
        if i == 0 or i == 3:
            ax.set_ylabel('Disorder Propensity', fontsize=19)

        # Plot aromatic residues as dots at y=0.4
        sequence = row['AA']
        aromatic_positions = [idx for idx, aa in enumerate(sequence) if aa in AROMATIC]
        ax.scatter(aromatic_positions, [0.4]*len(aromatic_positions), color='purple', marker='o', s=20, label='Aromatic Residues', zorder=5)

        # Calculate net charge per block of 5 residues
        charges = [CHARGE.get(aa, 0) for aa in sequence]
        window_size = 10
        net_charges = []
        net_charge_positions = []
        for idx in range(len(charges) - window_size + 1):
            block_charge = sum(charges[idx:idx + window_size])
            net_charges.append(block_charge)
            # Position the net charge at the center of the block
            net_charge_positions.append(idx + window_size // 2)

        # Scale down the net charges
        scaled_net_charges = [charge * 0.02 for charge in net_charges]  # Adjust scaling factor as needed
        # Shift the scaled net charges to fit within y-axis limits
        shifted_net_charges = [sc + 0.8 for sc in scaled_net_charges]  # Shift to center around 0.5

        # Plot net charge per block on the same y-axis
        ax.plot(net_charge_positions, shifted_net_charges, color='orange', linewidth=1.5, linestyle='-', label='Net Charge per 10-Residue Block', alpha=0.7)
        ax.plot(net_charge_positions, np.ones(len(shifted_net_charges)) * .8, color='black', linewidth=1.5, linestyle='--', alpha=0.3)
        # Adjust y-limits if necessary
        ax.set_ylim(0, 1)

        # Adjust legend to prevent duplicate entries
        if legends_pos is not None: 
            if i == legends_pos:
                handles, labels = ax.get_legend_handles_labels()
                by_label = dict(zip(labels, handles))
                ax.legend(by_label.values(), by_label.keys(), loc='upper right')
        
        ax.tick_params(axis='both', which='major', labelsize=13)
        ax.spines[['top', 'right']].set_visible(False)

    plt.tight_layout()
    plt.show()