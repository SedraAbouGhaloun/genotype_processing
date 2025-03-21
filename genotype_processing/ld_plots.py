import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def plot_ld_rank_scatter(ld, chromosome, r2_threshold=0.05, window_size=100):
    ld_filtered = ld[
        (ld['CHR_A'] == chromosome) & 
        (ld['CHR_B'] == chromosome) & 
        (ld['R2'] >= r2_threshold)
    ]

    if ld_filtered.empty:
        print(f"No LD pairs found for chromosome {chromosome} with R² >= {r2_threshold}.")
        return

    unique_positions = pd.unique(ld_filtered[['BP_A', 'BP_B']].values.ravel())
    unique_positions.sort()
    position_to_rank = {pos: rank for rank, pos in enumerate(unique_positions, start=1)}

    ld_filtered['Rank_POS1'] = ld_filtered['BP_A'].map(position_to_rank)
    ld_filtered['Rank_POS2'] = ld_filtered['BP_B'].map(position_to_rank)

    max_dist_row = ld_filtered.loc[ld_filtered['dist'].idxmax()]
    center_rank1 = position_to_rank[max_dist_row['BP_A']]
    center_rank2 = position_to_rank[max_dist_row['BP_B']]

    lower_bound = max(min(center_rank1, center_rank2) - window_size // 2, 1)
    upper_bound = min(max(center_rank1, center_rank2) + window_size // 2, len(unique_positions))

    all_ranks = list(range(lower_bound, upper_bound + 1))

    ld_window = ld_filtered[
        (ld_filtered['Rank_POS1'] >= lower_bound) & 
        (ld_filtered['Rank_POS1'] <= upper_bound) & 
        (ld_filtered['Rank_POS2'] >= lower_bound) & 
        (ld_filtered['Rank_POS2'] <= upper_bound) &
        (ld_filtered['Rank_POS1'] <= ld_filtered['Rank_POS2'])
    ]

    self_correlation = pd.DataFrame({
        'Rank_POS1': all_ranks,
        'Rank_POS2': all_ranks,
        'R2': [1.0] * len(all_ranks)
    })

    ld_window = pd.concat([ld_window, self_correlation], ignore_index=True)

    all_pairs = pd.DataFrame([(x, y) for x in all_ranks for y in all_ranks if x <= y], columns=['Rank_POS1', 'Rank_POS2'])

    plt.figure(figsize=(10, 8))
    plt.scatter(all_pairs['Rank_POS1'], all_pairs['Rank_POS2'], color='grey', s=1, alpha=0.5)
    plt.scatter(ld_window['Rank_POS1'], ld_window['Rank_POS2'], c=ld_window['R2'], cmap='Reds', s=20, alpha=1, marker='s')

    plt.colorbar(label='R²')
    plt.title(f"LD Scatter Plot for Chromosome {chromosome} (Rank-Based Window)")
    plt.xlabel("Rank Position (POS1)")
    plt.ylabel("Rank Position (POS2)")
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.show()

def plot_ld_rank_scatter(ld, chromosome, r2_threshold=0.05, window_size=100):
    ld_filtered = ld[
        (ld['CHR_A'] == chromosome) & 
        (ld['CHR_B'] == chromosome) & 
        (ld['R2'] >= r2_threshold)
    ]

    if ld_filtered.empty:
        print(f"No LD pairs found for chromosome {chromosome} with R² >= {r2_threshold}.")
        return

    unique_positions = pd.unique(ld_filtered[['BP_A', 'BP_B']].values.ravel())
    unique_positions.sort()
    position_to_rank = {pos: rank for rank, pos in enumerate(unique_positions, start=1)}

    ld_filtered['Rank_POS1'] = ld_filtered['BP_A'].map(position_to_rank)
    ld_filtered['Rank_POS2'] = ld_filtered['BP_B'].map(position_to_rank)

    max_dist_row = ld_filtered.loc[ld_filtered['dist'].idxmax()]
    center_rank1 = position_to_rank[max_dist_row['BP_A']]
    center_rank2 = position_to_rank[max_dist_row['BP_B']]

    lower_bound = max(min(center_rank1, center_rank2) - window_size // 2, 1)
    upper_bound = min(max(center_rank1, center_rank2) + window_size // 2, len(unique_positions))

    all_ranks = list(range(lower_bound, upper_bound + 1))

    ld_window = ld_filtered[
        (ld_filtered['Rank_POS1'] >= lower_bound) & 
        (ld_filtered['Rank_POS1'] <= upper_bound) & 
        (ld_filtered['Rank_POS2'] >= lower_bound) & 
        (ld_filtered['Rank_POS2'] <= upper_bound) &
        (ld_filtered['Rank_POS1'] <= ld_filtered['Rank_POS2'])
    ]

    self_correlation = pd.DataFrame({
        'Rank_POS1': all_ranks,
        'Rank_POS2': all_ranks,
        'R2': [1.0] * len(all_ranks)
    })

    ld_window = pd.concat([ld_window, self_correlation], ignore_index=True)

    all_pairs = pd.DataFrame([(x, y) for x in all_ranks for y in all_ranks if x <= y], columns=['Rank_POS1', 'Rank_POS2'])

    plt.figure(figsize=(10, 8))
    plt.scatter(all_pairs['Rank_POS1'], all_pairs['Rank_POS2'], color='grey', s=1, alpha=0.5)
    plt.scatter(ld_window['Rank_POS1'], ld_window['Rank_POS2'], c=ld_window['R2'], cmap='Reds', s=20, alpha=1, marker='s')

    plt.colorbar(label='R²')
    plt.title(f"LD Scatter Plot for Chromosome {chromosome} (Rank-Based Window)")
    plt.xlabel("Rank Position (POS1)")
    plt.ylabel("Rank Position (POS2)")
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.show()





