# source: https://github.com/MichaelPudjihartono/regionperm.git

import pandas as pd
import numpy as np
import pybedtools
import matplotlib.pyplot as plt
import seaborn as sns

from tqdm import tqdm
import time


def add_fragment_length_column(df):
    df['fragment_length'] = df.apply(lambda row: row.iloc[2]-row.iloc[1], axis=1)
    return df

def add_feature_id_column(df):
    df['feature_id'] = df.apply(lambda row: str(row.iloc[0])+':'+str(row.iloc[1])+'-'+str(row.iloc[2]), axis=1)
    return df

def get_universe_intersect(universe_region_set, B_region_set):
    universe_region_set_BEDtool = pybedtools.BedTool.from_dataframe(universe_region_set)
    B_region_set_BEDtool = pybedtools.BedTool.from_dataframe(B_region_set)

    universe_intersect = universe_region_set_BEDtool.intersect(B_region_set_BEDtool, wa=True, wb=True)
    universe_intersect = universe_intersect.to_dataframe()
    universe_test_stat = len(set(universe_intersect.iloc[:,-1])) #The last column correspond to 'feature_id'
    #Thus, test statistic counts the total number of unique B region set that intersects the universe region set

    return universe_intersect, universe_test_stat

def get_intersect(region_set, universe_intersect):
    left_on = list(region_set.columns[0:3])
    right_on = list(universe_intersect.columns[0:3])

    intersect = pd.merge(region_set, universe_intersect, how='inner', left_on=left_on, right_on=right_on)
    test_stat = len(set(intersect.iloc[:,-1])) #The last column correspond to 'feature_id'
    #Thus, test statistic counts the total number of unique B region set that intersects the specified region set

    return test_stat

def create_complement_region_set(A_region_set, universe_region_set): #This creates the complement of A_region_set
    left_on = list(A_region_set.columns[0:3])
    right_on = list(universe_region_set.columns[0:3])

    complement_region_set = pd.merge(A_region_set, universe_region_set, how='outer', left_on=left_on, right_on=right_on)
    complement_region_set = complement_region_set[complement_region_set.isnull().any(axis=1)]

    return complement_region_set

def permutation_simulation(A_region_set, universe_region_set, universe_intersect, 
                           original_test_stat, num_iterations, match_by):
    
    perm_test_stats = [original_test_stat]
    more_than_original_count = 1
    less_than_original_count = 1

    bar_format = '{desc}: {percentage:3.0f}% |{bar}| {n_fmt}/{total_fmt} {unit}'

    for iteration in tqdm(
        range(num_iterations), 
        desc='Permutation simulations in progress', 
        unit='simulation', 
        bar_format = bar_format, 
        ncols=100):

        iteration_region_set = universe_region_set.sample(frac=1).reset_index().drop('index', axis=1)
        iteration_region_set['fragment_length_cumsum'] = iteration_region_set['fragment_length'].cumsum()

        if match_by == 'count':
            total_count = A_region_set.shape[0]
            up_to_index = total_count
        elif match_by == 'length':
            total_length = A_region_set['fragment_length'].sum()
            up_to_index = iteration_region_set[iteration_region_set['fragment_length_cumsum'] >= total_length].index[0] + 1
        
        iteration_region_set = iteration_region_set.iloc[:up_to_index]
        iteration_test_stat = get_intersect(iteration_region_set, universe_intersect)

        perm_test_stats.append(iteration_test_stat)

        if iteration_test_stat >= original_test_stat:
            more_than_original_count += 1
        elif iteration_test_stat <= original_test_stat:
            less_than_original_count += 1

    return perm_test_stats, more_than_original_count, less_than_original_count

def generate_perm_test_result(perm_test_stats, original_test_stat, complement_test_stat, 
                              universe_test_stat, more_than_original_count, less_than_original_count, num_iterations):
    
    alternative = 'greater' if more_than_original_count<less_than_original_count else 'less'

    if alternative == 'greater':
        p_val = more_than_original_count/len(perm_test_stats)
    elif alternative == 'less':
        p_val = less_than_original_count/len(perm_test_stats)
    
    #Calculate z-score
    mean = sum(perm_test_stats) / len(perm_test_stats)
    std_dev = (sum((x - mean) ** 2 for x in perm_test_stats) / len(perm_test_stats)) ** 0.5
    z_score = (original_test_stat - mean) / std_dev

    #Generate perm_test_result df
    perm_test_result = pd.DataFrame({
                                     'p_value':[p_val],
                                     'z_score':[z_score],
                                     'n_iterations':[num_iterations],
                                     'alternative':[alternative],
                                     'original_evaluation':[original_test_stat],
                                     'complement_evaluation':[complement_test_stat],
                                     'universe_evaluation':[universe_test_stat]})
    
    return perm_test_result

def generate_perm_test_figure(perm_test_stats, perm_test_result):
    fig, ax = plt.subplots(figsize=(8, 6))

    #Plot histogram
    sns.histplot(data=perm_test_stats, bins=100, kde=False, color='#a5aab0', stat='density', alpha=0.5, ax=ax)

    #Set axis labels and texts
    ax.set_xlabel('Number of overlaps', weight='bold')
    ax.set_ylabel('Density', weight='bold')

    text = f"p-value: {perm_test_result.loc[0,'p_value']}\nz-score: {perm_test_result.loc[0,'z_score']}\nn-iterations: {perm_test_result.loc[0, 'n_iterations']}"
    ax.text(x=0.5, y=1.05, s=text, transform=ax.transAxes, ha='center')

    #Add vertical line to indicate significance_cutoff and observed value
    percentile = 95 if perm_test_result.loc[0, 'alternative'] == 'greater' else 5
    significance_cutoff = np.percentile(perm_test_stats, percentile)
    ax.axvline(x=significance_cutoff, ymax=0.7, color= 'red', linestyle='--')
    ax.text(x=significance_cutoff, y=ax.get_ylim()[1]*0.72, s=f'{int(significance_cutoff)}',
           color='red', va='center', ha='center', fontsize=8, weight='bold')
    
    ax.axvline(x=perm_test_result.loc[0, 'original_evaluation'], ymax=0.7, color='green', linestyle='--')
    ax.text(x=perm_test_result.loc[0, 'original_evaluation'], y=ax.get_ylim()[1]*0.72, s=f"{int(perm_test_result.loc[0, 'original_evaluation'])}",
           color='green', va='center', ha='center', fontsize=8, weight='bold')
    
    ax.legend(['\u03B1 = 0.05', 'Observed'])

    #Plot KDE line
    sns.kdeplot(data=perm_test_stats, color='black', ax=ax)
    ax.axvline(x=significance_cutoff, ymax=0.7, color= 'red', linestyle='--')
    ax.axvline(x=perm_test_result.loc[0, 'original_evaluation'], ymax=0.7, color='green', linestyle='--')

    plt.tight_layout()



def run(A_region_set, B_region_set, universe_region_set, num_iterations, match_by):
    '''Genomic region association analysis with permutation tests.

    Assess the association between a set of genomic regions and other genomic features using permutation tests.

    Parameters:
    A_region_set: the set of regions to randomize in BED format. 
    B_region_set: the set of genomic features that will be analyzed for their association with the regions in A_region_set in BED format.
    universe_region_set: the total valid regions from which subsequent random iterations will be sampled from in BED format. A_region_set is a subset of the universe_region_set.
    num_iterations: number of iterations for permutation test.
    match_by: for each iteration, randomize regions by matching the count or length of the original A region set. choices=['count', 'length']
    '''

    start_time = time.time()
    
    A_region_set = add_fragment_length_column(A_region_set)
    B_region_set = add_feature_id_column(B_region_set)

    universe_region_set = add_fragment_length_column(universe_region_set)

    universe_intersect, universe_test_stat = get_universe_intersect(universe_region_set, B_region_set)

    original_test_stat = get_intersect(A_region_set, universe_intersect)

    complement_region_set = create_complement_region_set(A_region_set, universe_region_set)
    complement_test_stat = get_intersect(complement_region_set, universe_intersect)

    perm_test_stats, more_than_original_count, less_than_original_count = permutation_simulation(A_region_set, universe_region_set, universe_intersect, 
                                                                                                 original_test_stat, num_iterations, match_by)
    
    perm_test_result = generate_perm_test_result(perm_test_stats, original_test_stat, complement_test_stat, 
                                                 universe_test_stat, more_than_original_count, less_than_original_count, num_iterations)
    
    generate_perm_test_figure(perm_test_stats, perm_test_result)
    
    print(f'Original evaluation: {original_test_stat}')
    print(f"Alternative: {perm_test_result.loc[0,'alternative']}")
    print(f"P-value: {perm_test_result.loc[0,'p_value']}")
    print(f"Z-score: {perm_test_result.loc[0,'z_score']}\n")

    print('Done.')
    print(f'Time elapsed: {(time.time() - start_time) / 60:.2f} minutes.')

    return perm_test_result
