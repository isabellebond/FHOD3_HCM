import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from funcs.stats import anova_with_tukey
from funcs.plotting import create_box_swarm_plot

def format_data(df):
    """
    Format the data for analysis.
    
    Parameters:
    karytyping_data (DataFrame): The input data containing karyotyping results.
    
    Returns:
    DataFrame: The formatted data ready for analysis.
    """
    # Define genotype order
    df['Genotype'] = pd.Categorical(
        df['Genotype'], 
        ordered=True, 
        categories=['Control', 'Wildtype', 'S527del', 'Y528C']
    ) 
    # Convert 'Copy number' to numeric
    df['Copy number'] = pd.to_numeric(df['Copy number'], errors='coerce')

    #Multiple copy number by 2 as KOLF is XY genotype
    df.loc[df['Chromosome'] == 'X', 'Copy number'] = df.loc[df['Chromosome'] == 'X', 'Copy number']*2

    df.dropna(subset=['Copy number'], inplace=True)

    return df

def main():
    # Load the data
    file_path = 'data/karyotyping/karyotyping_results.csv' 
    data = pd.read_csv(file_path)
    data = format_data(data)
    # Print the first few rows of the data
    print(data.head())

    ######################################################################
    ################# Statistical Analysis of Karyotyping ################
    ######################################################################
    # Analysis suggested as per https://cdn.stemcell.com/media/files/pis/10000000353-PIS_01.pdf

    genotypes = ['Wildtype','S527del', 'Y528C']
    experiments = [1,2]
    anova_results = []

    for experiment in experiments:
        for chromosome in data['Chromosome'].unique():
            chromosome_data = data.loc[((data['Experiment'] == experiment)) & (data['Chromosome'] == chromosome)]
            if len(chromosome_data) > 0:
                results = anova_with_tukey(chromosome_data, 'Genotype', 'Copy number')
                results['Chromosome'] = chromosome
                results['Experiment'] = experiment
                anova_results.append(results)
                
    results = pd.concat(anova_results)
    print(results.columns)
    results = results.loc[results['group1'] == 'Control']
    results = results[['Chromosome', 'Experiment', 'group2', 'ANOVA F', 'ANOVA P', 'ANOVA Reject', 'Tukey Reject']]
    

    results.to_csv('results/statistical_analyses/karyotyping.csv')

    ######################################################################
    ################# Karyotyping Plots ##################################
    ######################################################################
    chromosome_order = [str(i) for i in range(1, 23)] + ['X', 'Y']
    data['Chromosome'] = pd.Categorical(data['Chromosome'].astype(str), categories=chromosome_order, ordered=True)
    data = data.loc[data['Experiment'] == 1]

    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(7.09, 7.09))
    axes = axes.flatten()
    i=0
    # Iterate through chromosomes and plot
    for chromosome in chromosome_order:
        if chromosome in data['Chromosome'].astype(str).unique():
            chrom_data = data[data['Chromosome'] == chromosome]
            ax = axes[i]
            create_box_swarm_plot(
                    chrom_data, 'Genotype', 'Copy number', axis=ax,
                    order=['Control', 'Wildtype', 'S527del', 'Y528C'],
                    y_lim=(0.5, 3.5), y_line= 2, swarm_color='dimgrey',
                    title=f'Chromosome {chromosome}-{chrom_data["Arm"].iloc[0]}',
                    marker_size=4, marker_column='Day'
                )
            
            if chromosome == '4':
                ax.set_title(f'Chromosome {chromosome}-{chrom_data["Arm"].iloc[0]} (Control)')

            ax.tick_params(axis='x', rotation=45, labelrotation=45, labelright=False)
            ax.set_xlabel('')
            ax.set_ylabel('Copy Number')
            i += 1

    plt.tight_layout()
    plt.savefig('results/figures/karyotyping.png', dpi=300)

if __name__ == "__main__":
    main()