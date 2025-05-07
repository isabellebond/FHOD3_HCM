import pandas as pd
from scipy.stats import f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd



def anova_with_tukey(df, group_column, data_column):
    """
    Perform one-way ANOVA and Tukey's post hoc test.

    Parameters:
    df (DataFrame): The input data containing groups and values to compare.
    group_column (str): Column name for the grouping variable (e.g., 'Genotype').
    data_column (str): Column name for the numerical variable to compare (e.g., 'Copy number').
    p_threshold (float): P-value threshold for significance (default is 0.05).

    Returns:
    DataFrame: A DataFrame containing the ANOVA F-statistic, P-value, and Tukey's HSD results (if significant).
    """


    # Perform one-way ANOVA
    groups = [group[data_column].dropna() for _, group in df.groupby(group_column)]
    stat, p = f_oneway(*groups)

    # Initialize results DataFrame
    results = pd.DataFrame({
        'F-statistic': [stat],
        'P-value': [p]
    })

    # If ANOVA is significant, perform Tukey's post hoc test
    tukey = pairwise_tukeyhsd(df[data_column], df[group_column])
    tukey_df = pd.DataFrame(data=tukey.summary().data[1:], columns=tukey.summary().data[0])

    tukey_df.rename(columns={'reject': 'Tukey Reject'}, inplace=True)
    tukey_df['ANOVA P'] = p
    tukey_df['ANOVA F'] = stat
    tukey_df['ANOVA Reject'] = tukey_df['ANOVA P'].apply(lambda x: True if x < 0.05 else False)

    return tukey_df

def kruskal_with_dunn(control_data, outcome_data, comparison_column, p_threshold = 0.05):
    """
    Perform Kruskal-Wallis test and Dunn's post hoc test.

    Parameters:
    control_data (DataFrame): Data for the control group.
    outcome_data (DataFrame): Data for the outcome group.
    comparison_column (str): Column name to compare (e.g., 'Copy number').

    Returns:
    dict: A dictionary containing the Kruskal-Wallis H-statistic, P-value, and Dunn's test results (if significant).
    """
    # Perform Kruskal-Wallis test
    stat, p = kruskal(control_data[comparison_column], outcome_data[comparison_column])
    
    # Initialize results DataFrame
    results = pd.DataFrame({
        'H-statistic': [stat],
        'P-value': [p]
    })

    print(results)

    # If Kruskal-Wallis is significant, perform Dunn's post hoc test
    if p < p_threshold:
        combined_data = pd.concat([control_data, outcome_data])
        dunn = sp.posthoc_dunn(combined_data[comparison_column], combined_data['Genotype'])
        dunn_df = pd.DataFrame(data=dunn.summary().data[1:], columns=dunn.summary().data[0])
        results = pd.concat([results, dunn_df], axis=1)

    return results