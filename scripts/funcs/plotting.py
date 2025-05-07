import seaborn as sns
import matplotlib.pyplot as plt

def create_box_swarm_plot(df, x_name, y_name, axis, order=['Wildtype', 'S527del', 'Y528C'], palette='pastel', swarm_color='dimgrey', marker_size=5, title=None, y_line=None, y_lim=None, marker_column=None):
    """
    Create a box plot with a swarm plot overlay on a specified axis.

    Parameters:
    df (DataFrame): The input data.
    x_name (str): Column name for the x-axis (categorical variable).
    y_name (str): Column name for the y-axis (numerical variable).
    axis (matplotlib.axes.Axes): The axis to draw the plot on.
    order (list): Order of categories for the x-axis.
    palette (str): Color palette for the box plot.
    swarm_color (str): Color for the swarm plot points.
    marker_size (int): Size of the swarm plot markers.
    title (str): Title of the plot (optional).
    y_line (float): Horizontal line to add to the plot (optional).
    y_lim (tuple): Limits for the y-axis (optional).
    marker_column (str): Column name to assign different markers based on unique values (optional).

    Returns:
    None
    """

    # Create box plot
    sns.boxplot(
        data=df, x=x_name, y=y_name, order=order,
        palette=palette, showfliers=False, ax=axis,
        linewidth=0.5
    )

    # Overlay swarm plot
    if marker_column and marker_column in df.columns:
        # Assign different markers for unique values in the marker_column
        unique_markers = ['o', 'D', 's', '^', 'v', 'P', '*', 'X']  # Add more markers if needed
        marker_map = {value: unique_markers[i % len(unique_markers)] for i, value in enumerate(df[marker_column].unique())}

        for value, marker in marker_map.items():
            subset = df[df[marker_column] == value]
            sns.swarmplot(
                data=subset, x=x_name, y=y_name, order=order,
                color=swarm_color, alpha=0.75, ax=axis, size=marker_size, marker=marker,
                label=value  # Add legend labels for markers
            )
    else:
        sns.swarmplot(
            data=df, x=x_name, y=y_name, order=order,
            color=swarm_color, alpha=0.75, ax=axis, size=marker_size
        )

    # Add horizontal line if specified
    if y_line is not None:
        axis.axhline(y=y_line, linestyle='dotted', color='silver', zorder=0)

    # Set y-axis limits if specified
    if y_lim is not None:
        axis.set_ylim(y_lim)

    # Set title if specified
    if title:
        axis.set_title(title)

    return
