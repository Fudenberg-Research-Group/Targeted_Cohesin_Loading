from scipy.ndimage import zoom
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

def make_equipotential_lines(heatmap_data, line_values, expand_cof=10, label=None, legend=None, color='red', **kwargs):
    """
    Plots equipotential contour lines based on a given heatmap.

    Parameters:
    -----------
    heatmap_data : np.ndarray
        A 2D numpy array representing the heatmap data.
    line_values : list
        A list of contour levels to be plotted.
    expand_cof : float, optional (default=1)
        Expansion coefficient to scale the heatmap before plotting.
    label : str, optional
        Label for the plotted lines (used in the legend).
    legend : str, optional
        Not currently used in the function.
    color : str, optional (default='red')
        Color of the contour lines.
    **kwargs : dict
        Additional keyword arguments passed to `plt.contour`.

    Returns:
    --------
    line : matplotlib.lines.Line2D or None
        A reference to the plotted line if `label` is provided, otherwise None.
    """

    expanded_d = zoom(heatmap_data, zoom=expand_cof, order=3)
    heatmap_shape_x = heatmap_data.shape[1]
    heatmap_shape_y = heatmap_data.shape[0]
    cntr1 = plt.contour(
        np.linspace(0, heatmap_shape_x, int(heatmap_shape_x * expand_cof)),
        np.linspace(0, heatmap_shape_y, int(heatmap_shape_y * expand_cof)),
        expanded_d, levels=line_values, colors=color, **kwargs
    )

    line = None
    if label is not None:
        line, = plt.plot([], [], color=color, linestyle='dashed', alpha=0.83, label=label)

    return line
