#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  1 15:14:02 2025

@author: matthew
"""

#%% select_points()



def select_points(x, y, *, marker='o', size=30, color='tab:blue'):
    """
    Plot (x, y) as a scatter.  The user draws a freehand polygon
    with the mouse; the function returns the indices of the points
    that fall inside that polygon.

    Returns
    -------
    idx : ndarray
        1-D array of integer indices into (x, y) that were selected.
    """
    
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.widgets import PolygonSelector
    from matplotlib.path import Path
    
    x, y = np.asarray(x), np.asarray(y)
    fig, ax = plt.subplots()
    sc = ax.scatter(x, y, c=color, s=size, marker=marker)
    ax.set_title("Draw polygon → press 'enter' to accept, 'escape' to cancel")

    selected_idx = []

    def onselect(verts):
        """Callback from PolygonSelector when the user finishes."""
        nonlocal selected_idx
        path = Path(verts)
        points = np.column_stack((x, y))
        selected = path.contains_points(points)
        selected_idx = np.where(selected)[0]

        # highlight selection
        sc.set_facecolors(['red' if s else color for s in selected])
        fig.canvas.draw_idle()

    selector = PolygonSelector(ax, onselect, useblit=True)

    def accept(event):
        if event.key == "enter":
            plt.close(fig)           # close and return indices
        elif event.key == "escape":
            nonlocal selected_idx
            selected_idx = np.array([], dtype=int)  # nothing
            plt.close(fig)

    fig.canvas.mpl_connect("key_press_event", accept)
    plt.show(block=True)

    return selected_idx


#%% timeline_and_gap_stack()

def timeline_and_gap_stack(ifg_list):
    """
    Figure with
      • top: horizontal timeline of unique dates
      • bottom: stacked dots showing *all* positive day-gaps between dates
    """
    import matplotlib.pyplot as plt
    from matplotlib.dates import AutoDateLocator, DateFormatter
    from datetime import datetime
    from itertools import combinations
    import numpy as np
    from collections import Counter
    
    # 1 ── unique dates -----------------------------------------------------
    unique = {d for pair in ifg_list for d in pair.split('_')}
    dates = sorted(datetime.strptime(d, "%Y%m%d") for d in unique)

    # 2 ── all positive gaps (days) ----------------------------------------
    ordinals = [d.toordinal() for d in dates]
    gaps = [later - earlier for earlier, later in combinations(ordinals, 2)]

    # 3 ── prepare stacked-dot coordinates ---------------------------------
    counts = Counter(gaps)                  # how many times each gap appears
    xs, ys = [], []
    for gap, freq in sorted(counts.items()):
        xs.extend([gap] * freq)             # same x for duplicates
        ys.extend(range(freq))              # 0,1,2,… stack upward

    # 4 ── build figure -----------------------------------------------------
    fig, (ax_tl, ax_stack) = plt.subplots(
        2, 1, figsize=(11, 4.5),
        gridspec_kw=dict(height_ratios=[1, 2]),
        sharex=False
    )

    # ── top: timeline ------------------------------------------------------
    ax_tl.plot(
        dates, np.zeros_like(dates),
        marker='|', ms=16, lw=0, color='tab:blue'
    )
    ax_tl.set_yticks([])
    ax_tl.set_title("Unique interferogram acquisition dates")
    ax_tl.margins(y=0.25)
    ax_tl.xaxis.set_major_locator(AutoDateLocator())
    ax_tl.xaxis.set_major_formatter(DateFormatter('%Y-%m'))

    # ── bottom: stacked dots of day-gaps -----------------------------------
    ax_stack.scatter(xs, ys, s=18, color='tab:green')
    ax_stack.set_xlabel("Days between any two acquisitions")
    ax_stack.set_ylabel("Count (stacked)")
    ax_stack.set_title("Stacked distribution of all positive date gaps")
    ax_stack.grid(axis='x', alpha=0.2)

    fig.tight_layout()
    plt.show()



#%%


def maps_tcs_rescale(maps, tcs):
    """
    A function to rescale spaital maps to have unit range and rescale each's time cource (tc)
    so that there is no change to the product of the two matrices.  
    
    
    input:
        maps | array | spatial maps as rows (e.g. 2x1775)
        tcs  | array | time courses as columns (e.g. 15x2)
    Output:
        maps_scaled | array | spatial maps as rows with each row having unit range
        tcs_scaled | array | TCs scaled so that new maps x new tcs equals maps x tcs
    
    2017/05/15 | written
    2024_02_29 | MEG | Checked that no change to product of two matrices.  
    
    """
    import numpy as np
    
    def rescale_unit_range(signals):
        """
        rescale a matrix of row vectors so that each row has a range of 1.  
        Also record the scaling factor required to rescale each row vector
        signals are rows
        Input:
            signals: array with each signal as a new row 
        Output:
            signals_rescale | array | each row has a range of 1
            signals_factor | array | factor that each row is dvided by to adjust range
        """
        import numpy as np
        
        # initiate outputs
        signals_rescale = np.ones(signals.shape)                                    
        signals_factor = np.ones((np.size(signals, axis=0) , 1))                    
           
        # iterate over each signal
        for i in np.arange(np.size(signals, axis = 0)):
            # get the range of the signal
            signals_factor[i,0] = (np.max(signals[i,:])-np.min(signals[i,:]))
            # divide by range, so it's now 1 for the signal.  
            signals_rescale[i,:] = signals[i,:]/signals_factor[i,0]
            
        return signals_rescale, signals_factor

    
    maps_scaled , scaling = rescale_unit_range(maps)
    tcs_scaled = tcs * np.ravel(scaling)
    return maps_scaled, tcs_scaled

