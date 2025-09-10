#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  1 15:14:02 2025

@author: matthew
"""


#%% interactive_mask_explorer()




def interactive_mask_explorer(stack, dates, *, cmap="viridis"):
    """
    Interactive viewer for a 3-D array (t, y, x) with a time slider.

    Parameters
    ----------
    stack : ndarray or np.ma.MaskedArray, shape (nt, ny, nx)
        Time series of images.
    dates : sequence of str / datetime / numpy.datetime64
        One label per time step (len == nt). Shown above the image.
    cmap : str, optional
        Matplotlib colormap name.
    vmin, vmax : float, optional
        Color limits. If None, set from robust percentiles across the stack.
    """
    from matplotlib.widgets import Slider
    import numpy as np
    import numpy.ma as ma
    import matplotlib.pyplot as plt
    
    if stack.ndim != 3:
        raise ValueError("stack must be 3-D (time, y, x)")
    nt, ny, nx = stack.shape
    if len(dates) != nt:
        raise ValueError(f"len(dates)={len(dates)} must equal nt={nt}")

    # Robust color limits (ignore masked/NaN)
    # if vmin is None or vmax is None:
    #     arr = stack.filled(np.nan) if ma.isMaskedArray(stack) else stack.astype(float)
    #     vmin = np.nanpercentile(arr, 1)  if vmin is None else vmin
    #     vmax = np.nanpercentile(arr, 99) if vmax is None else vmax
        
    vmin=np.nanmin(stack)
    vmax=np.nanmax(stack)

    # Helper to stringify dates
    def date_str(d):
        # datetime/date objects get strftime if available
        if hasattr(d, "strftime"):
            return d.strftime("%Y-%m-%d")
        # numpy datetime64 -> python datetime (safe-ish)
        try:
            return str(np.datetime_as_string(np.datetime64(d), unit='D'))
        except Exception:
            return str(d)

    # Start at the last frame
    idx0 = nt - 1
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    plt.subplots_adjust(bottom=0.18)  # room for slider

    im = ax.imshow(stack[idx0], cmap=cmap, vmin=vmin, vmax=vmax, origin="upper")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.set_xlabel("x (col)")
    ax.set_ylabel("y (row)")
    ax.set_title(f"{date_str(dates[idx0])}   (t={idx0}/{nt-1})")

    # Slider (0 … nt-1, integer steps)
    ax_s = fig.add_axes([0.12, 0.07, 0.76, 0.04])
    s = Slider(ax_s, "time", 0, nt-1, valinit=idx0, valstep=1, valfmt="%0.0f")

    def update(_):
        i = int(s.val)
        im.set_data(stack[i])
        ax.set_title(f"{date_str(dates[i])}   (t={i}/{nt-1})")
        fig.canvas.draw_idle()

    s.on_changed(update)

    # Nice-to-have: arrow keys to step frames
    def on_key(event):
        if event.key in ("left", "right"):
            i = int(s.val) + (-1 if event.key == "left" else 1)
            i = max(0, min(nt-1, i))
            if i != int(s.val):
                s.set_val(i)

    fig.canvas.mpl_connect("key_press_event", on_key)
    plt.tight_layout()
    return fig, ax, s




#%%
                
def interactive_tc_correlation_explorer(
    x,
    y,
    imgs,
    imgs_reco,
    A,
    mask,
    *,
    date_pairs,
    sources=None,
    thresh_frac: float = 0.01,
    scatter_kw: dict | None = None,
    cmap: str = "viridis",
    resid_cmap: str = "RdBu_r",
    highlight_kw: dict | None = None,
    ):
    """
    Interactive PCA/ICA explorer with an extra subplot that shows each sample's
    (date₁, date₂) position.

    Parameters
    ----------
    x, y           : 1-D arrays – coordinates for the primary scatter plot
    date_pairs     : list[str]  – strings "yyyymmdd_yyyymmdd", *one per row*
    imgs, imgs_reco
    A, mask, sources, thresh_frac, scatter_kw, cmap, resid_cmap
                    : same meaning as in the original function
    highlight_kw   : dict – style overrides for the highlight marker
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from matplotlib.widgets import Cursor
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    from datetime import datetime
    from licsalert.aux import col_to_ma
    
    
    
    # ── sanity checks ───────────────────────────────────────────────────────
    x = np.asarray(x)
    y = np.asarray(y)
    imgs       = np.asarray(imgs)
    imgs_reco  = np.asarray(imgs_reco)
    A          = np.asarray(A)
    date_pairs = np.asarray(date_pairs)

    if len(date_pairs) != x.size:
        raise ValueError("`date_pairs` must have the same length as `x`/`y`.")
    if imgs.shape != imgs_reco.shape:
        raise ValueError("imgs and imgs_reco must have the same shape")
    if not (imgs.shape[0] == x.size == y.size == A.shape[0]):
        raise ValueError("First dimension of x, y, imgs and A must match")
    if imgs.ndim != 2:
        raise ValueError("imgs must be a 2-D array of row vectors (N, ny*nx)")
    if mask.ndim != 2:
        raise ValueError("mask must be a 2-D array (ny, nx)")

    ny, nx      = mask.shape
    n_sources   = 0 if sources is None else sources.shape[0]
    highlight_kw = {"s": 60, "c": "crimson", "marker": "X", "zorder": 4,
                    **(highlight_kw or {})}

    # ── parse dates once ----------------------------------------------------
    def _split(p):
        a, b = p.split("_")
        return (datetime.strptime(a, "%Y%m%d"),
                datetime.strptime(b, "%Y%m%d"))
    date1, date2 = zip(*map(_split, date_pairs))
    date1, date2 = np.array(date1), np.array(date2)

    dmin, dmax = min(date1.min(), date2.min()), max(date1.max(), date2.max())

    # ── figure & layout -----------------------------------------------------
    fig = plt.figure(figsize=(15, 6))

    if n_sources:
        # Two-row grid: sources on top (full width) + bottom row with 3 panels
        gs_root = fig.add_gridspec(
            2, 1, height_ratios=[0.6, 2.4], hspace=0.15
        )
        gs_src   = gs_root[0].subgridspec(1, n_sources, wspace=0.05)
        gs_bottom = gs_root[1].subgridspec(
            1, 3, width_ratios=[1.3, 1.3, 2.6], wspace=0.30
        )
        # -- top sources -----------------------------------------------------
        for i in range(n_sources):
            ax = fig.add_subplot(gs_src[0, i])
            ax.set_xticks([]), ax.set_yticks([])
            ax.set_title(f"S{i}", fontsize=9)
            im_s = col_to_ma(sources[i], mask)
            im = ax.imshow(im_s, cmap=cmap, origin="upper")
            cax = inset_axes(
                ax, width="6%", height="70%", loc="lower right",
                bbox_to_anchor=(0, 0, 1, 1), bbox_transform=ax.transAxes,
                borderpad=0,
            )
            fig.colorbar(im, cax=cax)
    else:
        gs_bottom = fig.add_gridspec(
            1, 3, width_ratios=[1.3, 1.3, 2.6], wspace=0.30
        )

    # LEFT: main scatter -----------------------------------------------------
    ax_sc = fig.add_subplot(gs_bottom[0, 0])
    skw = {"s": 22, "c": "tab:blue", "alpha": 0.8, **(scatter_kw or {})}
    ax_sc.scatter(x, y, **skw)
    ax_sc.axhline(0, color="k", linewidth=1)
    ax_sc.set_title("Hover a point to preview")
    Cursor(ax_sc, useblit=True, color="grey", linewidth=0.5)

    # MIDDLE: date-pair scatter ---------------------------------------------
    ax_dt = fig.add_subplot(gs_bottom[0, 1])
    ax_dt.scatter(date1, date2, s=20, c="tab:green", alpha=0.75)
    ax_dt.set_xlim(dmin, dmax)
    ax_dt.set_ylim(dmin, dmax)      # invert y-axis (earlier dates at top)
    year_locator = mdates.YearLocator()
    year_fmt     = mdates.DateFormatter("%Y")
    for axis in (ax_dt.xaxis, ax_dt.yaxis):
        axis.set_major_locator(year_locator)
        axis.set_major_formatter(year_fmt)
    ax_dt.grid(True, linestyle=":", lw=0.6)
    ax_dt.set_xlabel("Date 1")
    ax_dt.set_ylabel("Date 2")
    ax_dt.set_title("Date pairs (yearly ticks, shared limits)")
    # highlight handle – created on first hover
    hl_dt = None

    # RIGHT: per-point previews ---------------------------------------------
    gs_prev = gs_bottom[0, 2].subgridspec(2, 2, wspace=0.07, hspace=0.25)
    ax_orig = fig.add_subplot(gs_prev[0, 0])
    ax_reco = fig.add_subplot(gs_prev[0, 1])
    ax_resi = fig.add_subplot(gs_prev[1, 0])
    ax_bar  = fig.add_subplot(gs_prev[1, 1])

    for a in (ax_orig, ax_reco, ax_resi):
        a.set_xticks([]), a.set_yticks([])

    ax_orig.set_title("Original", fontsize=10)
    ax_reco.set_title("Reconstruction", fontsize=10)
    ax_resi.set_title("Residual", fontsize=10)

    # placeholders (for fast updates) ---------------------------------------
    im_orig = im_reco = im_resi = None
    current_idx = None
    pick_radius = thresh_frac * np.hypot(np.ptp(x), np.ptp(y))

    # ── hover callback ─────────────────────────────────────────────────────-
    def on_move(event):
        nonlocal current_idx, im_orig, im_reco, im_resi, hl_dt
        if event.inaxes is not ax_sc or event.xdata is None:
            return

        dist = np.hypot(x - event.xdata, y - event.ydata)
        idx  = np.argmin(dist)
        if dist[idx] > pick_radius or idx == current_idx:
            return
        current_idx = idx

        # -- update the three image panels ----------------------------------
        img_o = col_to_ma(imgs[idx], mask)
        img_r = col_to_ma(imgs_reco[idx], mask)
        img_s = img_o - img_r

        if im_orig is None:
            im_orig = ax_orig.imshow(img_o, cmap=cmap, origin="upper")
            im_reco = ax_reco.imshow(img_r, cmap=cmap, origin="upper")
            im_resi = ax_resi.imshow(img_s, cmap=resid_cmap, origin="upper")
            for ax_i, im in zip(
                (ax_orig, ax_reco, ax_resi),
                (im_orig, im_reco, im_resi),
            ):
                cax = inset_axes(
                    ax_i, width="4%", height="70%", loc="lower right",
                    bbox_to_anchor=(0, 0, 1, 1), bbox_transform=ax_i.transAxes,
                    borderpad=0,
                )
                fig.colorbar(im, cax=cax)
        else:
            im_orig.set_data(img_o)
            im_reco.set_data(img_r)
            im_resi.set_data(img_s)
            im_orig.set_clim(np.nanmin(img_o), np.nanmax(img_o))
            im_reco.set_clim(np.nanmin(img_r), np.nanmax(img_r))
            im_resi.set_clim(np.nanmin(img_s), np.nanmax(img_s))

        # -- update the weights bar plot ------------------------------------
        ax_bar.cla()
        ax_bar.bar(np.arange(A.shape[1]), A[idx], color="tab:orange")
        ax_bar.set_title("Mixing weights", fontsize=10)
        ax_bar.set_xlabel("Component", fontsize=8)
        ax_bar.set_ylabel("Weight", fontsize=8)
        ax_bar.tick_params(axis="both", labelsize=7)
        ax_bar.set_ylim(np.nanmin(A), np.nanmax(A))

        # -- highlight the corresponding point in date scatter -------------
        if hl_dt is None:
            hl_dt = ax_dt.scatter(
                date1[idx], date2[idx], **highlight_kw
            )
        else:
            hl_dt.set_offsets([date1[idx], date2[idx]])

        fig.canvas.draw_idle()

    fig.canvas.mpl_connect("motion_notify_event", on_move)
    fig.tight_layout()
    return fig


#%% interactive_ts_viewer

def interactive_ts_viewer(ts_r3):
    """
    Display the last frame of a 3-D array.  Clicking a pixel shows its time-series.

    Parameters
    ----------
    ts_r3 : (nt, ny, nx) ndarray or MaskedArray
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import numpy.ma as ma 
    
    if ts_r3.ndim != 3:
        raise ValueError("ts_r3 must be 3-D (time, y, x)")

    nt, ny, nx = ts_r3.shape
    last_frame = ts_r3[-1]

    # ---- figure with two axes ------------------------------------------------
    fig, (ax_img, ax_ts) = plt.subplots(
        1, 2, figsize=(10, 4), gridspec_kw=dict(width_ratios=[1, 1.2])
    )

    # image
    im = ax_img.imshow(last_frame, origin="upper", cmap="viridis")
    ax_img.set_title(f"last frame (t = {nt-1})")
    ax_img.set_xlabel("x (col)")
    ax_img.set_ylabel("y (row)")
    fig.colorbar(im, ax=ax_img, fraction=0.046, pad=0.04)

    # time-series placeholder
    (line,) = ax_ts.plot([], [], marker="o")
    ax_ts.set_xlim(0, nt - 1)
    ax_ts.set_xlabel("time index")
    ax_ts.set_ylabel("value")
    ax_ts.set_title("click a pixel → its time-series")

    # ---- click handler -------------------------------------------------------
    def onclick(event):
        # disregard clicks outside the image axes
        if event.inaxes is not ax_img or event.xdata is None or event.ydata is None:
            return

        col = int(round(event.xdata))
        row = int(round(event.ydata))
        if not (0 <= row < ny and 0 <= col < nx):
            return

        ts = ts_r3[:, row, col]
        ts = ts.filled(np.nan) if ma.isMaskedArray(ts) else ts

        line.set_data(np.arange(nt), ts)
        ax_ts.set_ylim(np.nanmin(ts) - 0.1, np.nanmax(ts) + 0.1)
        ax_ts.set_title(f"time-series at (row={row}, col={col})")
        fig.canvas.draw_idle()

    cid = fig.canvas.mpl_connect("button_press_event", onclick)
    plt.tight_layout()
    return cid  # returning the connection id lets you disconnect if desired



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

