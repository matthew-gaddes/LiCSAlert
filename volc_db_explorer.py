#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 21 15:51:54 2025

@author: matthew
"""

import pandas as pd
import numpy as np
import pdb

import cartopy.crs as ccrs
import mplcursors
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.lines import Line2D


#%% Plot volcano map


def plot_volcano_map(
    df: pd.DataFrame,
    *,
    include_priorities=None,          # e.g. ["A1", "A2"]
    figsize=(20, 11),
    projection=ccrs.Robinson(),
    point_size=40,
    interactive=True
    ):
    """
    Plot volcano locations on a global map, coloured by priority.

    Parameters
    ----------
    df : DataFrame with columns 'lat', 'lon', 'name', 'priority', 'alt'
    include_priorities : list or None
        Subset of priority classes to plot (e.g. ["A1", "A2"]).
        None → plot all.
    figsize, projection, point_size : map aesthetics.

    Returns
    -------
    fig, ax
    """

    # ------------------------------------------------------------------
    # 0. Priority palette & dtype
    # ------------------------------------------------------------------
    priority_order  = ["A1", "A2", "B1", "B2", "C", "None"]
    priority_colors = ["#d73027", "#fc8d59", "#fee08b",
                       "#d9ef8b", "#91cf60", "#cccccc"]

    # Make a safe copy; normalise nulls to the string "None"
    df = df.copy()
    df["priority"] = (
        df["priority"]
          .where(~df["priority"].isna(), "None")
          .astype(str).str.strip()
          .astype(pd.CategoricalDtype(priority_order, ordered=True))
    )

    # Filter to requested priorities
    if include_priorities is not None:
        include_priorities = [str(p).strip() for p in include_priorities]
        mask = df["priority"].isin(include_priorities)
        df = df.loc[mask]

    # Mapping objects
    pri_to_idx = {p: i for i, p in enumerate(priority_order)}
    cmap = ListedColormap(priority_colors)
    norm = BoundaryNorm(np.arange(-0.5, len(priority_order) + 0.5), cmap.N)

    # ------------------------------------------------------------------
    # 1. Map setup
    # ------------------------------------------------------------------
    fig = plt.figure(figsize=figsize)
    ax  = fig.add_subplot(1, 1, 1, projection=projection)
    ax.set_global()
    ax.stock_img()
    ax.coastlines()

    # ------------------------------------------------------------------
    # 2. Scatter volcanoes
    # ------------------------------------------------------------------
    sc = ax.scatter(
        df["lon"], df["lat"],
        c=df["priority"].map(pri_to_idx),
        cmap=cmap, norm=norm,
        s=point_size, edgecolor="k", linewidth=0.3,
        transform=ccrs.PlateCarree()
    )

    # ------------------------------------------------------------------
    # 3. Interactive tooltip
    # ------------------------------------------------------------------
    if interactive:
        cursor = mplcursors.cursor(sc, hover=True)
        
        #approach 1: explicit registration
        # cursor.connect("add")(_on_add)   # explicit registration
    
        # approach 2: use a decorator to call following function from cursor.connect
        @cursor.connect("add")
        def _on_add(sel):
            row = df.iloc[sel.index]
            sel.annotation.set_text(
                f"{row['name']}\nPriority: {row['priority']}\n"
            )
            sel.annotation.get_bbox_patch().set_alpha(0.8)

    # ------------------------------------------------------------------
    # 4. Legend (full palette, but greyed handles for filtered-out classes)
    # ------------------------------------------------------------------
    handles = []
    for p, col in zip(priority_order, priority_colors):
        visible = p in df["priority"].unique()
        handles.append(
            Line2D([0], [0], marker='o', color='w',
                   markerfacecolor=col if visible else "#ffffff",
                   markeredgecolor='k', markersize=8,
                   label=p if p != "None" else "None/NaN",
                   alpha=1.0 if visible else 0.25)
        )
    ax.legend(handles=handles, title="Priority", loc="lower left")

    ax.set_title("Global Volcano Priorities", pad=20)
    return fig, ax



#%% Load volcano database
# ------------------------------------------------------------------
# 1. Load the DataFrame from the pickle
# ------------------------------------------------------------------
df = pd.read_pickle("comet_volcano_database.pkl")

# ------------------------------------------------------------------
# 2. Headline info
# ------------------------------------------------------------------
print("\nShape (rows, cols):", df.shape)
print("\nColumn names:", df.columns.tolist())
print("\nData types:\n", df.dtypes)

# basic memory footprint
print("\nMemory usage:")
df.info(memory_usage="deep", show_counts=True)

# ------------------------------------------------------------------
# 3. Quick numeric summary
# ------------------------------------------------------------------
print("\nDescriptive stats (numeric columns):")
print(df.describe(include=[np.number]).T)

# ------------------------------------------------------------------
# 4. Peek at the data
# ------------------------------------------------------------------
print("\nFirst five rows:")
print(df.head())

print("\nRandom sample:")
print(df.sample(5, random_state=0))

# ------------------------------------------------------------------
# 5. Missing-value overview
# ------------------------------------------------------------------
missing = df.isna().mean().sort_values(ascending=False)
print("\nMissing-value ratio by column (top 10):")
print(missing.head(10))

# ------------------------------------------------------------------
# 6. Simple correlation matrix (numeric cols),
#    limited to the first 10 columns for readability
# ------------------------------------------------------------------
numeric_cols = df.select_dtypes("number").columns[:10]
print("\nCorrelation (first 10 numeric columns):")
print(df[numeric_cols].corr().round(2))

# ------------------------------------------------------------------
# 7. If you prefer an interactive report (optional):
# ------------------------------------------------------------------
# pip install pandas-profiling or ydata-profiling first
# from pandas_profiling import ProfileReport
# profile = ProfileReport(df, title="VTB Data Overview")
# profile.to_file("vtb_profile.html")
# print("Interactive profile saved to vtb_profile.html")



#%% Interactive  map for different priorities



fig, ax = plot_volcano_map(
    df, 
    include_priorities=['A1'],
    interactive=True)
fig.suptitle(f"Included priorities: A1", fontsize=16, y=0.95)
plt.show()



#%% map for different priorities (interactive off)


# priority_steps = [
#     ["A1"],
#     ["A1", "A2"],
#     ["A1", "A2", "B1"],
#     ["A1", "A2", "B1", "B2"],
#     ["A1", "A2", "B1", "B2", "C"],
#     None                   # final map with every class including None/NaN
# ]

# for step in priority_steps:
#     fig, ax = plot_volcano_map(
#         df, 
#         include_priorities=step,
#         interactive=True)
#     fig.suptitle(f"Included priorities: {step or 'All'}", fontsize=16, y=0.95)
#     plt.show()