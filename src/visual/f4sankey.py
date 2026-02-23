import plotly.graph_objects as go
from visual.utils import *

# Load data
data_dir = "/home/wjiang49/group/wjiang49/data/traceCB/EAS_eQTLGen/coloc"
save_path = "/home/wjiang49/group/wjiang49/data/traceCB/EAS_eQTLGen/results/coloc"

prefix = "bcx"
# Read all files ending with _coloc.csv
coloc_files = glob.glob(f"{data_dir}/{prefix}*_coloc.csv")
print(f"Found {len(coloc_files)} coloc files")

# Merge all data (filter out empty DataFrames to avoid FutureWarning)
df_list = []
for f in coloc_files:
    tmp_df = pd.read_csv(f)
    if not tmp_df.empty:
        df_list.append(tmp_df)
df = pd.concat(df_list, ignore_index=True)
print(f"Total rows before filtering: {len(df)}")

# Threshold is 0.7
threshold = 0.7

# Filter: Keep genes with h3 or h4 > 0.7 in at least one condition
has_signal = (
    (df["pp_h3_original"] > threshold)
    | (df["p_original"] > threshold)
    | (df["pp_h3_traceC"] > threshold)
    | (df["p_traceC"] > threshold)
    | (df["pp_h3_traceCB"] > threshold)
    | (df["p_traceCB"] > threshold)
)
df = df[has_signal]
print(f"Total rows after filtering (at least one h3/h4 > 0.7): {len(df)}")

# Calculate the count of h3, h4, and others (1-h3-h4) for each method
# h3 > 0.7
h3_original = (df["pp_h3_original"] > threshold).sum()
h3_traceC = (df["pp_h3_traceC"] > threshold).sum()
h3_traceCB = (df["pp_h3_traceCB"] > threshold).sum()

# h4 (p) > 0.7
h4_original = (df["p_original"] > threshold).sum()
h4_traceC = (df["p_traceC"] > threshold).sum()
h4_traceCB = (df["p_traceCB"] > threshold).sum()

# Other (neither h3 nor h4 > 0.7)
other_original = (
    (df["pp_h3_original"] <= threshold) & (df["p_original"] <= threshold)
).sum()
other_traceC = ((df["pp_h3_traceC"] <= threshold) & (df["p_traceC"] <= threshold)).sum()
other_traceCB = (
    (df["pp_h3_traceCB"] <= threshold) & (df["p_traceCB"] <= threshold)
).sum()

print("=== Statistics ===")
print(f"Original: H3={h3_original}, H4={h4_original}, Other={other_original}")
print(f"traceC:   H3={h3_traceC}, H4={h4_traceC}, Other={other_traceC}")
print(f"traceCB:  H3={h3_traceCB}, H4={h4_traceCB}, Other={other_traceCB}")


# Calculate state categories (state of each row)
def get_category(df, h3_col, h4_col, threshold):
    """Return category for each row: 0=COLOC(H4), 1=Independent(H3), 2=Other"""
    is_h3 = df[h3_col] > threshold
    is_h4 = df[h4_col] > threshold
    # Priority: H4(COLOC) > H3(Independent) > Other
    category = pd.Series(2, index=df.index)  # Default Other
    category[is_h3] = 1  # Independent (Original H3)
    category[is_h4] = 0  # COLOC (Original H4)
    return category


cat_original = get_category(df, "pp_h3_original", "p_original", threshold)
cat_traceC = get_category(df, "pp_h3_traceC", "p_traceC", threshold)
cat_traceCB = get_category(df, "pp_h3_traceCB", "p_traceCB", threshold)

# Calculate transition matrix: Original -> traceC
# source: Original's COLOC(0), Independent(1), Other(2)
# target: traceC's COLOC(3), Independent(4), Other(5)
flow_orig_to_traceC = pd.crosstab(cat_original, cat_traceC)
print("\n=== Flow: Original -> traceC ===")
print(flow_orig_to_traceC)

# Calculate transition matrix: traceC -> traceCB
# source: traceC's COLOC(3), Independent(4), Other(5)
# target: traceCB's COLOC(6), Independent(7), Other(8)
flow_traceC_to_traceCB = pd.crosstab(cat_traceC, cat_traceCB)
print("\n=== Flow: traceC -> traceCB ===")
print(flow_traceC_to_traceCB)

# Calculate count for each node
orig_counts = cat_original.value_counts().sort_index()
traceC_counts = cat_traceC.value_counts().sort_index()
traceCB_counts = cat_traceCB.value_counts().sort_index()

# Define nodes (with count labels)
# Left column: Original's 3 categories (0,1,2)
# Middle column: traceC's 3 categories (3,4,5)
# Right column: traceCB's 3 categories (6,7,8)
category_names = ["COLOC", "Independent", "Undetermined"]
nodes = [
    f"COLOC ({orig_counts.get(0, 0)})",
    f"Independent ({orig_counts.get(1, 0)})",
    f"Undetermined ({orig_counts.get(2, 0)})",  # Original: 0, 1, 2
    f"COLOC ({traceC_counts.get(0, 0)})",
    f"Independent ({traceC_counts.get(1, 0)})",
    f"Undetermined ({traceC_counts.get(2, 0)})",  # traceC: 3, 4, 5
    f"COLOC ({traceCB_counts.get(0, 0)})",
    f"Independent ({traceCB_counts.get(1, 0)})",
    f"Undetermined ({traceCB_counts.get(2, 0)})",  # traceCB: 6, 7, 8
]

# Build links
sources = []
targets = []
values = []
colors_list = []

# Define colors for different categories
coloc_color = "rgb(124, 200, 124)"  # Green (COLOC)
independent_color = "rgb(150, 181, 248)"  # Blue (Independent)
undetermined_color = "rgb(235, 160, 152)"  # Pink (Undetermined)
category_colors = [coloc_color, independent_color, undetermined_color]


# Link colors (semi-transparent)
def color_to_rgba(color, alpha=0.5):
    """Convert hex or rgb color to rgba string"""
    if color.startswith("rgb("):
        # rgb(r, g, b) format
        rgb_values = color[4:-1].split(",")
        r, g, b = [int(v.strip()) for v in rgb_values]
    else:
        # hex format
        hex_color = color.lstrip("#")
        r, g, b = tuple(int(hex_color[i : i + 2], 16) for i in (0, 2, 4))
    return f"rgba({r}, {g}, {b}, {alpha})"


# Flow: Original -> traceC
for src_cat in range(3):  # COLOC, Independent, Undetermined
    for tgt_cat in range(3):
        if (
            src_cat in flow_orig_to_traceC.index
            and tgt_cat in flow_orig_to_traceC.columns
        ):
            val = flow_orig_to_traceC.loc[src_cat, tgt_cat]
            if val > 0:
                sources.append(src_cat)  # Original nodes: 0, 1, 2
                targets.append(tgt_cat + 3)  # traceC nodes: 3, 4, 5
                values.append(val)
                colors_list.append(color_to_rgba(category_colors[src_cat]))

# Flow: traceC -> traceCB
for src_cat in range(3):  # COLOC, Independent, Undetermined
    for tgt_cat in range(3):
        if (
            src_cat in flow_traceC_to_traceCB.index
            and tgt_cat in flow_traceC_to_traceCB.columns
        ):
            val = flow_traceC_to_traceCB.loc[src_cat, tgt_cat]
            if val > 0:
                sources.append(src_cat + 3)  # traceC nodes: 3, 4, 5
                targets.append(tgt_cat + 6)  # traceCB nodes: 6, 7, 8
                values.append(val)
                colors_list.append(color_to_rgba(category_colors[src_cat]))

links = {
    "source": sources,
    "target": targets,
    "value": values,
}

node_colors = [
    coloc_color,
    independent_color,
    undetermined_color,  # Original
    coloc_color,
    independent_color,
    undetermined_color,  # traceC
    coloc_color,
    independent_color,
    undetermined_color,  # traceCB
]

link_colors = colors_list

# Create Sankey diagram
fig = go.Figure(
    data=[
        go.Sankey(
            arrangement="snap",
            node=dict(
                pad=20,
                thickness=25,
                line=dict(color="white", width=1),
                label=nodes,
                color=node_colors,
                x=[
                    0.01,
                    0.01,
                    0.01,
                    0.5,
                    0.5,
                    0.5,
                    0.99,
                    0.99,
                    0.99,
                ],  # Manually position 3 columns
                y=[
                    0.15,
                    0.5,
                    0.85,
                    0.15,
                    0.5,
                    0.85,
                    0.15,
                    0.5,
                    0.85,
                ],  # Manually position 3 rows
            ),
            link=dict(
                source=links["source"],
                target=links["target"],
                value=links["value"],
                color=link_colors,
            ),
        )
    ]
)

# Update layout
yloc = 1.04
fig.update_layout(
    title={
        "text": "Colocalization Pattern Flow",
        "x": 0.5,
        "xanchor": "center",
        "font": {"size": 16, "family": "Arial", "color": "#333"},
    },
    font=dict(size=13, family="Arial", color="#333"),
    plot_bgcolor="white",
    paper_bgcolor="white",
    width=550,
    height=380,
    margin=dict(l=5, r=5, t=70, b=10),
    annotations=[
        dict(
            text="Original",
            x=-0.01,
            y=yloc,
            xref="paper",
            yref="paper",
            showarrow=False,
            font=dict(size=14, color="#555", family="Arial"),
        ),
        dict(
            text="traceC",
            x=0.5,
            y=yloc,
            xref="paper",
            yref="paper",
            showarrow=False,
            font=dict(size=14, color="#555", family="Arial"),
        ),
        dict(
            text="traceCB",
            x=1.01,
            y=yloc,
            xref="paper",
            yref="paper",
            showarrow=False,
            font=dict(size=14, color="#555", family="Arial"),
        ),
    ],
)

# Save as HTML
fig.write_html(f"{save_path}/{prefix}_sankey_by_category.html")
print(f"Sankey diagram saved to {save_path}/{prefix}_sankey_by_category.html")

# Save as static image
try:
    fig.write_image(
        f"{save_path}/{prefix}_sankey_by_category.pdf", width=550, height=380, scale=2
    )
    print(f"PDF saved to {save_path}/{prefix}_sankey_by_category.pdf")
except Exception as e:
    print(f"Could not save images: {e}")

# Show plot (uncomment if interactive view is needed)
# fig.show()
