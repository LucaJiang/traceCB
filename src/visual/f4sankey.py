import plotly.graph_objects as go
from visual.utils import *

# 加载数据
data_dir = "/home/wjiang49/group/wjiang49/data/traceCB/EAS_eQTLGen/coloc"
save_path = "/home/wjiang49/group/wjiang49/data/traceCB/EAS_eQTLGen/results/coloc"

prefix = "bcx"
# 读取所有以_coloc.csv结尾的文件
coloc_files = glob.glob(f"{data_dir}/{prefix}*_coloc.csv")
print(f"Found {len(coloc_files)} coloc files")

# 合并所有数据（过滤掉空的DataFrame以避免FutureWarning）
df_list = []
for f in coloc_files:
    tmp_df = pd.read_csv(f)
    if not tmp_df.empty:
        df_list.append(tmp_df)
df = pd.concat(df_list, ignore_index=True)
print(f"Total rows before filtering: {len(df)}")

# 阈值为0.7
threshold = 0.7

# 过滤：只保留在至少一个条件下有h3或h4>0.7的gene
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

# 计算每个方法的h3, h4, 和 其他（1-h3-h4）的数量
# h3 > 0.7
h3_original = (df["pp_h3_original"] > threshold).sum()
h3_traceC = (df["pp_h3_traceC"] > threshold).sum()
h3_traceCB = (df["pp_h3_traceCB"] > threshold).sum()

# h4 (p) > 0.7
h4_original = (df["p_original"] > threshold).sum()
h4_traceC = (df["p_traceC"] > threshold).sum()
h4_traceCB = (df["p_traceCB"] > threshold).sum()

# 其他 (既不是h3也不是h4 > 0.7)
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


# 计算状态类别（每行数据的状态）
def get_category(df, h3_col, h4_col, threshold):
    """返回每行的类别: 0=COLOC(H4), 1=Independent(H3), 2=Other"""
    is_h3 = df[h3_col] > threshold
    is_h4 = df[h4_col] > threshold
    # 优先级: H4(COLOC) > H3(Independent) > Other
    category = pd.Series(2, index=df.index)  # 默认Other
    category[is_h3] = 1  # Independent (原H3)
    category[is_h4] = 0  # COLOC (原H4)
    return category


cat_original = get_category(df, "pp_h3_original", "p_original", threshold)
cat_traceC = get_category(df, "pp_h3_traceC", "p_traceC", threshold)
cat_traceCB = get_category(df, "pp_h3_traceCB", "p_traceCB", threshold)

# 计算转换矩阵: Original -> traceC
# source: Original的COLOC(0), Independent(1), Other(2)
# target: traceC的COLOC(3), Independent(4), Other(5)
flow_orig_to_traceC = pd.crosstab(cat_original, cat_traceC)
print("\n=== Flow: Original -> traceC ===")
print(flow_orig_to_traceC)

# 计算转换矩阵: traceC -> traceCB
# source: traceC的COLOC(3), Independent(4), Other(5)
# target: traceCB的COLOC(6), Independent(7), Other(8)
flow_traceC_to_traceCB = pd.crosstab(cat_traceC, cat_traceCB)
print("\n=== Flow: traceC -> traceCB ===")
print(flow_traceC_to_traceCB)

# 计算每个节点的数量
orig_counts = cat_original.value_counts().sort_index()
traceC_counts = cat_traceC.value_counts().sort_index()
traceCB_counts = cat_traceCB.value_counts().sort_index()

# 定义节点（带数量标签）
# 左列: Original的3个类别 (0,1,2)
# 中列: traceC的3个类别 (3,4,5)
# 右列: traceCB的3个类别 (6,7,8)
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

# 构建链接
sources = []
targets = []
values = []
colors_list = []

# 为不同类别定义颜色
coloc_color = "rgb(124, 200, 124)"  # 绿色 (COLOC)
independent_color = "rgb(150, 181, 248)"  # 蓝色 (Independent)
undetermined_color = "rgb(235, 160, 152)"  # 粉红色 (Undetermined)
category_colors = [coloc_color, independent_color, undetermined_color]


# 链接颜色（半透明）
def color_to_rgba(color, alpha=0.5):
    """Convert hex or rgb color to rgba string"""
    if color.startswith("rgb("):
        # rgb(r, g, b) 格式
        rgb_values = color[4:-1].split(",")
        r, g, b = [int(v.strip()) for v in rgb_values]
    else:
        # hex 格式
        hex_color = color.lstrip("#")
        r, g, b = tuple(int(hex_color[i : i + 2], 16) for i in (0, 2, 4))
    return f"rgba({r}, {g}, {b}, {alpha})"


# Original -> traceC 的流动
for src_cat in range(3):  # COLOC, Independent, Undetermined
    for tgt_cat in range(3):
        if (
            src_cat in flow_orig_to_traceC.index
            and tgt_cat in flow_orig_to_traceC.columns
        ):
            val = flow_orig_to_traceC.loc[src_cat, tgt_cat]
            if val > 0:
                sources.append(src_cat)  # Original节点: 0, 1, 2
                targets.append(tgt_cat + 3)  # traceC节点: 3, 4, 5
                values.append(val)
                colors_list.append(color_to_rgba(category_colors[src_cat]))

# traceC -> traceCB 的流动
for src_cat in range(3):  # COLOC, Independent, Undetermined
    for tgt_cat in range(3):
        if (
            src_cat in flow_traceC_to_traceCB.index
            and tgt_cat in flow_traceC_to_traceCB.columns
        ):
            val = flow_traceC_to_traceCB.loc[src_cat, tgt_cat]
            if val > 0:
                sources.append(src_cat + 3)  # traceC节点: 3, 4, 5
                targets.append(tgt_cat + 6)  # traceCB节点: 6, 7, 8
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

# 创建Sankey图
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
                x=[0.01, 0.01, 0.01, 0.5, 0.5, 0.5, 0.99, 0.99, 0.99],  # 手动定位3列
                y=[0.15, 0.5, 0.85, 0.15, 0.5, 0.85, 0.15, 0.5, 0.85],  # 手动定位3行
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

# 更新布局
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

# 保存为HTML
fig.write_html(f"{save_path}/{prefix}_sankey_by_category.html")
print(f"Sankey diagram saved to {save_path}/{prefix}_sankey_by_category.html")

# 保存为静态图片
try:
    fig.write_image(
        f"{save_path}/{prefix}_sankey_by_category.pdf", width=550, height=380, scale=2
    )
    print(f"PDF saved to {save_path}/{prefix}_sankey_by_category.pdf")
except Exception as e:
    print(f"Could not save images: {e}")

# 显示图表（如果需要交互查看可以取消注释）
# fig.show()
