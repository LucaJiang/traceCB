import pandas as pd
import plotly.graph_objects as go
from visual.utils import *

# 加载数据
data_path = (
    "/Users/lucajiang/learn/CityU/traceCB/data/coloc/coloc_statistics_summary.csv"
)
save_path = "/Users/lucajiang/learn/CityU/traceCB/data/img/eas_eqtlgen"

# 读取数据
df = pd.read_csv(data_path)

# 计算每列的总和
original_sum = df["p_original_gt0.7"].sum()
traceC_sum = df["p_traceC_gt0.7"].sum()
traceCB_sum = df["p_traceCB_gt0.7"].sum()

print(f"Original: {original_sum}")
print(f"traceC: {traceC_sum}")
print(f"traceCB: {traceCB_sum}")

# 定义节点 - 使用meta_data中的方法名称
nodes = meta_data["method_name"]  # ["Original", "traceC", "traceCB"]

# 定义链接（流）
links = {
    "source": [0, 1],  # Original -> traceC, traceC -> traceCB
    "target": [1, 2],  # traceC, traceCB
    "value": [
        traceC_sum,  # Original -> traceC 的流量
        traceCB_sum,  # traceC -> traceCB 的流量
    ],
}

# 从meta_data获取颜色
colors = meta_data["Colors"]  # 使用大写C的Colors
node_colors = [colors[method] for method in meta_data["method_name"]]


# 链接颜色（半透明）- 从RGB转换为RGBA
def hex_to_rgba(hex_color, alpha=0.4):
    """Convert hex color to rgba string"""
    hex_color = hex_color.lstrip("#")
    r, g, b = tuple(int(hex_color[i : i + 2], 16) for i in (0, 2, 4))
    return f"rgba({r}, {g}, {b}, {alpha})"


link_colors = [
    hex_to_rgba(colors[meta_data["method_name"][1]], 0.4),  # traceC
    hex_to_rgba(colors[meta_data["method_name"][2]], 0.4),  # traceCB
]

# 创建Sankey图
fig = go.Figure(
    data=[
        go.Sankey(
            node=dict(
                pad=15,
                thickness=20,
                line=dict(color="black", width=0.5),
                label=nodes,
                color=node_colors,
                customdata=[original_sum, traceC_sum, traceCB_sum],
                hovertemplate="%{label}<br>Count: %{customdata}<extra></extra>",
            ),
            link=dict(
                source=links["source"],
                target=links["target"],
                value=links["value"],
                color=link_colors,
                hovertemplate="%{source.label} → %{target.label}<br>Count: %{value}<extra></extra>",
            ),
        )
    ]
)

# 更新布局
fig.update_layout(
    title={
        "text": "Colocalization Evidence Flow",
        "x": 0.5,
        "xanchor": "center",
        "font": {"size": 20, "family": "Arial, sans-serif"},
    },
    font=dict(size=14, family="Arial, sans-serif"),
    plot_bgcolor="white",
    paper_bgcolor="white",
    width=800,
    height=500,
    margin=dict(l=20, r=20, t=60, b=20),
)

# 保存为HTML
fig.write_html(f"{save_path}/sankey_p_h4.html")
print(f"Sankey diagram saved to {save_path}/sankey_p_h4.html")

# 保存为静态图片（需要安装kaleido）
try:
    fig.write_image(f"{save_path}/sankey_p_h4.pdf", width=800, height=500)
    print(f"PDF saved to {save_path}/sankey_p_h4.pdf")
except Exception as e:
    print(f"Could not save PDF: {e}")
    print("Install kaleido with: pip install kaleido")

# 显示图表
fig.show()

# 同样为 pp_h3 创建 Sankey 图
h3_original_sum = df["pp_h3_original_gt0.7"].sum()
h3_traceC_sum = df["pp_h3_traceC_gt0.7"].sum()
h3_traceCB_sum = df["pp_h3_traceCB_gt0.7"].sum()

print(f"\nH3 Original: {h3_original_sum}")
print(f"H3 traceC: {h3_traceC_sum}")
print(f"H3 traceCB: {h3_traceCB_sum}")

links_h3 = {
    "source": [0, 1],
    "target": [1, 2],
    "value": [h3_traceC_sum, h3_traceCB_sum],
}

fig_h3 = go.Figure(
    data=[
        go.Sankey(
            node=dict(
                pad=15,
                thickness=20,
                line=dict(color="black", width=0.5),
                label=nodes,
                color=node_colors,
                customdata=[h3_original_sum, h3_traceC_sum, h3_traceCB_sum],
                hovertemplate="%{label}<br>Count: %{customdata}<extra></extra>",
            ),
            link=dict(
                source=links_h3["source"],
                target=links_h3["target"],
                value=links_h3["value"],
                color=link_colors,
                hovertemplate="%{source.label} → %{target.label}<br>Count: %{value}<extra></extra>",
            ),
        )
    ]
)

fig_h3.update_layout(
    title={
        "text": "Colocalization Evidence Flow: PP.H3 > 0.7",
        "x": 0.5,
        "xanchor": "center",
        "font": {"size": 20, "family": "Arial, sans-serif"},
    },
    font=dict(size=14, family="Arial, sans-serif"),
    plot_bgcolor="white",
    paper_bgcolor="white",
    width=800,
    height=500,
    margin=dict(l=20, r=20, t=60, b=20),
)

fig_h3.write_html(f"{save_path}/sankey_pp_h3.html")
print(f"H3 Sankey diagram saved to {save_path}/sankey_pp_h3.html")

try:
    fig_h3.write_image(f"{save_path}/sankey_pp_h3.pdf", width=800, height=500)
    print(f"H3 PDF saved to {save_path}/sankey_pp_h3.pdf")
except Exception as e:
    print(f"Could not save H3 PDF: {e}")

fig_h3.show()
