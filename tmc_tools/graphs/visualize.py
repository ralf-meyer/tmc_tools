import networkx as nx
from tmc_tools.constants import atomic_numbers, jmol_colors


def draw_graph(graph, **kwargs):
    node_color = [
        jmol_colors[atomic_numbers[sym]] for _, sym in graph.nodes().data("symbol")
    ]
    nx.draw(graph, node_color=node_color, edgecolors="k", **kwargs)
