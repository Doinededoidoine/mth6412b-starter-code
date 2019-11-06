include("node.jl")
include("edge.jl")
include("graph.jl")
include("read_stsp.jl")
using Plots

""" Analyse une instance de TSP symétrique dont les poids sont donnés au format
 EXPLICIT et construit un objet de type Graph.
"""
function main1(filename::String)
  header = read_header(filename)
  graph_edges = read_edges(header, filename)
  graph_nodes = read_nodes(header, filename)
  nodes_empty = Vector{Node}()
  edges_empty = Vector{Edge}()
  graph = Graph{Array{Float64, 1}}("graph", nodes_empty, edges_empty)
  for i = 1 : length(graph_nodes)
    new_node = Node{Array{Float64, 1}}(string(collect(keys(graph_nodes))[i]), graph_nodes[i], Inf)
    add_node!(graph, new_node)
  end
  for edge in graph_edges
    s_node = Node{Array{Float64, 1}}("", [], Inf)
    d_node = Node{Array{Float64, 1}}("", [], Inf)
    for node in nodes(graph)
      if name(node) == string(edge[1])
        s_node = node
      end
      if name(node) == string(edge[2])
        d_node = node
      end
    end
    new_edge = Edge{Array{Float64, 1}}("", s_node, d_node, parse(Int64, edge[3]))
    add_edge!(graph, new_edge)
  end
  # On affiche puis on renvoie le graphe.
  show(graph)
  graph
end
