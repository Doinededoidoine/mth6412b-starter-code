import Base.show

include("node.jl")
include("edge.jl")

"""Type abstrait dont d'autres types de graphes dériveront."""
abstract type AbstractGraph{T} end

"""Type representant un graphe comme un ensemble de noeuds.

Exemple :

    node1 = Node("Joe", 3.14)
    node2 = Node("Steve", exp(1))
    node3 = Node("Jill", 4.12)
    G = Graph("Ick", [node1, node2, node3])

Attention, tous les noeuds doivent avoir des données de même type.
"""
mutable struct Graph{T} <: AbstractGraph{T}
  name::String
  nodes::Vector{Node{T}}
  edges::Vector{Edge{T}}
end

"""Ajoute un noeud au graphe."""
function add_node!(graph::Graph{T}, node::Node{T}) where T
  push!(graph.nodes, node)
  graph
end

"""Ajoute une arête au graphe."""
function add_edge!(graph::Graph{T}, edge::Edge{T}) where T
  push!(graph.edges, edge)
  graph
end

"""Retire une arête au graphe."""
function pop_edge!(graph::Graph{T}, i::Int64) where T
  deleteat!(graph.edges, i)
  graph
end

# on présume que tous les graphes dérivant d'AbstractGraph
# posséderont des champs `name` et `nodes`.

"""Renvoie le nom du graphe."""
name(graph::AbstractGraph) = graph.name

"""Renvoie la liste des noeuds du graphe."""
nodes(graph::AbstractGraph) = graph.nodes

"""Renvoie le nombre de noeuds du graphe."""
nb_nodes(graph::AbstractGraph) = length(graph.nodes)

"""Renvoie la liste des arcs du graphe."""
edges(graph::AbstractGraph) = graph.edges

"""Renvoie le nombre d'arcs du graphe."""
nb_edges(graph::AbstractGraph) = length(graph.edges)

"""Affiche un graphe"""
function show(graph::Graph)
  println("Graph ", name(graph), " has ", nb_nodes(graph), " nodes and ", nb_edges(graph), " edges.")
  for node in nodes(graph)
    show(node)
  end
  for edge in edges(graph)
    show(edge)
  end
end

"""Trouver le poids de l'arc reliant ces deux noeuds"""
function findweight(graph::AbstractGraph, node1::AbstractNode, node2::AbstractNode)
  try
    return weight(filter(edge -> name(s_node(edge)) == name(node1) && name(d_node(edge)) == name(node2) || name(s_node(edge)) == name(node2) && name(d_node(edge)) == name(node1), edges(graph))[1])
  catch e
    error("Aucune arête entre les noeuds ", name(node1), " et ", name(node2), " n'a été trouvée dans le graphe ", name(graph))
  end
end

"""Trouver le poids du graphe"""
function graphweight(graph::AbstractGraph)
  weightgraph = 0
  for edge in edges(graph)
    weightgraph += weight(edge)
  end
  return weightgraph
end

"""Afficher si un noeud se trouve dans un graphe"""
function find_node(graph::AbstractGraph, node::AbstractNode)
  return name(node) in name.(nodes(graph))
end

"""Afficher le noeud le plus lourd du graph"""
function find_heaviest_node(graph::AbstractGraph)
  edge_max = []
  weight_max = -Inf
  for edge in edges(graph)
    if weight(edge) > weight_max
      edge_max = edge
      weight_max = weight(edge_max)
    end
  end
  return s_node(edge_max)
end

"""Afficher le noeud le plus leger du graph"""
function find_lightest_node(graph::AbstractGraph)
  edge_min = []
  weight_min = Inf
  for edge in edges(graph)
    if weight(edge) < weight_min
      edge_min = edge
      weight_min = weight(edge_min)
    end
  end
  return s_node(edge_min)
end
