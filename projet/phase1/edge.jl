import Base.show

include("node.jl")

"""Type abstrait dont d'autres types d'arcs dériveront."""
abstract type AbstractEdge{T} end

"""Type représentant les arcs d'un graphe.

Exemple:

        arc = Edge("James", [π, exp(1)])
        arc = Edge("Kirk", "guitar")
        arc = Edge("Lars", 2)

"""
mutable struct Edge{T} <: AbstractEdge{T}
  name::String
  s_node::AbstractNode{T}
  d_node::AbstractNode{T}
  weight::Float64
end

# on présume que tous les arcs dérivant d'AbstractEdge
# posséderont des champs `name`, `node1`, `node2` et `weight`.

"""Renvoie le nom de l'arc."""
name(edge::AbstractEdge) = edge.name

"""Renvoie le noeud de départ de l'arc."""
s_node(edge::AbstractEdge) = edge.s_node

"""Renvoie le noeud d'arrivée de l'arc."""
d_node(edge::AbstractEdge) = edge.d_node

"""Renvoie le poids de l'arc."""
weight(edge::AbstractEdge) = edge.weight

function set_weight!(edge::AbstractEdge, new_weight::Float64)
  edge.weight = new_weight
  edge
end

"""Affiche un arc."""
function show(edge::AbstractEdge)
  println("Edge ", name(edge), " between nodes ", name(s_node(edge)), " and ", name(d_node(edge)), ", weight: ", weight(edge))
end
