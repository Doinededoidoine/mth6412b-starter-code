import Base.show
import Base.copy

include(joinpath(@__DIR__, "..", "phase1", "graph.jl"))

abstract type AbstractDegreeTable{T} end

"""Structure de données associant à chaque noeud d'un graphe son degré."""
mutable struct DegreeTable{T} <: AbstractDegreeTable{T}
    nodes::Vector{Node{T}}
    degrees::Vector{Int64}
end

DegreeTable{T}() where T = DegreeTable(T[])

"""Renvoie la liste des noeuds."""
nodes(degree_table::AbstractDegreeTable) = degree_table.nodes

"""Renvoie la liste des degrés des noeuds."""
degrees(degree_table::AbstractDegreeTable) = degree_table.degrees

"""Initialise la table des degrés pour un graphe donné."""
function init_degree_table(graph::AbstractGraph{T}) where T
    return DegreeTable{T}(nodes(graph), zeros(Int, nb_nodes(graph)))
end

"""Attribue à chaque noeud son degré dans le graphe en entrée."""
function d_table(graph::AbstractGraph)
    degree_table = init_degree_table(graph)
    for edge in edges(graph)
        degrees(degree_table)[findfirst(n -> n == name(s_node(edge)), name.(nodes(degree_table)))] += 1
        degrees(degree_table)[findfirst(n -> n == name(d_node(edge)), name.(nodes(degree_table)))] += 1
    end
    degree_table
end
