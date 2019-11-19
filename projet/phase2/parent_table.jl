import Base.show
import Base.copy

include(joinpath(@__DIR__, "..", "phase1", "graph.jl"))

abstract type AbstractParentTable{T} end

"""Structure de données associant à chaque noeud d'un graphe son noeud parent."""
mutable struct ParentTable{T} <: AbstractParentTable{T}
    enfants::Vector{Node{T}}
    parents::Vector{Node{T}}
    ranks::Vector{Int64}
end

"""Renvoie la liste des noeuds."""
enfants(parent_table::AbstractParentTable) = parent_table.enfants

"""Renvoie la liste des parents de chaque noeud."""
parents(parent_table::AbstractParentTable) = parent_table.parents

"""Renvoie le rang associé à un noeud"""
ranks(parent_table::AbstractParentTable) = parent_table.ranks

"""Renvoie l'index d'un noeud donné dans l'attribut "enfants" de la table parent_table."""
function child_index(parent_table::AbstractParentTable, node::AbstractNode)
    try
        findfirst(x -> name(x) == name(node), enfants(parent_table))
    catch e
        error("Aucune occurrence de ce noeud n'a été trouvée dans la table : ", node)
    end
end

"""Renvoie le parent d'un noeud donné selon la table parent_table."""
function parent(parent_table::AbstractParentTable, node::AbstractNode)
    parents(parent_table)[child_index(parent_table, node)]
end

"""Attribue un noeud parent au noeud donné."""
function set_parent!(parent_table::AbstractParentTable, node_index::Int64, parent::AbstractNode)
        parents(parent_table)[node_index] = parent
        parent_table
end

"""Si le nouveau poids "new_weight" est inférieur au poids actuel du noeud "node",
alors cette méthode lui attribue ce nouveau poids et lui affecte le noeud "parent"
comme parent dans la table "parent_table".
Sinon, rien ne se passe.
"""
function update_parent_and_weight!(parent_table::AbstractParentTable, node::AbstractNode, parent::AbstractNode, new_weight::Float64)
    if new_weight < min_weight(node)
        set_min_weight!(node, new_weight)
        node_index = child_index(parent_table, node)
        set_parent!(parent_table, node_index, parent)
    end
end

"""Renvoie le rang associé à un noeud donné selon la table parent_table."""
function rank(parent_table::AbstractParentTable, node_index::Int64)
    rank = ranks(parent_table)[node_index]
end

"""Incrémente de 1 le rang d'un noeud donné."""
function increase_rank!(parent_table::AbstractParentTable, node_index::Int64)
    ranks(parent_table)[node_index] += 1
end

"""Renvoie la racine d'un noeud donné selon la table parent_table, c'est-à-dire
le plus ancien parent de ce noeud. Au passage, compresse le chemin vers la racine.
"""
function root!(parent_table::AbstractParentTable, node_index::Int64)
    parent = parents(parent_table)[node_index]
    if name(enfants(parent_table)[node_index]) !== name(parent)
        parent = root!(parent_table, child_index(parent_table, parent))
        set_parent!(parent_table, node_index, parent)
    end
    parent
end

"""Réunit deux composantes connexes en joignant leurs deux racines si elles sont distinctes."""
function unite!(parent_table::AbstractParentTable, node_index1::Int64, node_index2::Int64)
    root1 = root!(parent_table, node_index1)
    root2 = root!(parent_table, node_index2)
    r1 = child_index(parent_table, root1)
    r2 = child_index(parent_table, root2)
    if rank(parent_table, r1) < rank(parent_table, r2)
        set_parent!(parent_table, r2, root1)
    elseif rank(parent_table, r1) > rank(parent_table, r2)
        set_parent!(parent_table, r1, root2)
    else
        set_parent!(parent_table, r2, root1)
        increase_rank!(parent_table, r1)
    end
end

"""Construit un objet ParentTable aux dimensions d'un graphe donné.
Kruskal : initialement, chaque noeud est son propre parent.
"""
function init_parent_table_kruskal(graph::Graph{T}) where T
    graph_nodes_copy = copy(nodes(graph))
    parent_table = ParentTable{T}(nodes(graph), graph_nodes_copy, zeros(length(nodes(graph))))
end

"""Construit un objet ParentTable aux dimensions d'un graphe donné.
Prim : initialement, le parent de chaque noeud est 'nothing'.
"""
function init_parent_table_prim(graph::Graph{T}) where T
    graph_nodes_copy = copy(nodes(graph))
    parent_table = ParentTable{T}(nodes(graph), graph_nodes_copy, Vector{Node}([]))
end

"""Construit un objet ParentTable aux dimensions d'un graphe donné.
RSL : initialement, le parent de chaque noeud est un noeud appelé 'init'.
"""
function init_parent_table_RSL(graph::Graph{T}) where T
    graph_nodes_copy = copy(nodes(graph))
    parent_table = ParentTable{T}(nodes(graph), fill(Node{T}("init", [], Inf), length(nodes(graph))), Vector{Node}([]))
end
