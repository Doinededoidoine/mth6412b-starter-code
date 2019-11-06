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
ranks(parent_table::AbstractParentTable, node::AbstractNode) = parent_table.ranks

"""Renvoie le parent d'un noeud donné selon la table parent_table."""
function parent(parent_table::AbstractParentTable, node::AbstractNode)
    for (i, enfant) in enumerate(enfants(parent_table))
        if name(enfant) == name(node)
            parent = parents(parent_table)[i]
            return parent
        end
    end
    error("Aucun parent n'a été trouvé pour ce noeud : ", node)
end

"""Attribue un noeud parent au noeud donné."""
function set_parent!(parent_table::AbstractParentTable, node::AbstractNode, parent::AbstractNode)
    for (i, enfant) in enumerate(enfants(parent_table))
        if name(enfant) == name(node)
            parents(parent_table)[i] = parent
            return parent_table
        end
    end
    error("Ce noeud n'existe pas : ", node)
end

"""Si le nouveau poids "new_weight" est inférieur au poids actuel du noeud "node",
alors cette méthode lui attribue ce nouveau poids et lui affecte le noeud "parent"
comme parent dans la table "parent_table".
Sinon, rien ne se passe.
"""
function update_parent_and_weight!(parent_table::AbstractParentTable, node::AbstractNode, parent::AbstractNode, new_weight::Float64)
    if new_weight < min_weight(node)
        set_min_weight!(node, new_weight)
        set_parent!(parent_table, node, parent)
    end
end

"""Renvoie le rang associé à un noeud donné selon la table parent_table."""
function rank(parent_table::AbstractParentTable, node::AbstractNode)
    for (i, enfant) in enumerate(enfants(parent_table))
        if name(enfant) == name(node)
            rank = ranks(parent_table, node)[i]
            return rank
        end
    end
    error("Aucun rang n'a été trouvé pour ce noeud : ", node)
end

"""Attribue un rang au noeud donné."""
function increase_rank!(parent_table::AbstractParentTable, node::AbstractNode)
    for (i, enfant) in enumerate(enfants(parent_table))
        if name(enfant) == name(node)
            ranks(parent_table, node)[i] += 1
            return parent_table
        end
    end
    error("Ce noeud n'existe pas : ", node)
end

"""Renvoie la racine d'un noeud donné selon la table parent_table, c'est-à-dire
le plus ancien parent de ce noeud. Au passage, compresse le chemin vers la racine.
"""
function root!(parent_table::AbstractParentTable, node::AbstractNode)
    parent_ = parent(parent_table, node)
    if name(node) !== name(parent_)
        parent_ = root!(parent_table, parent_)
        set_parent!(parent_table, node, parent_)
    end
    parent_
end

"""Réunit deux composantes connexes en joignant leurs deux racines si elles sont distinctes."""
function unite!(parent_table::AbstractParentTable, node1::AbstractNode, node2::AbstractNode)
    root1 = root!(parent_table, node1)
    root2 = root!(parent_table, node2)
    if rank(parent_table, root1) < rank(parent_table, root2)
        set_parent!(parent_table, root2, root1)
    elseif rank(parent_table, root1) > rank(parent_table, root2)
        set_parent!(parent_table, root1, root2)
    else
        set_parent!(parent_table, root2, root1)
        increase_rank!(parent_table, root1)
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
