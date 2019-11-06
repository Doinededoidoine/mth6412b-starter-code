

include(joinpath(@__DIR__, "..", "phase1", "graph.jl"))
include(joinpath(@__DIR__, "..", "phase2", "parent_table.jl"))

"""Renvoie un arbre de recouvrement minimal du graphe symétrique en entrée en
utilisant l'algorithme de Kruskal. La méthode renvoie un objet de type Graph.
ATTENTION : dans le graphe en entrée, les noeuds doivent tous avoir un attribut "name" différent.
"""
function main2(graph::AbstractGraph{T}) where T
    parent_table = init_parent_table_kruskal(graph)
    # on trie les arêtes par poids croissant.
    edges_copy = copy(edges(graph))
    sort!(edges_copy, by=weight)
    # On initialise l'arbre de recouvrement minimum en lui ajoutant tous les sommets du graphe.
    min_tree = Graph{T}("min_tree", [], [])
    for node in nodes(graph)
        add_node!(min_tree, node)
    end
    # On parcourt les arêtes du graphe.
    for edge in edges_copy
        # On récupère les extrémités de l'arête.
        node1 = s_node(edge)
        node2 = d_node(edge)
        # Si les deux noeuds ont des racines différentes, on ajoute
        # l'arête à l'arbre de recouvrement minimum puis on réunit les
        # composantes connexes en désignant l'une des deux racines comme
        # le nouveau parent de l'autre racine.
        if root!(parent_table, node1) != root!(parent_table, node2)
            add_edge!(min_tree, edge)
            unite!(parent_table, node1, node2)
        end
    end
    # On affiche puis on renvoie l'arbre de recouvrement minimum.
    show(min_tree)
    min_tree
end
