using Test

include(joinpath(@__DIR__, "..", "phase1", "graph.jl"))
include(joinpath(@__DIR__, "..", "phase1", "main1.jl"))
include(joinpath(@__DIR__, "main2.jl"))

# Construction des graphes

# Graphe à 3 noeuds
node1 = Node{Int64}("1", 1)
node2 = Node{Int64}("2", 2)
node3 = Node{Int64}("3", 3)
edge1 = Edge{Int64}("1", node1, node2, 1)
edge2 = Edge{Int64}("2", node1, node3, 2)
edge3 = Edge{Int64}("3", node2, node3, 3)
graph_3n = Graph{Int64}("graph_3n",[],[])
add_edge!(graph_3n, edge1)
add_edge!(graph_3n, edge2)
add_edge!(graph_3n, edge3)
add_node!(graph_3n, node1)
add_node!(graph_3n, node2)
add_node!(graph_3n, node3)
parent_table_3n = init_parent_table_kruskal_kruskal(graph_3n)


# Graphe à 4 noeuds
graph_4n = Graph{Int64}("graph_4n", [Node{Int64}("1", 1), Node{Int64}("2", 2), Node{Int64}("3", 3), Node{Int64}("4", 4)],
[Edge{Int64}("1", Node{Int64}("1", 1), Node{Int64}("2", 2), 6),
Edge{Int64}("2", Node{Int64}("1", 1), Node{Int64}("3", 3), 5),
Edge{Int64}("3", Node{Int64}("1", 1), Node{Int64}("4", 4), 4),
Edge{Int64}("4", Node{Int64}("2", 2), Node{Int64}("3", 3), 3),
Edge{Int64}("5", Node{Int64}("2", 2), Node{Int64}("4", 4), 2),
Edge{Int64}("6", Node{Int64}("3", 3), Node{Int64}("4", 4), 1)])
parent_table_4n = init_parent_table_kruskal_kruskal(graph_4n)


# Graphe à 6 noeuds
graph_6n = Graph{Int64}("graph_6n", [Node{Int64}("1", 1), Node{Int64}("2", 2), Node{Int64}("3", 3), Node{Int64}("4", 4), Node{Int64}("5", 5), Node{Int64}("6", 6)],
[Edge{Int64}("", Node{Int64}("1", 1), Node{Int64}("2", 2), 10),
Edge{Int64}("", Node{Int64}("1", 1), Node{Int64}("3", 3), 9),
Edge{Int64}("", Node{Int64}("1", 1), Node{Int64}("4", 4), 8),
Edge{Int64}("", Node{Int64}("2", 2), Node{Int64}("3", 3), 7),
Edge{Int64}("", Node{Int64}("3", 3), Node{Int64}("4", 4), 6),
Edge{Int64}("", Node{Int64}("4", 4), Node{Int64}("5", 5), 5),
Edge{Int64}("", Node{Int64}("2", 2), Node{Int64}("5", 5), 4),
Edge{Int64}("", Node{Int64}("3", 3), Node{Int64}("6", 6), 3)])
parent_table_6n = init_parent_table_kruskal_kruskal(graph_6n)


# Graphe à 10 noeuds
graph_10n = Graph{Int64}("graph10n", [], [])
for i = 1 : 10
    add_node!(graph_10n, Node{Int64}(string(i), i))
end
for e = 1 : 9
    add_edge!(graph_10n, Edge{Int64}("", nodes(graph_10n)[e], nodes(graph_10n)[e + 1], 100 - e))
end
parent_table_10n = init_parent_table_kruskal_kruskal(graph_10n)


# Graphe à n noeuds : bayg29.tsp
graph_path1 = joinpath(@__DIR__, "..", "..", "instances/stsp/bayg29.tsp")
graph_Nn1 = main1(graph_path1)
parent_table_Nn1 = init_parent_table_kruskal_kruskal(graph_Nn1)

# Graphe à n noeuds : bays29.tsp
graph_path2 = joinpath(@__DIR__, "..", "..", "instances/stsp/bays29.tsp")
graph_Nn2 = main1(graph_path2)
parent_table_Nn2 = init_parent_table_kruskal_kruskal(graph_Nn2)

# Graphe à n noeuds : dantzig42.tsp
graph_path3 = joinpath(@__DIR__, "..", "..", "instances/stsp/dantzig42.tsp")
graph_Nn3 = main1(graph_path3)
parent_table_Nn3 = init_parent_table_kruskal_kruskal(graph_Nn3)


# Fonctions

"""Renvoie un booléen. Indique si la méthode parent() renvoie bien le parent associé au noeud donné dans la table parent_table"""
function is_parent_ok(parent_table::AbstractParentTable)
    parent_ok = true
    for i = 1 : length(parents(parent_table))
        if parents(parent_table)[i] != parent(parent_table, enfants(parent_table)[i])
            parent_ok = false
        end
    end
    parent_ok
end

"""Renvoie un booléen. Indique si tous les noeuds d'un graph donné ont la même racine."""
function is_root_unique(parent_table::AbstractParentTable)
    nodes = enfants(parent_table)
    for i = 1 : length(nodes)
        if root(parent_table, nodes[i]) != root(parent_table, nodes[1])
            return false
        end
    end
    return true
end

# Tests

# Au départ, chaque noeud est son propre parent
@test enfants(parent_table_3n) == parents(parent_table_3n)
@test enfants(parent_table_4n) == parents(parent_table_4n)
@test enfants(parent_table_6n) == parents(parent_table_6n)
@test enfants(parent_table_10n) == parents(parent_table_10n)
@test enfants(parent_table_Nn1) == parents(parent_table_Nn1)
@test enfants(parent_table_Nn2) == parents(parent_table_Nn2)
@test enfants(parent_table_Nn3) == parents(parent_table_Nn3)

# Test de la méthode parent()
@test is_parent_ok(parent_table_3n)
@test is_parent_ok(parent_table_4n)
@test is_parent_ok(parent_table_6n)
@test is_parent_ok(parent_table_10n)
@test is_parent_ok(parent_table_Nn1)
@test is_parent_ok(parent_table_Nn2)
@test is_parent_ok(parent_table_Nn3)


# Appliquons l'algorithme de Kruskal aux graphes
min_tree_3n = main2(graph_3n)
min_tree_4n = main2(graph_4n)
min_tree_6n = main2(graph_6n)
min_tree_10n = main2(graph_10n)
min_tree_Nn1 = main2(graph_Nn1)
min_tree_Nn2 = main2(graph_Nn2)
min_tree_Nn3 = main2(graph_Nn3)


# L'arbre de recouvrement minimal devrait avoir exactmeent une arête de moins que de sommets
@test length(edges(min_tree_3n)) == length(nodes(min_tree_3n)) - 1
@test length(edges(min_tree_4n)) == length(nodes(min_tree_4n)) - 1
@test length(edges(min_tree_6n)) == length(nodes(min_tree_6n)) - 1
@test length(edges(min_tree_10n)) == length(nodes(min_tree_10n)) - 1
@test length(edges(min_tree_Nn1)) == length(nodes(min_tree_Nn1)) - 1
@test length(edges(min_tree_Nn2)) == length(nodes(min_tree_Nn2)) - 1
@test length(edges(min_tree_Nn3)) == length(nodes(min_tree_Nn3)) - 1
