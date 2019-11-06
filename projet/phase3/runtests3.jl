using Test

include(joinpath(@__DIR__, "..", "phase1", "graph.jl"))
include(joinpath(@__DIR__, "..", "phase1", "main1.jl"))
include(joinpath(@__DIR__, "main3.jl"))

# Construction des graphes

# Graphe à 3 noeuds
node1 = Node{Float64}("1", 1.0, Inf)
node2 = Node{Float64}("2", 2.0, Inf)
node3 = Node{Float64}("3", 3.0, Inf)
edge1 = Edge{Float64}("1", node1, node2, 1)
edge2 = Edge{Float64}("2", node1, node3, 2)
edge3 = Edge{Float64}("3", node2, node3, 3)
graph_3n = Graph{Float64}("graph_3n",[],[])
add_edge!(graph_3n, edge1)
add_edge!(graph_3n, edge2)
add_edge!(graph_3n, edge3)
add_node!(graph_3n, node1)
add_node!(graph_3n, node2)
add_node!(graph_3n, node3)
parent_table_kruskal_3n = init_parent_table_kruskal(graph_3n)
parent_table_prim_3n = init_parent_table_prim(graph_3n)


# Graphe à 4 noeuds
graph_4n = Graph{Float64}("graph_4n", [Node{Float64}("1", 1.0, Inf), Node{Float64}("2", 2.0, Inf), Node{Float64}("3", 3.0, Inf), Node{Float64}("4", 4.0, Inf)],
[Edge{Float64}("1", Node{Float64}("1", 1.0, Inf), Node{Float64}("2", 2.0, Inf), 6),
Edge{Float64}("2", Node{Float64}("1", 1.0, Inf), Node{Float64}("3", 3.0, Inf), 5),
Edge{Float64}("3", Node{Float64}("1", 1.0, Inf), Node{Float64}("4", 4.0, Inf), 4),
Edge{Float64}("4", Node{Float64}("2", 2.0, Inf), Node{Float64}("3", 3.0, Inf), 3),
Edge{Float64}("5", Node{Float64}("2", 2.0, Inf), Node{Float64}("4", 4.0, Inf), 2),
Edge{Float64}("6", Node{Float64}("3", 3.0, Inf), Node{Float64}("4", 4.0, Inf), 1)])
parent_table_kruskal_4n = init_parent_table_kruskal(graph_4n)
parent_table_prim_4n = init_parent_table_prim(graph_4n)


# Graphe à 6 noeuds
graph_6n = Graph{Float64}("graph_6n", [Node{Float64}("1", 1.0, Inf), Node{Float64}("2", 2.0, Inf), Node{Float64}("3", 3.0, Inf), Node{Float64}("4", 4.0, Inf), Node{Float64}("5", 5.0, Inf), Node{Float64}("6", 6.0, Inf)],
[Edge{Float64}("", Node{Float64}("1", 1.0, Inf), Node{Float64}("2", 2.0, Inf), 10),
Edge{Float64}("", Node{Float64}("1", 1.0, Inf), Node{Float64}("3", 3.0, Inf), 9),
Edge{Float64}("", Node{Float64}("1", 1.0, Inf), Node{Float64}("4", 4.0, Inf), 8),
Edge{Float64}("", Node{Float64}("2", 2.0, Inf), Node{Float64}("3", 3.0, Inf), 7),
Edge{Float64}("", Node{Float64}("3", 3.0, Inf), Node{Float64}("4", 4.0, Inf), 6),
Edge{Float64}("", Node{Float64}("4", 4.0, Inf), Node{Float64}("5", 5.0, Inf), 5),
Edge{Float64}("", Node{Float64}("2", 2.0, Inf), Node{Float64}("5", 5.0, Inf), 4),
Edge{Float64}("", Node{Float64}("3", 3.0, Inf), Node{Float64}("6", 6.0, Inf), 3)])
parent_table_kruskal_6n = init_parent_table_kruskal(graph_6n)
parent_table_prim_6n = init_parent_table_prim(graph_6n)


# Graphe à 9 noeuds : exemple du cours
graph_9n = Graph{Float64}("graph_9n", [], [])
for i = 1 : 9
    node = Node{Float64}(string(i), i, Inf)
    add_node!(graph_9n, node)
end
nodes_9n = nodes(graph_9n)
add_edge!(graph_9n, Edge{Float64}("", nodes_9n[1], nodes_9n[2], 4))
add_edge!(graph_9n, Edge{Float64}("", nodes_9n[1], nodes_9n[8], 8))
add_edge!(graph_9n, Edge{Float64}("", nodes_9n[2], nodes_9n[3], 8))
add_edge!(graph_9n, Edge{Float64}("", nodes_9n[2], nodes_9n[8], 11))
add_edge!(graph_9n, Edge{Float64}("", nodes_9n[3], nodes_9n[4], 7))
add_edge!(graph_9n, Edge{Float64}("", nodes_9n[3], nodes_9n[6], 4))
add_edge!(graph_9n, Edge{Float64}("", nodes_9n[3], nodes_9n[9], 2))
add_edge!(graph_9n, Edge{Float64}("", nodes_9n[4], nodes_9n[5], 9))
add_edge!(graph_9n, Edge{Float64}("", nodes_9n[4], nodes_9n[6], 14))
add_edge!(graph_9n, Edge{Float64}("", nodes_9n[5], nodes_9n[6], 10))
add_edge!(graph_9n, Edge{Float64}("", nodes_9n[6], nodes_9n[7], 2))
add_edge!(graph_9n, Edge{Float64}("", nodes_9n[7], nodes_9n[8], 1))
add_edge!(graph_9n, Edge{Float64}("", nodes_9n[7], nodes_9n[9], 6))
add_edge!(graph_9n, Edge{Float64}("", nodes_9n[8], nodes_9n[9], 7))
parent_table_kruskal_9n = init_parent_table_kruskal(graph_9n)
parent_table_prim_9n = init_parent_table_prim(graph_9n)


# Graphe à 10 noeuds
graph_10n = Graph{Float64}("graph10n", [], [])
for i = 1 : 10
    add_node!(graph_10n, Node{Float64}(string(i), i, Inf))
end
for e = 1 : 9
    add_edge!(graph_10n, Edge{Float64}("", nodes(graph_10n)[e], nodes(graph_10n)[e + 1], 100 - e))
end
parent_table_kruskal_10n = init_parent_table_kruskal(graph_10n)
parent_table_prim_10n = init_parent_table_prim(graph_10n)


# Graphe à n noeuds : bayg29.tsp
graph_path1 = joinpath(@__DIR__, "..", "..", "instances/stsp/bayg29.tsp")
graph_Nn1 = main1(graph_path1)
parent_table_kruskal_Nn1 = init_parent_table_kruskal(graph_Nn1)
parent_table_prim_Nn1 = init_parent_table_prim(graph_Nn1)

# Graphe à n noeuds : bays29.tsp
graph_path2 = joinpath(@__DIR__, "..", "..", "instances/stsp/bays29.tsp")
graph_Nn2 = main1(graph_path2)
parent_table_kruskal_Nn2 = init_parent_table_kruskal(graph_Nn2)
parent_table_prim_Nn2 = init_parent_table_prim(graph_Nn2)

# Graphe à n noeuds : dantzig42.tsp
graph_path3 = joinpath(@__DIR__, "..", "..", "instances/stsp/dantzig42.tsp")
graph_Nn3 = main1(graph_path3)
parent_table_kruskal_Nn3 = init_parent_table_kruskal(graph_Nn3)
parent_table_prim_Nn3 = init_parent_table_prim(graph_Nn3)

# Graphe à n noeuds sans coordonnées : fri26.tsp
graph_path4 = joinpath(@__DIR__, "..", "..", "instances/stsp/fri26.tsp")
graph_Nn4 = main1(graph_path4)
parent_table_kruskal_Nn4 = init_parent_table_kruskal(graph_Nn4)
parent_table_prim_Nn4 = init_parent_table_prim(graph_Nn4)


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
        if root!(parent_table, nodes[i]) != root!(parent_table, nodes[1])
            return false
        end
    end
    return true
end

# Tests

# Test de la méthode parent()
@test is_parent_ok(parent_table_kruskal_3n)
@test is_parent_ok(parent_table_kruskal_4n)
@test is_parent_ok(parent_table_kruskal_6n)
@test is_parent_ok(parent_table_kruskal_9n)
@test is_parent_ok(parent_table_kruskal_10n)
@test is_parent_ok(parent_table_kruskal_Nn1)
@test is_parent_ok(parent_table_kruskal_Nn2)
@test is_parent_ok(parent_table_kruskal_Nn3)
@test is_parent_ok(parent_table_kruskal_Nn4)
@test is_parent_ok(parent_table_prim_3n)
@test is_parent_ok(parent_table_prim_4n)
@test is_parent_ok(parent_table_prim_6n)
@test is_parent_ok(parent_table_prim_9n)
@test is_parent_ok(parent_table_prim_10n)
@test is_parent_ok(parent_table_prim_Nn1)
@test is_parent_ok(parent_table_prim_Nn2)
@test is_parent_ok(parent_table_prim_Nn3)
@test is_parent_ok(parent_table_prim_Nn4)


# Appliquons l'algorithme de Kruskal aux graphes
min_tree_kruskal_3n = kruskal(graph_3n)
min_tree_kruskal_4n = kruskal(graph_4n)
min_tree_kruskal_6n = kruskal(graph_6n)
min_tree_kruskal_9n = kruskal(graph_9n)
min_tree_kruskal_10n = kruskal(graph_10n)
min_tree_kruskal_Nn1 = kruskal(graph_Nn1)
min_tree_kruskal_Nn2 = kruskal(graph_Nn2)
min_tree_kruskal_Nn3 = kruskal(graph_Nn3)
min_tree_kruskal_Nn4 = kruskal(graph_Nn4)

# Appliquons l'algorithme de Prim aux graphes (le noeud de départ est choisi arbitrairement)
node0_3n = nodes(graph_3n)[2]
min_tree_prim_3n = prim(graph_3n, node0_3n)
node0_4n = nodes(graph_4n)[2]
min_tree_prim_4n = prim(graph_4n, node0_4n)
node0_6n = nodes(graph_6n)[2]
min_tree_prim_6n = prim(graph_6n, node0_6n)
node0_9n = nodes(graph_9n)[2]
min_tree_prim_9n = prim(graph_9n, node0_9n)
node0_10n = nodes(graph_10n)[2]
min_tree_prim_10n = prim(graph_10n, node0_10n)
node0_Nn1 = nodes(graph_Nn1)[2]
min_tree_prim_Nn1 = prim(graph_Nn1, node0_Nn1)
node0_Nn2 = nodes(graph_Nn2)[2]
min_tree_prim_Nn2 = prim(graph_Nn2, node0_Nn2)
node0_Nn3 = nodes(graph_Nn3)[2]
min_tree_prim_Nn3 = prim(graph_Nn3, node0_Nn3)
node0_Nn4 = nodes(graph_Nn4)[2]
min_tree_prim_Nn4 = prim(graph_Nn4, node0_Nn4)


# L'arbre de recouvrement minimal devrait avoir exactement une arête de moins que de sommets
@test length(edges(min_tree_kruskal_3n)) == length(nodes(min_tree_kruskal_3n)) - 1
@test length(edges(min_tree_kruskal_4n)) == length(nodes(min_tree_kruskal_4n)) - 1
@test length(edges(min_tree_kruskal_6n)) == length(nodes(min_tree_kruskal_6n)) - 1
@test length(edges(min_tree_kruskal_9n)) == length(nodes(min_tree_kruskal_9n)) - 1
@test length(edges(min_tree_kruskal_10n)) == length(nodes(min_tree_kruskal_10n)) - 1
@test length(edges(min_tree_kruskal_Nn1)) == length(nodes(min_tree_kruskal_Nn1)) - 1
@test length(edges(min_tree_kruskal_Nn2)) == length(nodes(min_tree_kruskal_Nn2)) - 1
@test length(edges(min_tree_kruskal_Nn3)) == length(nodes(min_tree_kruskal_Nn3)) - 1
@test length(edges(min_tree_kruskal_Nn4)) == length(nodes(min_tree_kruskal_Nn4)) - 1
@test length(edges(min_tree_prim_3n)) == length(nodes(min_tree_prim_3n)) - 1
@test length(edges(min_tree_prim_4n)) == length(nodes(min_tree_prim_4n)) - 1
@test length(edges(min_tree_prim_6n)) == length(nodes(min_tree_prim_6n)) - 1
@test length(edges(min_tree_prim_9n)) == length(nodes(min_tree_prim_9n)) - 1
@test length(edges(min_tree_prim_10n)) == length(nodes(min_tree_prim_10n)) - 1
@test length(edges(min_tree_prim_Nn1)) == length(nodes(min_tree_prim_Nn1)) - 1
@test length(edges(min_tree_prim_Nn2)) == length(nodes(min_tree_prim_Nn2)) - 1
@test length(edges(min_tree_prim_Nn3)) == length(nodes(min_tree_prim_Nn3)) - 1
@test length(edges(min_tree_prim_Nn4)) == length(nodes(min_tree_prim_Nn4)) - 1
