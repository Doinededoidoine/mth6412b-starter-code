include(joinpath(@__DIR__, "..", "phase1", "main1.jl"))
include(joinpath(@__DIR__, "..", "phase2", "parent_table.jl"))
include(joinpath(@__DIR__, "..", "phase3", "priority_queue.jl"))
include(joinpath(@__DIR__, "..", "phase3", "main3.jl"))
include(joinpath(@__DIR__, "..", "phase4", "degree_table.jl"))


"""Première partie de l'algorithme de Rosenkrantz, Stearns et Lewis. Par défaut,
le premier noeud du graphe est choisi comme point de départ. La méthode prend en
argument un arbre de recouvrement et renvoie un tableau de type ParentTable
contenant les relations de parenté entre les différents noeuds.
ATTENTION : dans le graphe en entrée, les noeuds doivent tous avoir un attribut "name" différent.
"""
function RSL_struct(graph::AbstractGraph{T}, starting_node::AbstractNode = nodes(graph)[1]) where T
    parent_table = init_parent_table_RSL(graph)
    q = PriorityQueue{Node{T}}([])
    node = Node{T}(name(starting_node), data(starting_node))
    node_temp = []
    set_parent!(parent_table, child_index(parent_table, node), node)
    push!(q, node)
    #l'algorithme s'arrête lorsque la pile q est vide ce qui indique que les parents de tous les noeuds ont été attribués
    while !(is_empty(q))
        node_temp = []
        #on parcout toutes les arêtes du graphe de recouvrement minimal
        for edge in edges(graph)
            i1 = 0
            #tous les noeuds adjacents à la variable node n'ayant pas déja un parent s'en font attribuer un et sont ajouté à 'q'
            if (name(s_node(edge)) == name(node)) && (name(parent(parent_table, d_node(edge))) == "init")
                i1 = child_index(parent_table, d_node(edge))
                set_parent!(parent_table, i1, node)
                push!(q, d_node(edge))
                node_temp = d_node(edge)

            elseif (name(d_node(edge)) == name(node)) && (name(parent(parent_table, s_node(edge))) == "init")
                i1 = child_index(parent_table, s_node(edge))
                set_parent!(parent_table, i1, node)
                push!(q, s_node(edge))
                node_temp = s_node(edge)
            end
        end

        #si aucun noeud adjacent n'ayant pas déja un parent n'est trouvé, la variable node prend la valeur du dernier noeud rentré dans 'q'
        if node_temp == []
            node = popfirst!(q)
        #si un noeud a été trouvé, la variable node prend la valeur de ce noeud
        else
            node = node_temp
        end
    end
    #un tableau contenant les liens de parenté des différents noeuds est retourné
    return(parent_table)
end

"""Deuxième partie de l'algorithme de Rosenkrantz, Stearns et Lewis.
Renvoie une tournée optimisée passant par tous les points du graphe symétrique
en entrée. Par défaut, le premier noeud du graphe est choisi comme point de départ.
La méthode prend en argument un arbre de recouvrement, le graphe complet auquel il
correspond, une table des liens de parenté entre les noeuds (calculée par défaut
par la fonction RSL_struct) et un noeud de départ optionnel. Elle renvoie un objet
de type Graph contenant cette tournée.
ATTENTION : dans le graphe en entrée, les noeuds doivent tous avoir un attribut "name" différent.
"""
function RSL(graph::AbstractGraph{T}, graph_base::AbstractGraph{T}, starting_node::AbstractNode = nodes(graph)[1], parent_table::AbstractParentTable{T} = RSL_struct(graph, starting_node)) where T
    RSL_tree = Graph{T}("RSL_tree", [], [])
    node = Node{T}(name(starting_node), data(starting_node))
    add_node!(RSL_tree, node)
    node_temp = []
    j = 0
    #l'algorithme s'arrête lorsqu'il y a une arête de moins que de noeuds
    while ((nb_edges(RSL_tree)) < (nb_nodes(graph_base)-1))
        node_temp = []
        weight_temp = Inf
        #on parcourt toutes les arêtes du graphe de recouvrement minimal
        for edge in edges(graph)
            #nous cherchons l'arête ayant le poids le plus faible reliant la variable 'node a un autre noeud n'étant pas déjà dans le graphe 'RSL_tree'. La variable 'node_temp' prend la valeur de cet autre noeud.
            if ((name(s_node(edge)) == name(node)) && (name(parent(parent_table, d_node(edge))) == name(node)) && !(find_node(RSL_tree, d_node(edge))))
                if weight(edge) < weight_temp
                    node_temp = d_node(edge)
                    weight_temp = weight(edge)
                end
            elseif ((name(d_node(edge)) == name(node)) && (name(parent(parent_table, s_node(edge))) == name(node)) && !(find_node(RSL_tree, s_node(edge))))
                if weight(edge) < weight_temp
                    node_temp = s_node(edge)
                    weight_temp = weight(edge)
                end
            end
        end
        #si la variable 'node_temp' resort de la loop avec une valeur nulle, nous sommes à une feuille. Nous remontons donc dans l'arbre.
        if (node_temp == [])
            j = j + 1
            node = nodes(RSL_tree)[end-j]
        #sinon, la variable 'node' prend la valeur de 'node_temp' et ce noeud est ajouté à l'arbre 'RSL_tree'
        else
            j = 0
            node = node_temp
            add_node!(RSL_tree, node)
            #Nous ajoutons une arête reliant les deux derniers noeuds ajouté dans le graphe 'RSL_tree'
            add_edge!(RSL_tree, Edge("", nodes(RSL_tree)[end], nodes(RSL_tree)[end - 1], findweight(graph_base, nodes(RSL_tree)[end], nodes(RSL_tree)[end - 1])))
        end
    end
    #Nous ajoutons l'arête finale reliant le premier noeud au dernier noeud au graphe
    add_edge!(RSL_tree, Edge("", nodes(RSL_tree)[end], nodes(RSL_tree)[1], findweight(graph_base, nodes(RSL_tree)[end], nodes(RSL_tree)[1])))
    return(RSL_tree)
end

"""Renvoie la valeur du coefficient d'un noeud selon la table n_pi en entrée."""
function coef(coefs::Vector{Float64}, node::AbstractNode, nodes::Vector{Node{T}}) where T
    i = findfirst(n -> name(n) == name(node), nodes)
    coefs[i]
end

"""Calcule et renvoie la borne inférieure "w" que l'algorithme de Held et Karp cherche à maximiser."""
function lower_bound(one_tree::AbstractGraph{T}, degree_table::AbstractDegreeTable{T}, n_pi::Vector{Float64}) where T
    return graphweight(one_tree) + sum((degrees(degree_table)[i] - 2) * n_pi[i] for i = 1 : length(n_pi))
end

"""Ajoute au poids de chaque arête le coefficient de chacune de ses extrémités."""
function modify_costs!(graph::AbstractGraph, n_pi::Vector{Float64})
    for edge in edges(graph)
        set_weight!(edge, weight(edge) + coef(n_pi, s_node(edge), nodes(graph)) + coef(n_pi, d_node(edge), nodes(graph)))
    end
end

"""Soustrait au poids de chaque arête le coefficient de chacune de ses extrémités."""
function restore_costs!(graph::AbstractGraph, n_pi::Vector{Float64})
    for edge in edges(graph)
        set_weight!(edge, weight(edge) - coef(n_pi, s_node(edge), nodes(graph)) - coef(n_pi, d_node(edge), nodes(graph)))
    end
end

"""Renvoie une tournée optimisée passant par tous les points du graphe symétrique
en entrée en utilisant l'algorithme de Held et Karp. La méthode renvoie un objet
de type Graph contenant la meilleure tournée obtenue après "max_iter" itérations.
ATTENTION : dans le graphe en entrée, les noeuds doivent tous avoir un attribut "name" différent.
"""
function HK(graph::AbstractGraph{T}, algo::String = "prim", starting_node::AbstractNode = nodes(graph)[1], step::Float64 = 2.0, max_iter::Int64 = 20, max_gap::Float64 = 0.0) where T
    # Initialisation : on résoud le problème avec RSL pour commencer avec une bonne référence.
    one_tree = algo == "kruskal" ? kruskal(graph) : prim(graph, starting_node)
    ref_tour = RSL(one_tree, graph)
    HK_tour = Graph{T}("HK_tour", [], [])                       # Tournée courante
    min_tour = ref_tour                                         # Meilleure tournée
    degree_table = init_degree_table(graph)                     # Degré de chaque noeud
    n_pi = Vector{Float64}(step * ones(nb_nodes(graph)))        # Coefficient pi de chaque noeud
    w = -Inf                                                    # Borne inférieure courante
    max_w = w                                                   # Meilleure borne inférieure
    max_n_pi = n_pi                                             # Meilleur jeu de coefficients
    iter = 0                                                    # Itérateur
    gap = Inf                                                   # Ecart entre la borne inférieure et le coût
    while iter < max_iter && gap > max_gap                      # du meilleur cycle obtenu par la méthode RSL
        # On applique les coûts modifiés au graphe, puis on détermine l'arbre de
        # coût modifié minimum grâce à Prim ou Kruskal. On lui applique ensuite
        # RSL afin d'obtenir une tournée.
        modify_costs!(graph, n_pi)
        one_tree = algo == "kruskal" ? kruskal(graph) : prim(graph, starting_node)
        HK_tour = RSL(one_tree, graph)
        # On ajoute à l'arbre de recouvrement l'arête qu'il lui manque pour être un 1-arbre
        # de coût modifié minimum, puis on calcule le degré de chaque noeud dans ce 1-arbre.
        new_edge = sort(filter(e -> (name(s_node(e)) == name(nodes(graph)[1]) || name(d_node(e)) == name(nodes(graph)[1])) && !(e in edges(one_tree)), edges(graph)), by=weight)[1]
        add_edge!(one_tree, Edge{T}(name(new_edge), s_node(new_edge), d_node(new_edge), weight(new_edge)))
        degree_table = d_table(one_tree)
        # On rend au graphe, au 1-arbre et à la tournée leurs poids d'origine.
        restore_costs!(graph, n_pi)
        restore_costs!(one_tree, n_pi)
        restore_costs!(HK_tour, n_pi)
        # On calcule une borne inférieure du coût d'une tournée. Si elle est meilleure que la
        # précédente meilleure valeur obtenue, on la remplace et on met à jour les paramètres.
        w = lower_bound(one_tree, degree_table, n_pi)
        if w > max_w
            max_w = w
            max_n_pi = n_pi
            min_tour = HK_tour
            gap = graphweight(ref_tour) - max_w
        end
        # On met à jour les coefficients pi des noeuds.
        for i = 1 : length(n_pi)
            n_pi[i] += step * (degrees(degree_table)[i] - 2)
        end
        iter += 1
    end
    return min_tour
end

"""Fonction pratique permettant d'appliquer l'algorithme RSL sur un graphe
de recouvrement minimal créé à l'aide de l'algorithme de Prim appliqué sur une
instance donné en argument. Le noeud utilisé comme racine est le premier noeud
du graphique de recouvrement minimal. Si true est donné en argument, la fonction
va retourner l'image du cycle résultant de l'algorithme.
"""
function RSL_total_prim(chemin::String, image::Bool = false)
    graph = main1(chemin)
    graph_min = prim(graph)
    RSL_tree = RSL(graph_min, graph)
    if image
        display(plot_graph(RSL_tree))
    end
    println(graphweight(RSL_tree))
    RSL_tree
end

"""Fonction pratique permettant d'appliquer l'algorithme RSL sur un graphe
de recouvrement minimal créé à l'aide de l'algorithme de Kruskal appliqué sur une
instance donné en argument. Le noeud utilisé comme racine est le premier noeud
du graphique de recouvrement minimal. Si true est donné en argument, la fonction
va retourner l'image du cycle résultant de l'algorithme.
"""
function RSL_total_kruskal(chemin::String, image::Bool = false)
     graph = main1(chemin)
     graph_min = kruskal(graph)
     RSL_tree = RSL(graph_min, graph)
     if image
         display(plot_graph(RSL_tree))
     end
     println(graphweight(RSL_tree))
     RSL_tree
end

"""Fonction pratique permettant d'appliquer l'algorithme HK sur un graphe
de recouvrement minimal créé à l'aide de l'algorithme de Prim appliqué sur une
instance donné en argument. Le noeud utilisé comme racine est le premier noeud
du graphique de recouvrement minimal. Si true est donné en argument, la fonction
va retourner l'image du cycle résultant de l'algorithme.
"""
function HK_total_prim(chemin::String, image::Bool = false)
    graph = main1(chemin)
    HK_tree = HK(graph, "prim")
    if image
        display(plot_graph(HK_tree))
    end
    println(graphweight(HK_tree))
    HK_tree
end

"""Fonction pratique permettant d'appliquer l'algorithme HK sur un graphe
de recouvrement minimal créé à l'aide de l'algorithme de Kruskal appliqué sur une
instance donné en argument. Le noeud utilisé comme racine est le premier noeud
du graphique de recouvrement minimal. Si true est donné en argument, la fonction
va retourner l'image du cycle résultant de l'algorithme.
"""
function HK_total_kruskal(chemin::String, image::Bool = false)
    graph = main1(chemin)
    HK_tree = HK(graph, "kruskal")
    if image
        display(plot_graph(HK_tree))
    end
    println(graphweight(HK_tree))
    HK_tree
end


"""Fonction pratique permettant de crér un tableau avec les poids des différentes
solutions trouvées à l'aide de l'algorithme RSL. Ce tableau compare les résultats
trouvés avec les paramètres kruskal-premier à ceux trouvés avec les paramètres
prim-premier.
"""
function resultats()
    resultat = zeros(14,3)
    resultat[1,2]=graphweight(RSL_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp"))
    resultat[2,2]=graphweight(RSL_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp"))
    resultat[3,2]=graphweight(RSL_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp"))
    resultat[4,2]=graphweight(RSL_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp"))
    resultat[5,2]=graphweight(RSL_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp"))
    resultat[6,2]=graphweight(RSL_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp"))
    resultat[7,2]=graphweight(RSL_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp"))
    resultat[8,2]=graphweight(RSL_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp"))
    resultat[9,2]=graphweight(RSL_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp"))
    resultat[10,2]=graphweight(RSL_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp"))
    resultat[11,2]=graphweight(RSL_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp"))
    resultat[12,2]=graphweight(RSL_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp"))
    resultat[13,2]=graphweight(RSL_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp"))
    resultat[14,2]=graphweight(RSL_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp"))
    resultat[1,1]=graphweight(RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp"))
    resultat[2,1]=graphweight(RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp"))
    resultat[3,1]=graphweight(RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp"))
    resultat[4,1]=graphweight(RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp"))
    resultat[5,1]=graphweight(RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp"))
    resultat[6,1]=graphweight(RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp"))
    resultat[7,1]=graphweight(RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp"))
    resultat[8,1]=graphweight(RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp"))
    resultat[9,1]=graphweight(RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp"))
    resultat[10,1]=graphweight(RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp"))
    resultat[11,1]=graphweight(RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp"))
    resultat[12,1]=graphweight(RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp"))
    resultat[13,1]=graphweight(RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp"))
    resultat[14,1]=graphweight(RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp"))
    i = 1
    for i in 1:14
        resultat[i,3] = resultat[i,1] - resultat[i,2]
    end
    return resultat
end

"""Fonction pratique permettant de crér un tableau avec les poids des différentes
solutions trouvées à l'aide de l'algorithme HK. Ce tableau compare les résultats
trouvés avec les paramètres Kruskal à ceux trouvés avec les paramètres
Prim.
"""
function resultats_HK()
    resultat = zeros(14,3)
    resultat[1,2]=graphweight(HK_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp"))
    resultat[2,2]=graphweight(HK_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp"))
    resultat[3,2]=graphweight(HK_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp"))
    resultat[4,2]=graphweight(HK_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp"))
    resultat[5,2]=graphweight(HK_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp"))
    resultat[6,2]=graphweight(HK_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp"))
    resultat[7,2]=graphweight(HK_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp"))
    resultat[8,2]=graphweight(HK_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp"))
    resultat[9,2]=graphweight(HK_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp"))
    resultat[10,2]=graphweight(HK_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp"))
    resultat[11,2]=graphweight(HK_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp"))
    resultat[12,2]=graphweight(HK_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp"))
    resultat[13,2]=graphweight(HK_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp"))
    resultat[14,2]=graphweight(HK_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp"))
    resultat[1,1]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp"))
    resultat[2,1]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp"))
    resultat[3,1]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp"))
    resultat[4,1]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp"))
    resultat[5,1]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp"))
    resultat[6,1]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp"))
    resultat[7,1]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp"))
    resultat[8,1]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp"))
    resultat[9,1]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp"))
    resultat[10,1]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp"))
    resultat[11,1]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp"))
    resultat[12,1]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp"))
    resultat[13,1]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp"))
    resultat[14,1]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp"))
    i = 1
    for i in 1:14
        resultat[i,3] = resultat[i,1] - resultat[i,2]
    end
    return resultat
end

"""Fonction pratique permettant d'appliquer l'algorithme RSL sur un graphe
de recouvrement minimal créé à l'aide de l'algorithme de Prim appliqué sur une
instance donné en argument. Le noeud utilisé comme racine est le dernier noeud
du graphique de recouvrement minimal. Si true est donné en argument, la fonction
va retourner l'image du cycle résultant de l'algorithme.
"""
function RSL_total_prim_last(chemin::String, image::Bool = false)
     graph = main1(chemin)
     graph_min = prim(graph)
     RSL_tree = RSL(graph_min, graph, nodes(graph_min)[end])
     if image
         display(plot_graph(RSL_tree))
     end
     println(graphweight(RSL_tree))
     RSL_tree
end

"""Fonction pratique permettant d'appliquer l'algorithme RSL sur un graphe
de recouvrement minimal créé à l'aide de l'algorithme de Kruskal appliqué sur une
instance donné en argument. Le noeud utilisé comme racine est le dernier noeud
du graphique de recouvrement minimal. Si true est donné en argument, la fonction
va retourner l'image du cycle résultant de l'algorithme.
"""
function RSL_total_kruskal_last(chemin::String, image::Bool = false)
     graph = main1(chemin)
     graph_min = kruskal(graph)
     RSL_tree = RSL(graph_min, graph, nodes(graph_min)[end])
     if image
         display(plot_graph(RSL_tree))
     end
     println(graphweight(RSL_tree))
     RSL_tree
end

"""Fonction pratique permettant d'appliquer l'algorithme HK sur un graphe
de recouvrement minimal créé à l'aide de l'algorithme de Kruskal appliqué sur une
instance donné en argument. Le noeud utilisé comme racine est le dernier noeud
du graphique de recouvrement minimal. Si true est donné en argument, la fonction
va retourner l'image du cycle résultant de l'algorithme.
"""
function HK_total_prim_last(chemin::String, image::Bool = false)
     graph = main1(chemin)
     graph_min = prim(graph)
     HK_tree = HK(graph, "prim", nodes(graph_min)[end])
     if image
         display(plot_graph(HK_tree))
     end
     println(graphweight(HK_tree))
     HK_tree
end

"""Fonction pratique permettant d'appliquer l'algorithme HK sur un graphe
de recouvrement minimal créé à l'aide de l'algorithme de Kruskal appliqué sur une
instance donné en argument. Le noeud utilisé comme racine est le dernier noeud
du graphique de recouvrement minimal. Si true est donné en argument, la fonction
va retourner l'image du cycle résultant de l'algorithme.
"""
function HK_total_kruskal_last(chemin::String, image::Bool = false)
     graph = main1(chemin)
     graph_min = prim(graph)
     HK_tree = HK(graph, "kruskal", nodes(graph_min)[end])
     if image
         display(plot_graph(HK_tree))
     end
     println(graphweight(HK_tree))
     HK_tree
end

"""Fonction pratique permettant de créer un tableau avec les poids des différentes
solutions trouvées à l'aide de l'algorithme RSL. Ce tableau compare les résultats
trouvés avec les paramètres prim-dernier à ceux trouvés avec les paramètres
prim-premier.
"""
function resultats_prim_last_vs_first()
    resultat = zeros(14,3)
    resultat[1,1]=RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp")
    resultat[2,1]=RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp")
    resultat[3,1]=RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp")
    resultat[4,1]=RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp")
    resultat[5,1]=RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp")
    resultat[6,1]=RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp")
    resultat[7,1]=RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp")
    resultat[8,1]=RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp")
    resultat[9,1]=RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp")
    resultat[10,1]=RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp")
    resultat[11,1]=RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp")
    resultat[12,1]=RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp")
    resultat[13,1]=RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp")
    resultat[14,1]=RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp")
    resultat[1,2]=RSL_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp")
    resultat[2,2]=RSL_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp")
    resultat[3,2]=RSL_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp")
    resultat[4,2]=RSL_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp")
    resultat[5,2]=RSL_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp")
    resultat[6,2]=RSL_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp")
    resultat[7,2]=RSL_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp")
    resultat[8,2]=RSL_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp")
    resultat[9,2]=RSL_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp")
    resultat[10,2]=RSL_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp")
    resultat[11,2]=RSL_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp")
    resultat[12,2]=RSL_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp")
    resultat[13,2]=RSL_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp")
    resultat[14,2]=RSL_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp")
    i = 1
    for i in 1:14
        resultat[i,3] = resultat[i,1] - resultat[i,2]
    end
    return resultat
end

"""Fonction pratique permettant d'appliquer l'algorithme RSL sur un graphe
de recouvrement minimal créé à l'aide de l'algorithme de Prim appliqué sur une
instance donné en argument. Le noeud utilisé comme racine est le dernier noeud
du graphique de recouvrement minimal. Si true est donné en argument, la fonction
va retourner l'image du cycle résultant de l'algorithme.
"""
function RSL_total_prim_lightest(chemin::String, image::Bool = false)
     graph = main1(chemin)
     graph_min = prim(graph)
     RSL_tree = RSL(graph_min, graph, find_lightest_node(graph_min))
     if image
         display(plot_graph(RSL_tree))
     end
     println(graphweight(RSL_tree))
     RSL_tree
end

"""Fonction pratique permettant d'appliquer l'algorithme RSL sur un graphe
de recouvrement minimal créé à l'aide de l'algorithme de Prim appliqué sur une
instance donné en argument. Le noeud utilisé comme racine est un noeud appartenant
à l'arête la plus lourde du graphe de recouvrement minimal. Si true est donné en argument,
la fonction va retourner l'image du cycle résultant de l'algorithme.
"""
function RSL_total_prim_heaviest(chemin::String, image::Bool = false)
     graph = main1(chemin)
     graph_min = prim(graph)
     RSL_tree = RSL(graph_min, graph, find_heaviest_node(graph_min))
     if image
         display(plot_graph(RSL_tree))
     end
     println(graphweight(RSL_tree))
     RSL_tree
end

"""Fonction pratique permettant d'appliquer l'algorithme RSL sur un graphe
de recouvrement minimal créé à l'aide de l'algorithme de Kruskal appliqué sur une
instance donné en argument. Le noeud utilisé comme racine est un noeud appartenant
à l'arête la moins lourde du graphe de recouvrement minimal. Si true est donné en argument,
la fonction va retourner l'image du cycle résultant de l'algorithme.
"""
function RSL_total_kruskal_lightest(chemin::String, image::Bool = false)
     graph = main1(chemin)
     graph_min = kruskal(graph)
     RSL_tree = RSL(graph_min, graph, find_lightest_node(graph_min))
     if image
         display(plot_graph(RSL_tree))
     end
     println(graphweight(RSL_tree))
     RSL_tree
end

"""Fonction pratique permettant d'appliquer l'algorithme RSL sur un graphe
de recouvrement minimal créé à l'aide de l'algorithme de Kruskal appliqué sur une
instance donné en argument. Le noeud utilisé comme racine est un noeud appartenant
à l'arête la plus lourde du graphe de recouvrement minimal. Si true est donné en argument,
la fonction va retourner l'image du cycle résultant de l'algorithme.
"""
function RSL_total_kruskal_heaviest(chemin::String, image::Bool = false)
     graph = main1(chemin)
     graph_min = kruskal(graph)
     RSL_tree = RSL(graph_min, graph, find_heaviest_node(graph_min))
     if image
         display(plot_graph(RSL_tree))
     end
     println(graphweight(RSL_tree))
     RSL_tree
end

"""Fonction pratique permettant d'appliquer l'algorithme HK sur un graphe
de recouvrement minimal créé à l'aide de l'algorithme de Prim appliqué sur une
instance donné en argument. Le noeud utilisé comme racine est le dernier noeud
du graphique de recouvrement minimal. Si true est donné en argument, la fonction
va retourner l'image du cycle résultant de l'algorithme.
"""
function HK_total_prim_lightest(chemin::String, image::Bool = false)
     graph = main1(chemin)
     graph_min = prim(graph)
     HK_tree = HK(graph, "prim", find_lightest_node(graph_min))
     if image
         display(plot_graph(HK_tree))
     end
     println(graphweight(HK_tree))
     HK_tree
end

"""Fonction pratique permettant d'appliquer l'algorithme HK sur un graphe
de recouvrement minimal créé à l'aide de l'algorithme de Prim appliqué sur une
instance donné en argument. Le noeud utilisé comme racine est un noeud appartenant
à l'arête la plus lourde du graphe de recouvrement minimal. Si true est donné en argument,
la fonction va retourner l'image du cycle résultant de l'algorithme.
"""
function HK_total_prim_heaviest(chemin::String, image::Bool = false)
     graph = main1(chemin)
     graph_min = prim(graph)
     HK_tree = HK(graph, "prim", find_heaviest_node(graph_min))
     if image
         display(plot_graph(HK_tree))
     end
     println(graphweight(HK_tree))
     HK_tree
end

"""Fonction pratique permettant d'appliquer l'algorithme HK sur un graphe
de recouvrement minimal créé à l'aide de l'algorithme de Kruskal appliqué sur une
instance donné en argument. Le noeud utilisé comme racine est un noeud appartenant
à l'arête la moins lourde du graphe de recouvrement minimal. Si true est donné en argument,
la fonction va retourner l'image du cycle résultant de l'algorithme.
"""
function HK_total_kruskal_lightest(chemin::String, image::Bool = false)
     graph = main1(chemin)
     graph_min = prim(graph)
     HK_tree = HK(graph, "kruskal", find_lightest_node(graph_min))
     if image
         display(plot_graph(HK_tree))
     end
     println(graphweight(HK_tree))
     HK_tree
end

"""Fonction pratique permettant d'appliquer l'algorithme HK sur un graphe
de recouvrement minimal créé à l'aide de l'algorithme de Kruskal appliqué sur une
instance donné en argument. Le noeud utilisé comme racine est un noeud appartenant
à l'arête la plus lourde du graphe de recouvrement minimal. Si true est donné en argument,
la fonction va retourner l'image du cycle résultant de l'algorithme.
"""
function HK_total_kruskal_heaviest(chemin::String, image::Bool = false)
     graph = main1(chemin)
     graph_min = prim(graph)
     HK_tree = HK(graph, "kruskal", find_heaviest_node(graph_min))
     if image
         display(plot_graph(HK_tree))
     end
     println(graphweight(HK_tree))
     HK_tree
end

"""Fonction pratique permettant de créer un tableau avec les poids des différentes
solutions trouvées à l'aide de l'algorithme RSL. Ce tableau compare les résultats
trouvés avec les paramètres prim-lightest à ceux trouvés avec les paramètres
prim-heaviest.
"""
function resultats_prim_lightest_vs_heaviest()
    resultat = zeros(14,3)
    resultat[1,1]=RSL_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp")
    resultat[2,1]=RSL_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp")
    resultat[3,1]=RSL_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp")
    resultat[4,1]=RSL_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp")
    resultat[5,1]=RSL_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp")
    resultat[6,1]=RSL_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp")
    resultat[7,1]=RSL_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp")
    resultat[8,1]=RSL_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp")
    resultat[9,1]=RSL_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp")
    resultat[10,1]=RSL_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp")
    resultat[11,1]=RSL_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp")
    resultat[12,1]=RSL_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp")
    resultat[13,1]=RSL_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp")
    resultat[14,1]=RSL_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp")
    resultat[1,2]=RSL_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp")
    resultat[2,2]=RSL_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp")
    resultat[3,2]=RSL_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp")
    resultat[4,2]=RSL_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp")
    resultat[5,2]=RSL_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp")
    resultat[6,2]=RSL_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp")
    resultat[7,2]=RSL_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp")
    resultat[8,2]=RSL_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp")
    resultat[9,2]=RSL_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp")
    resultat[10,2]=RSL_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp")
    resultat[11,2]=RSL_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp")
    resultat[12,2]=RSL_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp")
    resultat[13,2]=RSL_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp")
    resultat[14,2]=RSL_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp")
    i = 1
    for i in 1:14
        resultat[i,3] = resultat[i,1] - resultat[i,2]
    end
    return resultat
end

"""Fonction pratique permettant de créer un tableau avec les poids des différentes
solutions trouvées à l'aide de l'algorithme RSL. Ce tableau compare les résultats
trouvés avec prim et 4 différentes valeurs de racines(lightest, heaviest, premier, dernier).
"""
function resultats_prim_first_vs_last_vs_lightest_vs_heaviest()
    resultat = zeros(14,4)
    resultat[1,1]=graphweight(RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp"))
    resultat[2,1]=graphweight(RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp"))
    resultat[3,1]=graphweight(RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp"))
    resultat[4,1]=graphweight(RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp"))
    resultat[5,1]=graphweight(RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp"))
    resultat[6,1]=graphweight(RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp"))
    resultat[7,1]=graphweight(RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp"))
    resultat[8,1]=graphweight(RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp"))
    resultat[9,1]=graphweight(RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp"))
    resultat[10,1]=graphweight(RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp"))
    resultat[11,1]=graphweight(RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp"))
    resultat[12,1]=graphweight(RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp"))
    resultat[13,1]=graphweight(RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp"))
    resultat[14,1]=graphweight(RSL_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp"))
    resultat[1,2]=graphweight(RSL_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp"))
    resultat[2,2]=graphweight(RSL_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp"))
    resultat[3,2]=graphweight(RSL_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp"))
    resultat[4,2]=graphweight(RSL_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp"))
    resultat[5,2]=graphweight(RSL_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp"))
    resultat[6,2]=graphweight(RSL_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp"))
    resultat[7,2]=graphweight(RSL_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp"))
    resultat[8,2]=graphweight(RSL_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp"))
    resultat[9,2]=graphweight(RSL_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp"))
    resultat[10,2]=graphweight(RSL_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp"))
    resultat[11,2]=graphweight(RSL_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp"))
    resultat[12,2]=graphweight(RSL_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp"))
    resultat[13,2]=graphweight(RSL_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp"))
    resultat[14,2]=graphweight(RSL_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp"))
    resultat[1,3]=graphweight(RSL_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp"))
    resultat[2,3]=graphweight(RSL_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp"))
    resultat[3,3]=graphweight(RSL_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp"))
    resultat[4,3]=graphweight(RSL_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp"))
    resultat[5,3]=graphweight(RSL_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp"))
    resultat[6,3]=graphweight(RSL_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp"))
    resultat[7,3]=graphweight(RSL_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp"))
    resultat[8,3]=graphweight(RSL_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp"))
    resultat[9,3]=graphweight(RSL_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp"))
    resultat[10,3]=graphweight(RSL_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp"))
    resultat[11,3]=graphweight(RSL_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp"))
    resultat[12,3]=graphweight(RSL_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp"))
    resultat[13,3]=graphweight(RSL_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp"))
    resultat[14,3]=graphweight(RSL_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp"))
    resultat[1,4]=graphweight(RSL_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp"))
    resultat[2,4]=graphweight(RSL_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp"))
    resultat[3,4]=graphweight(RSL_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp"))
    resultat[4,4]=graphweight(RSL_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp"))
    resultat[5,4]=graphweight(RSL_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp"))
    resultat[6,4]=graphweight(RSL_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp"))
    resultat[7,4]=graphweight(RSL_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp"))
    resultat[8,4]=graphweight(RSL_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp"))
    resultat[9,4]=graphweight(RSL_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp"))
    resultat[10,4]=graphweight(RSL_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp"))
    resultat[11,4]=graphweight(RSL_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp"))
    resultat[12,4]=graphweight(RSL_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp"))
    resultat[13,4]=graphweight(RSL_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp"))
    resultat[14,4]=graphweight(RSL_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp"))
    return resultat
end

"""Fonction pratique permettant de créer un tableau avec les poids des différentes
solutions trouvées à l'aide de l'algorithme RSL. Ce tableau compare les résultats
trouvés avec kruskal et 4 différentes valeurs de racines(lightest, heaviest, premier, dernier).
"""
function resultats_kruskal_first_vs_last_vs_lightest_vs_heaviest()
    resultat = zeros(14,4)
    resultat[1,1]=graphweight(RSL_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp"))
    resultat[2,1]=graphweight(RSL_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp"))
    resultat[3,1]=graphweight(RSL_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp"))
    resultat[4,1]=graphweight(RSL_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp"))
    resultat[5,1]=graphweight(RSL_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp"))
    resultat[6,1]=graphweight(RSL_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp"))
    resultat[7,1]=graphweight(RSL_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp"))
    resultat[8,1]=graphweight(RSL_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp"))
    resultat[9,1]=graphweight(RSL_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp"))
    resultat[10,1]=graphweight(RSL_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp"))
    resultat[11,1]=graphweight(RSL_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp"))
    resultat[12,1]=graphweight(RSL_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp"))
    resultat[13,1]=graphweight(RSL_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp"))
    resultat[14,1]=graphweight(RSL_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp"))
    resultat[1,2]=graphweight(RSL_total_kruskal_last("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp"))
    resultat[2,2]=graphweight(RSL_total_kruskal_last("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp"))
    resultat[3,2]=graphweight(RSL_total_kruskal_last("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp"))
    resultat[4,2]=graphweight(RSL_total_kruskal_last("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp"))
    resultat[5,2]=graphweight(RSL_total_kruskal_last("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp"))
    resultat[6,2]=graphweight(RSL_total_kruskal_last("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp"))
    resultat[7,2]=graphweight(RSL_total_kruskal_last("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp"))
    resultat[8,2]=graphweight(RSL_total_kruskal_last("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp"))
    resultat[9,2]=graphweight(RSL_total_kruskal_last("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp"))
    resultat[10,2]=graphweight(RSL_total_kruskal_last("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp"))
    resultat[11,2]=graphweight(RSL_total_kruskal_last("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp"))
    resultat[12,2]=graphweight(RSL_total_kruskal_last("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp"))
    resultat[13,2]=graphweight(RSL_total_kruskal_last("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp"))
    resultat[14,2]=graphweight(RSL_total_kruskal_last("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp"))
    resultat[1,3]=graphweight(RSL_total_kruskal_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp"))
    resultat[2,3]=graphweight(RSL_total_kruskal_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp"))
    resultat[3,3]=graphweight(RSL_total_kruskal_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp"))
    resultat[4,3]=graphweight(RSL_total_kruskal_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp"))
    resultat[5,3]=graphweight(RSL_total_kruskal_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp"))
    resultat[6,3]=graphweight(RSL_total_kruskal_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp"))
    resultat[7,3]=graphweight(RSL_total_kruskal_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp"))
    resultat[8,3]=graphweight(RSL_total_kruskal_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp"))
    resultat[9,3]=graphweight(RSL_total_kruskal_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp"))
    resultat[10,3]=graphweight(RSL_total_kruskal_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp"))
    resultat[11,3]=graphweight(RSL_total_kruskal_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp"))
    resultat[12,3]=graphweight(RSL_total_kruskal_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp"))
    resultat[13,3]=graphweight(RSL_total_kruskal_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp"))
    resultat[14,3]=graphweight(RSL_total_kruskal_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp"))
    resultat[1,4]=graphweight(RSL_total_kruskal_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp"))
    resultat[2,4]=graphweight(RSL_total_kruskal_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp"))
    resultat[3,4]=graphweight(RSL_total_kruskal_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp"))
    resultat[4,4]=graphweight(RSL_total_kruskal_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp"))
    resultat[5,4]=graphweight(RSL_total_kruskal_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp"))
    resultat[6,4]=graphweight(RSL_total_kruskal_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp"))
    resultat[7,4]=graphweight(RSL_total_kruskal_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp"))
    resultat[8,4]=graphweight(RSL_total_kruskal_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp"))
    resultat[9,4]=graphweight(RSL_total_kruskal_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp"))
    resultat[10,4]=graphweight(RSL_total_kruskal_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp"))
    resultat[11,4]=graphweight(RSL_total_kruskal_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp"))
    resultat[12,4]=graphweight(RSL_total_kruskal_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp"))
    resultat[13,4]=graphweight(RSL_total_kruskal_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp"))
    resultat[14,4]=graphweight(RSL_total_kruskal_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp"))
    return resultat
end

"""Fonction pratique permettant de créer un tableau avec les poids des différentes
solutions trouvées à l'aide de l'algorithme HK. Ce tableau compare les résultats
trouvés avec prim et 4 différentes valeurs de racines(lightest, heaviest, premier, dernier).
"""
function resultats_HK_prim_first_vs_last_vs_lightest_vs_heaviest()
    resultat = zeros(14,4)
    resultat[1,1]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp"))
    resultat[2,1]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp"))
    resultat[3,1]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp"))
    resultat[4,1]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp"))
    resultat[5,1]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp"))
    resultat[6,1]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp"))
    resultat[7,1]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp"))
    resultat[8,1]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp"))
    resultat[9,1]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp"))
    resultat[10,1]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp"))
    resultat[11,1]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp"))
    resultat[12,1]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp"))
    resultat[13,1]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp"))
    resultat[14,1]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp"))
    resultat[1,2]=graphweight(HK_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp"))
    resultat[2,2]=graphweight(HK_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp"))
    resultat[3,2]=graphweight(HK_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp"))
    resultat[4,2]=graphweight(HK_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp"))
    resultat[5,2]=graphweight(HK_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp"))
    resultat[6,2]=graphweight(HK_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp"))
    resultat[7,2]=graphweight(HK_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp"))
    resultat[8,2]=graphweight(HK_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp"))
    resultat[9,2]=graphweight(HK_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp"))
    resultat[10,2]=graphweight(HK_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp"))
    resultat[11,2]=graphweight(HK_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp"))
    resultat[12,2]=graphweight(HK_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp"))
    resultat[13,2]=graphweight(HK_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp"))
    resultat[14,2]=graphweight(HK_total_prim_last("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp"))
    resultat[1,3]=graphweight(HK_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp"))
    resultat[2,3]=graphweight(HK_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp"))
    resultat[3,3]=graphweight(HK_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp"))
    resultat[4,3]=graphweight(HK_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp"))
    resultat[5,3]=graphweight(HK_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp"))
    resultat[6,3]=graphweight(HK_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp"))
    resultat[7,3]=graphweight(HK_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp"))
    resultat[8,3]=graphweight(HK_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp"))
    resultat[9,3]=graphweight(HK_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp"))
    resultat[10,3]=graphweight(HK_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp"))
    resultat[11,3]=graphweight(HK_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp"))
    resultat[12,3]=graphweight(HK_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp"))
    resultat[13,3]=graphweight(HK_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp"))
    resultat[14,3]=graphweight(HK_total_prim_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp"))
    resultat[1,4]=graphweight(HK_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp"))
    resultat[2,4]=graphweight(HK_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp"))
    resultat[3,4]=graphweight(HK_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp"))
    resultat[4,4]=graphweight(HK_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp"))
    resultat[5,4]=graphweight(HK_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp"))
    resultat[6,4]=graphweight(HK_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp"))
    resultat[7,4]=graphweight(HK_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp"))
    resultat[8,4]=graphweight(HK_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp"))
    resultat[9,4]=graphweight(HK_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp"))
    resultat[10,4]=graphweight(HK_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp"))
    resultat[11,4]=graphweight(HK_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp"))
    resultat[12,4]=graphweight(HK_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp"))
    resultat[13,4]=graphweight(HK_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp"))
    resultat[14,4]=graphweight(HK_total_prim_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp"))
    return resultat
end

"""Fonction pratique permettant de créer un tableau avec les poids des différentes
solutions trouvées à l'aide de l'algorithme HK. Ce tableau compare les résultats
trouvés avec kruskal et 4 différentes valeurs de racines(lightest, heaviest, premier, dernier).
"""
function resultats_HK_kruskal_first_vs_last_vs_lightest_vs_heaviest()
    resultat = zeros(14,4)
    resultat[1,1]=graphweight(HK_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp"))
    resultat[2,1]=graphweight(HK_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp"))
    resultat[3,1]=graphweight(HK_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp"))
    resultat[4,1]=graphweight(HK_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp"))
    resultat[5,1]=graphweight(HK_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp"))
    resultat[6,1]=graphweight(HK_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp"))
    resultat[7,1]=graphweight(HK_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp"))
    resultat[8,1]=graphweight(HK_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp"))
    resultat[9,1]=graphweight(HK_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp"))
    resultat[10,1]=graphweight(HK_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp"))
    resultat[11,1]=graphweight(HK_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp"))
    resultat[12,1]=graphweight(HK_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp"))
    resultat[13,1]=graphweight(HK_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp"))
    resultat[14,1]=graphweight(HK_total_kruskal("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp"))
    resultat[1,2]=graphweight(HK_total_kruskal_last("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp"))
    resultat[2,2]=graphweight(HK_total_kruskal_last("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp"))
    resultat[3,2]=graphweight(HK_total_kruskal_last("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp"))
    resultat[4,2]=graphweight(HK_total_kruskal_last("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp"))
    resultat[5,2]=graphweight(HK_total_kruskal_last("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp"))
    resultat[6,2]=graphweight(HK_total_kruskal_last("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp"))
    resultat[7,2]=graphweight(HK_total_kruskal_last("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp"))
    resultat[8,2]=graphweight(HK_total_kruskal_last("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp"))
    resultat[9,2]=graphweight(HK_total_kruskal_last("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp"))
    resultat[10,2]=graphweight(HK_total_kruskal_last("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp"))
    resultat[11,2]=graphweight(HK_total_kruskal_last("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp"))
    resultat[12,2]=graphweight(HK_total_kruskal_last("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp"))
    resultat[13,2]=graphweight(HK_total_kruskal_last("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp"))
    resultat[14,2]=graphweight(HK_total_kruskal_last("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp"))
    resultat[1,3]=graphweight(HK_total_kruskal_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp"))
    resultat[2,3]=graphweight(HK_total_kruskal_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp"))
    resultat[3,3]=graphweight(HK_total_kruskal_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp"))
    resultat[4,3]=graphweight(HK_total_kruskal_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp"))
    resultat[5,3]=graphweight(HK_total_kruskal_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp"))
    resultat[6,3]=graphweight(HK_total_kruskal_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp"))
    resultat[7,3]=graphweight(HK_total_kruskal_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp"))
    resultat[8,3]=graphweight(HK_total_kruskal_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp"))
    resultat[9,3]=graphweight(HK_total_kruskal_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp"))
    resultat[10,3]=graphweight(HK_total_kruskal_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp"))
    resultat[11,3]=graphweight(HK_total_kruskal_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp"))
    resultat[12,3]=graphweight(HK_total_kruskal_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp"))
    resultat[13,3]=graphweight(HK_total_kruskal_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp"))
    resultat[14,3]=graphweight(HK_total_kruskal_lightest("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp"))
    resultat[1,4]=graphweight(HK_total_kruskal_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp"))
    resultat[2,4]=graphweight(HK_total_kruskal_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp"))
    resultat[3,4]=graphweight(HK_total_kruskal_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp"))
    resultat[4,4]=graphweight(HK_total_kruskal_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp"))
    resultat[5,4]=graphweight(HK_total_kruskal_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp"))
    resultat[6,4]=graphweight(HK_total_kruskal_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp"))
    resultat[7,4]=graphweight(HK_total_kruskal_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp"))
    resultat[8,4]=graphweight(HK_total_kruskal_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp"))
    resultat[9,4]=graphweight(HK_total_kruskal_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp"))
    resultat[10,4]=graphweight(HK_total_kruskal_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp"))
    resultat[11,4]=graphweight(HK_total_kruskal_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp"))
    resultat[12,4]=graphweight(HK_total_kruskal_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp"))
    resultat[13,4]=graphweight(HK_total_kruskal_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp"))
    resultat[14,4]=graphweight(HK_total_kruskal_heaviest("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp"))
    return resultat
end

"""Fonction pratique permettant d'appliquer l'algorithme HK sur un graphe
de recouvrement minimal créé à l'aide de l'algorithme de Prim appliqué sur une
instance donné en argument. Le noeud utilisé comme racine est le premier noeud
du graphique de recouvrement minimal. Si true est donné en argument, la fonction
va retourner l'image du cycle résultant de l'algorithme. Le pas utilisé pour
l'algorithme HK à une valeur de 1.
"""
function HK_total_prim_step1(chemin::String, image::Bool = false)
    graph = main1(chemin)
    graph_min = prim(graph)
    HK_tree = HK(graph, "prim",nodes(graph_min)[1],1.0)
    if image
        display(plot_graph(HK_tree))
    end
    println(graphweight(HK_tree))
    HK_tree
end

"""Fonction pratique permettant d'appliquer l'algorithme HK sur un graphe
de recouvrement minimal créé à l'aide de l'algorithme de Prim appliqué sur une
instance donné en argument. Le noeud utilisé comme racine est le premier noeud
du graphique de recouvrement minimal. Si true est donné en argument, la fonction
va retourner l'image du cycle résultant de l'algorithme. Le pas utilisé pour
l'algorithme HK à une valeur de 5.
"""
function HK_total_prim_step5(chemin::String, image::Bool = false)
    graph = main1(chemin)
    graph_min = prim(graph)
    HK_tree = HK(graph, "prim",nodes(graph_min)[1],5.0)
    if image
        display(plot_graph(HK_tree))
    end
    println(graphweight(HK_tree))
    HK_tree
end

"""Fonction pratique permettant d'appliquer l'algorithme HK sur un graphe
de recouvrement minimal créé à l'aide de l'algorithme de Prim appliqué sur une
instance donné en argument. Le noeud utilisé comme racine est le premier noeud
du graphique de recouvrement minimal. Si true est donné en argument, la fonction
va retourner l'image du cycle résultant de l'algorithme. Le pas utilisé pour
l'algorithme HK à une valeur de 10.
"""
function HK_total_prim_step10(chemin::String, image::Bool = false)
    graph = main1(chemin)
    graph_min = prim(graph)
    HK_tree = HK(graph, "prim",nodes(graph_min)[1],10.0)
    if image
        display(plot_graph(HK_tree))
    end
    println(graphweight(HK_tree))
    HK_tree
end

"""Fonction pratique permettant de créer un tableau avec les poids des différentes
solutions trouvées à l'aide de l'algorithme HK. Ce tableau compare les résultats
trouvés avec prim et 4 différentes valeurs de racines(lightest, heaviest, premier, dernier).
"""
function resultats_HK_prim_step1_vs_step2_vs_step5_vs_step10()
    resultat = zeros(14,4)
    resultat[1,1]=graphweight(HK_total_prim_step1("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp"))
    resultat[2,1]=graphweight(HK_total_prim_step1("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp"))
    resultat[3,1]=graphweight(HK_total_prim_step1("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp"))
    resultat[4,1]=graphweight(HK_total_prim_step1("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp"))
    resultat[5,1]=graphweight(HK_total_prim_step1("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp"))
    resultat[6,1]=graphweight(HK_total_prim_step1("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp"))
    resultat[7,1]=graphweight(HK_total_prim_step1("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp"))
    resultat[8,1]=graphweight(HK_total_prim_step1("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp"))
    resultat[9,1]=graphweight(HK_total_prim_step1("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp"))
    resultat[10,1]=graphweight(HK_total_prim_step1("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp"))
    resultat[11,1]=graphweight(HK_total_prim_step1("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp"))
    resultat[12,1]=graphweight(HK_total_prim_step1("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp"))
    resultat[13,1]=graphweight(HK_total_prim_step1("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp"))
    resultat[14,1]=graphweight(HK_total_prim_step1("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp"))
    resultat[1,2]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp"))
    resultat[2,2]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp"))
    resultat[3,2]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp"))
    resultat[4,2]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp"))
    resultat[5,2]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp"))
    resultat[6,2]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp"))
    resultat[7,2]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp"))
    resultat[8,2]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp"))
    resultat[9,2]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp"))
    resultat[10,2]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp"))
    resultat[11,2]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp"))
    resultat[12,2]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp"))
    resultat[13,2]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp"))
    resultat[14,2]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp"))
    resultat[1,3]=graphweight(HK_total_prim_step5("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp"))
    resultat[2,3]=graphweight(HK_total_prim_step5("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp"))
    resultat[3,3]=graphweight(HK_total_prim_step5("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp"))
    resultat[4,3]=graphweight(HK_total_prim_step5("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp"))
    resultat[5,3]=graphweight(HK_total_prim_step5("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp"))
    resultat[6,3]=graphweight(HK_total_prim_step5("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp"))
    resultat[7,3]=graphweight(HK_total_prim_step5("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp"))
    resultat[8,3]=graphweight(HK_total_prim_step5("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp"))
    resultat[9,3]=graphweight(HK_total_prim_step5("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp"))
    resultat[10,3]=graphweight(HK_total_prim_step5("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp"))
    resultat[11,3]=graphweight(HK_total_prim_step5("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp"))
    resultat[12,3]=graphweight(HK_total_prim_step5("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp"))
    resultat[13,3]=graphweight(HK_total_prim_step5("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp"))
    resultat[14,3]=graphweight(HK_total_prim_step5("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp"))
    resultat[1,4]=graphweight(HK_total_prim_step10("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp"))
    resultat[2,4]=graphweight(HK_total_prim_step10("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp"))
    resultat[3,4]=graphweight(HK_total_prim_step10("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp"))
    resultat[4,4]=graphweight(HK_total_prim_step10("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp"))
    resultat[5,4]=graphweight(HK_total_prim_step10("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp"))
    resultat[6,4]=graphweight(HK_total_prim_step10("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp"))
    resultat[7,4]=graphweight(HK_total_prim_step10("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp"))
    resultat[8,4]=graphweight(HK_total_prim_step10("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp"))
    resultat[9,4]=graphweight(HK_total_prim_step10("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp"))
    resultat[10,4]=graphweight(HK_total_prim_step10("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp"))
    resultat[11,4]=graphweight(HK_total_prim_step10("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp"))
    resultat[12,4]=graphweight(HK_total_prim_step10("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp"))
    resultat[13,4]=graphweight(HK_total_prim_step10("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp"))
    resultat[14,4]=graphweight(HK_total_prim_step10("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp"))
    return resultat
end

"""Fonction pratique permettant d'appliquer l'algorithme HK sur un graphe
de recouvrement minimal créé à l'aide de l'algorithme de Prim appliqué sur une
instance donné en argument. Le noeud utilisé comme racine est le premier noeud
du graphique de recouvrement minimal. Si true est donné en argument, la fonction
va retourner l'image du cycle résultant de l'algorithme. Le critère d'arrêt utilisé
est 5 itérations.
"""
function HK_total_prim_iter5(chemin::String, image::Bool = false)
    graph = main1(chemin)
    graph_min = prim(graph)
    HK_tree = HK(graph, "prim",nodes(graph_min)[1],2.0,5)
    if image
        display(plot_graph(HK_tree))
    end
    println(graphweight(HK_tree))
    HK_tree
end

"""Fonction pratique permettant d'appliquer l'algorithme HK sur un graphe
de recouvrement minimal créé à l'aide de l'algorithme de Prim appliqué sur une
instance donné en argument. Le noeud utilisé comme racine est le premier noeud
du graphique de recouvrement minimal. Si true est donné en argument, la fonction
va retourner l'image du cycle résultant de l'algorithme. Le critère d'arrêt utilisé
est 10 itérations.
"""
function HK_total_prim_iter10(chemin::String, image::Bool = false)
    graph = main1(chemin)
    graph_min = prim(graph)
    HK_tree = HK(graph, "prim",nodes(graph_min)[1],2.0,10)
    if image
        display(plot_graph(HK_tree))
    end
    println(graphweight(HK_tree))
    HK_tree
end

"""Fonction pratique permettant d'appliquer l'algorithme HK sur un graphe
de recouvrement minimal créé à l'aide de l'algorithme de Prim appliqué sur une
instance donné en argument. Le noeud utilisé comme racine est le premier noeud
du graphique de recouvrement minimal. Si true est donné en argument, la fonction
va retourner l'image du cycle résultant de l'algorithme. Le critère d'arrêt utilisé
est 50 itérations.
"""
function HK_total_prim_iter50(chemin::String, image::Bool = false)
    graph = main1(chemin)
    graph_min = prim(graph)
    HK_tree = HK(graph, "prim",nodes(graph_min)[1],2.0,50)
    if image
        display(plot_graph(HK_tree))
    end
    println(graphweight(HK_tree))
    HK_tree
end

"""Fonction pratique permettant de créer un tableau avec les poids des différentes
solutions trouvées à l'aide de l'algorithme HK. Ce tableau compare les résultats
trouvés avec prim et 4 différentes valeurs de racines(lightest, heaviest, premier, dernier).
"""
function resultats_HK_prim_iter5_vs_iter10_vs_iter20_vs_iter50()
    resultat = zeros(14,4)
    resultat[1,2]=graphweight(HK_total_prim_iter10("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp"))
    resultat[2,2]=graphweight(HK_total_prim_iter10("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp"))
    resultat[3,2]=graphweight(HK_total_prim_iter10("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp"))
    resultat[4,2]=graphweight(HK_total_prim_iter10("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp"))
    resultat[5,2]=graphweight(HK_total_prim_iter10("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp"))
    resultat[6,2]=graphweight(HK_total_prim_iter10("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp"))
    resultat[7,2]=graphweight(HK_total_prim_iter10("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp"))
    resultat[8,2]=graphweight(HK_total_prim_iter10("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp"))
    resultat[9,2]=graphweight(HK_total_prim_iter10("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp"))
    resultat[10,2]=graphweight(HK_total_prim_iter10("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp"))
    resultat[11,2]=graphweight(HK_total_prim_iter10("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp"))
    resultat[12,2]=graphweight(HK_total_prim_iter10("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp"))
    resultat[13,2]=graphweight(HK_total_prim_iter10("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp"))
    resultat[14,2]=graphweight(HK_total_prim_iter10("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp"))
    resultat[1,3]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp"))
    resultat[2,3]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp"))
    resultat[3,3]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp"))
    resultat[4,3]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp"))
    resultat[5,3]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp"))
    resultat[6,3]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp"))
    resultat[7,3]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp"))
    resultat[8,3]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp"))
    resultat[9,3]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp"))
    resultat[10,3]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp"))
    resultat[11,3]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp"))
    resultat[12,3]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp"))
    resultat[13,3]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp"))
    resultat[14,3]=graphweight(HK_total_prim("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp"))
    resultat[1,1]=graphweight(HK_total_prim_iter5("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp"))
    resultat[2,1]=graphweight(HK_total_prim_iter5("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp"))
    resultat[3,1]=graphweight(HK_total_prim_iter5("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp"))
    resultat[4,1]=graphweight(HK_total_prim_iter5("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp"))
    resultat[5,1]=graphweight(HK_total_prim_iter5("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp"))
    resultat[6,1]=graphweight(HK_total_prim_iter5("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp"))
    resultat[7,1]=graphweight(HK_total_prim_iter5("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp"))
    resultat[8,1]=graphweight(HK_total_prim_iter5("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp"))
    resultat[9,1]=graphweight(HK_total_prim_iter5("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp"))
    resultat[10,1]=graphweight(HK_total_prim_iter5("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp"))
    resultat[11,1]=graphweight(HK_total_prim_iter5("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp"))
    resultat[12,1]=graphweight(HK_total_prim_iter5("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp"))
    resultat[13,1]=graphweight(HK_total_prim_iter5("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp"))
    resultat[14,1]=graphweight(HK_total_prim_iter5("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp"))
    resultat[1,4]=graphweight(HK_total_prim_iter50("/home/jonathan/mth6412b-starter-code/instances/stsp/bayg29.tsp"))
    resultat[2,4]=graphweight(HK_total_prim_iter50("/home/jonathan/mth6412b-starter-code/instances/stsp/bays29.tsp"))
    resultat[3,4]=graphweight(HK_total_prim_iter50("/home/jonathan/mth6412b-starter-code/instances/stsp/brazil58.tsp"))
    resultat[4,4]=graphweight(HK_total_prim_iter50("/home/jonathan/mth6412b-starter-code/instances/stsp/brg180.tsp"))
    resultat[5,4]=graphweight(HK_total_prim_iter50("/home/jonathan/mth6412b-starter-code/instances/stsp/dantzig42.tsp"))
    resultat[6,4]=graphweight(HK_total_prim_iter50("/home/jonathan/mth6412b-starter-code/instances/stsp/fri26.tsp"))
    resultat[7,4]=graphweight(HK_total_prim_iter50("/home/jonathan/mth6412b-starter-code/instances/stsp/gr17.tsp"))
    resultat[8,4]=graphweight(HK_total_prim_iter50("/home/jonathan/mth6412b-starter-code/instances/stsp/gr21.tsp"))
    resultat[9,4]=graphweight(HK_total_prim_iter50("/home/jonathan/mth6412b-starter-code/instances/stsp/gr24.tsp"))
    resultat[10,4]=graphweight(HK_total_prim_iter50("/home/jonathan/mth6412b-starter-code/instances/stsp/gr48.tsp"))
    resultat[11,4]=graphweight(HK_total_prim_iter50("/home/jonathan/mth6412b-starter-code/instances/stsp/gr120.tsp"))
    resultat[12,4]=graphweight(HK_total_prim_iter50("/home/jonathan/mth6412b-starter-code/instances/stsp/hk48.tsp"))
    resultat[13,4]=graphweight(HK_total_prim_iter50("/home/jonathan/mth6412b-starter-code/instances/stsp/pa561.tsp"))
    resultat[14,4]=graphweight(HK_total_prim_iter50("/home/jonathan/mth6412b-starter-code/instances/stsp/swiss42.tsp"))
    return resultat
end
