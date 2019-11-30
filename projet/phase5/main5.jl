using Random, FileIO, Images, ImageView, ImageMagick
include(joinpath(@__DIR__, "..", "phase4", "main4.jl"))
include(joinpath(@__DIR__, "..", "phase5", "tools.jl"))

"Créer le fichier image corresondant à la tournée trouvée avec Prim"
function total_image_prim(chemin::String,input_name::String,name::String = "test.tour",output_name::String = "test.png") where T
	RSL_tree = RSL_total_prim(chemin)
	RSL_liste = liste_tour(RSL_tree)
	write_tour(name, RSL_liste, graphweight(RSL_tree))
	reconstruct_picture(name, input_name, output_name)
end

"Créer le fichier image corresondant à la tournée trouvée avec Kruskal"
function total_image_kruskal(chemin::String,input_name::String,name::String = "test.tour",output_name::String = "test.png") where T
	RSL_tree = RSL_total_kruskal(chemin)
	RSL_liste = liste_tour(RSL_tree)
	write_tour(name, RSL_liste, graphweight(RSL_tree))
	reconstruct_picture(name, input_name, output_name)
end

"""Crée une liste de noeud associé à une tournée"""
function liste_tour(graph::AbstractGraph{T}) where T
	liste = zeros(Int64,length(nodes(graph))+1)
	liste[1] = 0
	for (i,node) in enumerate(nodes(graph))
		liste[i+1] = (parse(Int64,name(nodes(graph)[i])))
	end
	liste
end
