import Base.show
import Base.isless
import Base.==

"""Type abstrait dont d'autres types de noeuds dériveront."""
abstract type AbstractNode{T} end

"""Type représentant les noeuds d'un graphe.

Exemple:

        noeud = Node("James", [π, exp(1)])
        noeud = Node("Kirk", "guitar")
        noeud = Node("Lars", 2)

        ATTENTION !
        Les algorithmes utilisant cette implémentation d'un noeud différencient
        parfois deux noeuds par leur attribut "name". Il est donc essentiel que
        tous les noeuds d'un graphe aient des attributs "name" différents pour
        que ces algorithmes fonctionnent correctement.

"""
mutable struct Node{T} <: AbstractNode{T}
  name::String
  data::T
  min_weight::Float64
end

"""Par défaut, un noeud instancié sans poids possède un poids infini."""
Node{T}(name::String, data::T) where T = Node{T}(name, data, Inf)

# on présume que tous les noeuds dérivant d'AbstractNode
# posséderont des champs `name`, `data` et `min_weight`.

"""Renvoie le nom du noeud."""
name(node::AbstractNode) = node.name

"""Renvoie les données contenues dans le noeud."""
data(node::AbstractNode) = node.data

"""Renvoie le poids minimal du noeud."""
min_weight(node::AbstractNode) = node.min_weight

"""Affecte une nouvelle valeur au poids d'un noeud."""
function set_min_weight!(node::AbstractNode, weight::Float64)
  node.min_weight = weight
  node
end

"""Renvoie un booléen : true si le poids du noeud 1 est inférieur au poids du noeud 2, false sinon."""
isless(node1::AbstractNode, node2::AbstractNode) = min_weight(node1) < min_weight(node2)

"""Renvoie un booléen : true si le poids du noeuds 1 est égal au poids du noeud 2, false sinon."""
==(node1::AbstractNode, node2::AbstractNode) = min_weight(node1) == min_weight(node2)

"""Affiche un noeud."""
function show(node::AbstractNode)
  println("Node ", name(node), ", data: ", data(node), ", min_weight : ", min_weight(node))
end
