import Base.length, Base.push!, Base.popfirst!
import Base.show

"""Type abstrait dont d'autres types de files dériveront."""
abstract type AbstractQueue{T} end


"""File de priorité."""
mutable struct PriorityQueue{T <: AbstractNode} <: AbstractQueue{T}
    items::Vector{T}
end

PriorityQueue{T}() where T = PriorityQueue(T[])

"""Retire et renvoie l'élément ayant le poids le plus faible."""
function poplast!(q::PriorityQueue)
    lowest = q.items[1]
    for item in q.items[2:end]
        if item < lowest
            lowest = item
        end
    end
    idx = findfirst(x -> x == lowest, q.items)
    deleteat!(q.items, idx)
    lowest
end

"""Retire et renvoie l'élément ayant le poids le plus élevé."""
function popfirst!(q::PriorityQueue)
    highest = q.items[1]
    for item in q.items[2:end]
        if item > highest
            highest = item
        end
    end
    idx = findfirst(x -> x == highest, q.items)
    deleteat!(q.items, idx)
    highest
end

"""Retire et renvoie l'objet du dessus de la queue."""
pop!(q::PriorityQueue) = pop!(q.items)

"""Consulte l'objet du dessus de la queue."""
top(q::PriorityQueue) = q.items[end]

"""Indique si la queue est vide."""
is_empty(q::PriorityQueue) = length(q.items) == 0

"""Ajoute item sur le dessus de la queue s."""
function push!(q::PriorityQueue, item::T) where T
    push!(q.items, item)
    q
end
