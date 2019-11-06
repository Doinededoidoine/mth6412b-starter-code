using Statistics
using Printf

include(joinpath(@__DIR__, "..", "phase1", "main1.jl"))
include("main3.jl")

# List of stsp instances
const STSP_DIR = joinpath(@__DIR__, "..", "..", "instances", "stsp")
const STSP = readdir(STSP_DIR)

const GRAPHS = [main1(joinpath(STSP_DIR, finst)) for finst in STSP]

"""
    benchmark_kruskal(graphs)
Run Kruskal's algorithm on each graph `g ∈ graphs`.
"""
function benchmark_kruskal(graphs)
    for g in graphs
        mst = kruskal(g)
    end
    return nothing
end

# First round for compilation
benchmark_kruskal(GRAPHS)

# Second round
# Here we record times
N = 16  # Number of runs. Higher yields more accurate results
T = Float64[]
for i in 1:N
    t = @elapsed benchmark_kruskal(GRAPHS)
    push!(T, t)
end

μ = mean(T)  # mean
σ = std(T)   # standard deviation

# Print results
@printf "Total time: %8.4f\n" sum(T)
@printf "Min   time: %8.4f\n" minimum(T)
@printf "CI-lo time: %8.4f\n" μ - (2.0 / sqrt(N)) * σ
@printf "Mean  time: %8.4f\n" μ
@printf "CI-up time: %8.4f\n" μ + (2.0 / sqrt(N)) * σ
@printf "Max   time: %8.4f\n" maximum(T)
