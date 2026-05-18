using Documenter
using TopologicalSkeleton

makedocs(
    sitename = "TopologicalSkeleton.jl",
    modules = [TopologicalSkeleton],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
    ]
)