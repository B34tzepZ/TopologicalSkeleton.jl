using Documenter
using TopologicalSkeleton

makedocs(
    sitename = "TopologicalSkeleton.jl",
    modules = [TopologicalSkeleton],
    format = Documenter.HTML(),
    pages = [
        "Installation" => "index.md",
        "Usage" => "usage.md",
        "Functions" => "functions.md",
    ]
)