module TopologicalSkeleton

    using LinearAlgebra
    using StaticArrays
    using Interpolations
    using VCFlowData
    using NCDatasets
    using VCDataSets
    using RK43
    using CairoMakie

    include("decoder.jl")
    include("abstracttopology.jl")
    include("interpolatedtopology.jl")
    include("plot.jl")

end
