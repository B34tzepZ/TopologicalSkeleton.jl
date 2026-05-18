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

end
