"""
    loadflow(flow::Function, xmin::Float64, xmax::Float64, ymin::Float64, ymax::Float64, nx::Int, ny::Int, patch::Bool)

Load a flow from a function, with specified domain and resolution. If `patch` is true, the v-components of the flow are set to zero at the top and bottom boundaries to enforce a "slip condition".
"""
function loadflow(flow::Function, xmin::Float64, xmax::Float64, ymin::Float64, ymax::Float64, nx::Int, ny::Int, patch::Bool)
    xs = range(xmin, xmax, length=nx)
    ys = range(ymin, ymax, length=ny)

    V = Array{SVector{2, Float64}}(undef, nx, ny)
    for i in 1:nx, j in 1:ny
        V[i, j] = flow(xs[i], ys[j])
    end

    if patch # NOTE: remove v-components at top and bottom layers to enforce "slip condition"
        vslice = @view V[:, 1]
        map!(v -> SVector(v[1], 0.0f0), vslice, vslice)

        vslice = @view V[:, end]
        map!(v -> SVector(v[1], 0.0f0), vslice, vslice)
    end

    itp = extrapolate(
        scale(interpolate(V, BSpline(Linear())), xs, ys),
        Flat()
    )
    return VCFlowData.InterpolatedFlow(itp)
end

# flow, domain, resolution
function loadflow(flow::Function, xmin::Float64, xmax::Float64, ymin::Float64, ymax::Float64, nx::Int, ny::Int; patch::Bool = false)
    return loadflow(flow, xmin, xmax, ymin, ymax, nx, ny, patch)
end

# flow, domain
function loadflow(flow::Function, xmin::Float64, xmax::Float64, ymin::Float64, ymax::Float64; nx::Int = Int((xmax - xmin) * 100) + 1, ny::Int = Int((ymax - ymin) * 100) + 1, patch::Bool = false)
    return loadflow(flow, xmin, xmax, ymin, ymax, nx, ny, patch)
end

# flow, resolution
function loadflow(flow::Function, nx::Int, ny::Int; xmin::Float64 = -10.0, xmax::Float64 = 10.0, ymin::Float64 = -10.0, ymax::Float64 = 10.0, patch::Bool = false)
    return loadflow(flow, xmin, xmax, ymin, ymax, nx, ny, patch)
end

# flow
function loadflow(flow::Function; xmin::Float64 = -10.0, xmax::Float64 = 10.0, ymin::Float64 = -10.0, ymax::Float64 = 10.0, nx::Int = Int((xmax - xmin) * 100) + 1, ny::Int = Int((ymax - ymin) * 100) + 1, patch::Bool = false)
    return loadflow(flow, xmin, xmax, ymin, ymax, nx, ny, patch)
end

# flow, domain, patch
function loadflow(flow::Function, xmin::Float64, xmax::Float64, ymin::Float64, ymax::Float64, patch::Bool; nx::Int = Int((xmax - xmin) * 100) + 1, ny::Int = Int((ymax - ymin) * 100) + 1)
    return loadflow(flow, xmin, xmax, ymin, ymax, nx, ny, patch)
end

# flow, resolution, patch
function loadflow(flow::Function, nx::Int, ny::Int, patch::Bool; xmin::Float64 = -10.0, xmax::Float64 = 10.0, ymin::Float64 = -10.0, ymax::Float64 = 10.0)
    return loadflow(flow, xmin, xmax, ymin, ymax, nx, ny, patch)
end

# flow, patch
function loadflow(flow::Function, patch::Bool; xmin::Float64 = -10.0, xmax::Float64 = 10.0, ymin::Float64 = -10.0, ymax::Float64 = 10.0, nx::Int = Int((xmax - xmin) * 100) + 1, ny::Int = Int((ymax - ymin) * 100) + 1)
    return loadflow(flow, xmin, xmax, ymin, ymax, nx, ny, patch)
end

"""
    loadflow(flow::Matrix, xmin::Float64, xmax::Float64, ymin::Float64, ymax::Float64, patch::Bool)
    
Load a flow from a matrix, with specified domain and patch.
"""
function loadflow(flow::Matrix, xmin::Float64, xmax::Float64, ymin::Float64, ymax::Float64, patch::Bool)
    nx, ny = size(flow)
    xs = range(xmin, xmax, length=nx)
    ys = range(ymin, ymax, length=ny)

    if patch # NOTE: remove v-components at top and bottom layers to enforce "slip condition"
        vslice = @view flow[:, 1]
        map!(v -> SVector(v[1], 0.0f0), vslice, vslice)

        vslice = @view flow[:, end]
        map!(v -> SVector(v[1], 0.0f0), vslice, vslice)
    end
    
    itp = extrapolate(
        scale(interpolate(flow, BSpline(Linear())), xs, ys),
        Flat()
    )
    return VCFlowData.InterpolatedFlow(itp)
end

function loadflow(flow::Matrix, xmin::Float64, xmax::Float64, ymin::Float64, ymax::Float64; patch::Bool = false)
    return loadflow(flow, xmin, xmax, ymin, ymax, patch)
end

"""
    loadflow(file::String, patch::Bool)

Load a flow from a NetCDF file using a string from VCDataSets, with optional patch.
"""
function loadflow(file::String, patch::Bool)
    NCDataset(VCDataSets.file(filename=file)) do ds
        xs = ds["xdim"][:] :: Vector{Float32}
        ys = ds["ydim"][:] :: Vector{Float32}

        u = ds["u"]
        v = ds["v"]

        V = SVector{2,Float64}.(
            u[:, :, 1],
            v[:, :, 1]
        ) :: Array{SVector{2,Float64}, 2}

        if patch # NOTE: remove v-components at top and bottom layers to enforce "slip condition"
            vslice = @view V[:, 1]
            map!(v -> SVector(v[1], 0.0f0), vslice, vslice)

            vslice = @view V[:, end]
            map!(v -> SVector(v[1], 0.0f0), vslice, vslice)
        end

        @assert size(V) == (length(xs), length(ys))
        rxs, rys = VCFlowData.asrange(xs), VCFlowData.asrange(ys)

        itp = extrapolate(
            scale(interpolate(V, BSpline(Linear())), rxs, rys), 
            Flat()
        )
        return VCFlowData.InterpolatedFlow(itp)
    end
end

function loadflow(file::String; patch::Bool=false)
    return loadflow(file, patch)
end

export loadflow