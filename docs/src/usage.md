# Usage

Load a flow via formula
```
using StaticArrays
formula(x, y) = @SVector [x^3 - x, (x - 0.5) * y]
flow = loadflow(formula, -2.0, 2.0, -2.0, 2.0, 401, 401)
```
Load a flow from a matrix
```
using StaticArrays
matrix = [@SVector rand(Float64, 2) for i in 1:401, j in 1:401]
flow = loadflow(matrix, -2.0, 2.0, -2.0, 2.0)
```
Load a flow from NetCDF file
```
flow = loadflow("pipedcylinder2d.nc")
```
Plot the resulting Vector field and the topological skeleton with critical points, boundary switch points and separatricees
```
plotskeleton(flow)
```
