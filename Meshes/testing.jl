#!/bin/julia
println("Testing Meshes")

using Meshes
using WGLMakie

points = Meshes.Point2[(0,0), (1,0), (0,1), (1,1), (0.25,0.5), (0.75,0.5)]
tris  = connect.([(1,5,3), (4,6,2)], Triangle)
quads = connect.([(1,2,6,5), (4,3,5,6)], Quadrangle)
#mesh = SimpleMesh(points, [tris; quads])
mesh = SimpleMesh(points, [tris; quads])
#col  = collect(elements(mesh))
viz(mesh,showfacets = true)

k = 2
p = partition(mesh,UniformPartition(k))
#viz(p[1])

g = CartesianGrid(100, 100)
#viz(g, showfacets=true)

k = 2      # No of partitions
p = partition(g,UniformPartition(k))
viz(p[1])

println("Press any key to continue")
x = readline()

viz(p[2])

top = g.topology

Ad = Adjacency{1}(top)


# Another example
points = Meshes.Point2[(0,0),(1,0),(0,1),(1,1),(0.25,0.5),(0.75,0.5)]

# connect the points into N-gon
connec = connect.([(1,2,6,5),(2,4,6),(4,3,5,6),(3,1,5)], Ngon)

# 2D mesh made of N-gon elements
mesh = SimpleMesh(points, connec)

viz(mesh, showfacets = true)



