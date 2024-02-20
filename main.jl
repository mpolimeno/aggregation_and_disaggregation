include("aggregation/RandomWalk.jl")
include("aggregation/Metrics.jl")
include("aggregation/BuildAggregate.jl")
include("aggregation/BaseCube.jl")
include("aggregation/BuildFaces.jl")
include("disaggregation/SingleLayer.jl")
using .RandomWalk
using .Metrics
using .BuildAggregate
using .BaseCube
using .BuildFaces
using .SingleLayer

using Random


# select seed for reproducibility
Random.seed!(1)

# select how many cubes will be in the aggregate
number_of_cubes::Integer = 2
# select dimensionality -> MUST BE 3
dimensionality::Integer = 3
if dimensionality != 3
    error("dimensionality MUST be 3")
end
# number of faces in a single cube
numberoffaces_inacube::Integer = 6;
if numberoffaces_inacube != 6
    error("numberoffaces_inacube MUST be 6")
end
# location of center of faces of base cube (centered at the origin)
centerof_baseface::Integer = 1;
if centerof_baseface != 1
    error("centerof_baseface MUST be 1")
end

# initialize array to hold position of cubes in aggregate
cubes::Matrix{Integer} = zeros(number_of_cubes,dimensionality)

### THE BASE CUBE ####
# define the position of the faces of the base cube of the aggregate
position_offaces::Matrix{Integer} = positionoffaces_of_basecube(numberoffaces_inacube,dimensionality,centerof_baseface)
# print position of faces of base cube
println("POSITION OF FACES OF BASE CUBE:")
for jj in axes(position_offaces,1)
    println(position_offaces[jj,:])
end
# define orientation of faces of the base cube of the aggregate
println("ORIENTATION OF FACES OF BASE CUBE:")
for kk in axes(position_offaces,1)
    orientationoffaces::String = orientationoffaces_ofbasecube(kk)
    println(orientationoffaces)
end

# implement Individually-added aggregation routine
# Reference: https://journals.aps.org/prfluids/abstract/10.1103/PhysRevFluids.5.044305
# constants for the routine
deltaR::Float64 = 10 # radius of sphere on whose surface each random walker is generated
steplength::Integer = 2
attaching_distance::Integer = steplength # must be equal to step-length for 3d-lattice random walk
if attaching_distance!=steplength
    error("attaching_distance MUST be equal to steplength: random walk is on a 3d-lattice")
end
# call function that returns the final aggregate
final_position::Matrix{Integer} = individuallyadded_aggregate!(cubes,number_of_cubes,dimensionality,deltaR,steplength,attaching_distance)

# print result nicely
println("POSITION OF CUBES IN AGGREGATE:")
for ii in axes(cubes,1)
    println(final_position[ii,:])
end

# print position of external faces
externalfaces::Matrix{Integer} = build_externalfacesofaggregate(numberoffaces_inacube,dimensionality,centerof_baseface,final_position)
println("POSITION OF EXTERNAL FACES OF AGGREGATE:")
for ii in axes(externalfaces,1)
    println(externalfaces[ii,:])
end

# print orientation of external faces
orientationof_externalfaces::Vector{String} = getorientation_externalfacesofaggregate(numberoffaces_inacube,dimensionality,centerof_baseface,final_position)
println("ORIENTATION OF EXTERNAL FACES OF AGGREGATE:")
for ii in eachindex(orientationof_externalfaces)
    println(orientationof_externalfaces[ii])
end

# single-layer
res = build_singlelayermatrix(externalfaces[4,:],externalfaces[4,:],orientationof_externalfaces[4],dimensionality)
for ii in axes(res,1)
    for jj in axes(res,2)
        print(res[ii,jj])
        print("\t")
    end
    println()
end