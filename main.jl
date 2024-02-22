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
using BenchmarkTools
using DelimitedFiles

# select seed for reproducibility
Random.seed!(1)

# initialize all constant types and parse them from the input files
number_of_cubes::Integer = 0
dimensionality::Integer = 0
sidelength::Integer = 0
numberoffaces_inacube::Integer = 0
centerof_baseface::Integer = 0
deltaR::Float64 = 0 # radius of sphere on whose surface each random walker is generated
steplength::Integer = 0
inputs = readdlm("inputs.txt",'=')
for ii in axes(inputs,1)
    for jj=1:(size(inputs,2)-1)
        if inputs[ii,jj] == "TOTAL_NUMBER_OF_CUBES"
            global number_of_cubes = inputs[ii,jj+1]
        elseif inputs[ii,jj] == "DIMENSIONALITY"
            global dimensionality = inputs[ii,jj+1]
        elseif inputs[ii,jj] == "SIDELENGTH"
            global sidelength = inputs[ii,jj+1]
        elseif inputs[ii,jj] == "FACES_IN_A_SINGLE_CUBE"
            global numberoffaces_inacube = inputs[ii,jj+1]
        elseif inputs[ii,jj] == "CENTER_OF_BASE_FACE"
            global centerof_baseface = inputs[ii,jj+1]
        elseif inputs[ii,jj] == "ADDITIONAL_RADIUS_OF_SPHERE_WHERE_RANDOM_WALKER_IS_GENERATED"
            global deltaR = inputs[ii,jj+1]
        elseif inputs[ii,jj] == "STEPLENGTH"
            global steplength = inputs[ii,jj+1]
        elseif inputs[ii,jj] == "ATTACHING_DISTANCE"
            global attaching_distance = inputs[ii,jj+1]
        end
    end 
end

if dimensionality != 3
    error("dimensionality MUST be 3")
end

if numberoffaces_inacube != 6
    error("numberoffaces_inacube MUST be 6")
end

if centerof_baseface != 1
    error("centerof_baseface MUST be 1")
end

if attaching_distance!=steplength
    error("attaching_distance MUST be equal to steplength: random walk is on a 3d-lattice")
end
# initialize array to hold position of cubes in aggregate
cubes::Matrix{Integer} = zeros(number_of_cubes,dimensionality)

### THE BASE CUBE ####
# define the position of the faces of the base cube of the aggregate
position_offaces::Matrix{Integer} = positionoffaces_of_basecube(numberoffaces_inacube,dimensionality,centerof_baseface)
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

############ implement Individually-added aggregation routine ########################
# Reference: https://journals.aps.org/prfluids/abstract/10.1103/PhysRevFluids.5.044305

# call function that returns the final aggregate
final_position::Matrix{Integer} = individuallyadded_aggregate!(cubes,number_of_cubes,dimensionality,deltaR,steplength,attaching_distance)
println("POSITION OF CUBES IN AGGREGATE:")
for ii in axes(cubes,1)
    println(final_position[ii,:])
end

# print position of external faces
externalfaces::Matrix{Integer} = build_externalfacesofaggregate(numberoffaces_inacube,dimensionality,sidelength,centerof_baseface,final_position)
println("POSITION OF EXTERNAL FACES OF AGGREGATE:")
for ii in axes(externalfaces,1)
    println(externalfaces[ii,:])
end

# print orientation of external faces
orientationof_externalfaces::Vector{String} = getorientation_externalfacesofaggregate(numberoffaces_inacube,dimensionality,sidelength,centerof_baseface,final_position)
println("ORIENTATION OF EXTERNAL FACES OF AGGREGATE:")
for ii in eachindex(orientationof_externalfaces)
    println(orientationof_externalfaces[ii])
end

# single-layer x_terms for one face (debugging purposes)
const_res, x_res = build_singlelayermatrix(externalfaces[4,:],externalfaces[6,:],orientationof_externalfaces[6],dimensionality)
println("CONSTANT TERMS")
for ii in axes(const_res,1)
    for jj in axes(const_res,2)
        print(const_res[ii,jj])
        print("\t")
    end
    println()
end
println("X TERMS")
for ii in axes(x_res,1)
    for jj in axes(x_res,2)
        print(x_res[ii,jj])
        print("\t")
    end
    println()
end

# full single layer potential 
println("SINGLE LAYER")
# single-layer
LHS_single = zeros(size(externalfaces,1)*dimensionality,size(externalfaces,1)*dimensionality)
@btime begin
Threads.@threads for ii in axes(externalfaces,1)
    integration_point::Vector{Integer} = externalfaces[ii,:]
    for jj in axes(externalfaces,1)
        constant_terms,x_terms = build_singlelayermatrix(integration_point,externalfaces[jj,:],orientationof_externalfaces[jj],dimensionality)
        LHS_single[dimensionality*(ii-1)+1:dimensionality*ii,dimensionality*(jj-1)+1:dimensionality*jj] = (constant_terms+x_terms);
    end
end
end

