include("RandomWalk.jl")
include("Metrics.jl")
include("BuildAggregate.jl")
using .RandomWalk
using .Metrics
using .BuildAggregate
using Random


# select seed for reproducibility
Random.seed!(1)

# select how many cubes will be in the aggregate
number_of_cubes::Integer = 10
# select dimensionality -> MUST BE 3
dimensionality::Integer = 3
if dimensionality != 3
    error("dimensionality MUST be 3")
end
# initialize array to hold position of cubes in aggregate
cubes::Matrix{Integer} = zeros(number_of_cubes,dimensionality)

# implement Individually-added aggregation routine
# Refeference: https://journals.aps.org/prfluids/abstract/10.1103/PhysRevFluids.5.044305
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
for ii in axes(cubes,1)
    println(final_position[ii,:])
end