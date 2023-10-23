module Metrics

using Statistics
using LinearAlgebra

export compute_centerofmass
export compute_maxdistance

function compute_centerofmass(centerofcubes_aggregate::VecOrMat{Integer},dimensionality::Integer)
    
    centerofmass_aggregate::Vector{Float64} = zeros(dimensionality)
    if size(centerofcubes_aggregate,1) == 1
        centerofmass_aggregate = vec(centerofcubes_aggregate) # vec converts Matrix type to vector type
    else
        centerofmass_aggregate = vec(mean(centerofcubes_aggregate,dims=1))
    end

    return centerofmass_aggregate
end

function compute_maxdistance(cubes_in_aggregate::VecOrMat{Integer},center_of_mass::Vector{Float64})
    
    distance_from_center::Vector{Float64} = zeros(size(cubes_in_aggregate,1));
    for any_cube in axes(cubes_in_aggregate,1)
        distance_from_center[any_cube] = norm(cubes_in_aggregate[any_cube,:]-center_of_mass)
    end
    max_distance::Float64 = maximum(distance_from_center)

    return max_distance
end

end