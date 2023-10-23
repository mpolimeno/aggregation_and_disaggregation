module RandomWalk 

export select_initialposition
export throwdice
export stepvector
export takestep!

function select_initialposition(centerofmass_aggregate::Vector{Float64},dimensionality::Integer,radius::Float64,randomnumber::Float64)

    polar_angle = acos(2*randomnumber-1)
    azimuthal_angle = randomnumber*2*pi

    random_initialposition::Vector{Float64} = zeros(dimensionality)
    random_initialposition[1] = centerofmass_aggregate[1] + radius*cos(azimuthal_angle)*sin(polar_angle)
    random_initialposition[2] = centerofmass_aggregate[2] + radius*sin(azimuthal_angle)*sin(polar_angle)
    random_initialposition[3] = centerofmass_aggregate[3] + radius*cos(polar_angle) 

    walkerposition::Vector{Integer} = zeros(dimensionality)
    for ii in eachindex(walkerposition)
        walkerposition[ii] = 2*round(random_initialposition[ii]/2)
    end

    return walkerposition
end

function throwdice(random_number_bw_0_and_1::Float64)

    direction::String = "any_direction"

    if random_number_bw_0_and_1 <= 1/6
        direction = "positive_x"
    elseif random_number_bw_0_and_1 > 1/6 && random_number_bw_0_and_1 <= 1/3
        direction = "negative_x"
    elseif random_number_bw_0_and_1 > 1/3 && random_number_bw_0_and_1 <= 1/2
        direction = "positive_y"
    elseif random_number_bw_0_and_1 > 1/2 && random_number_bw_0_and_1 <= 2/3
        direction = "negative_y"
    elseif random_number_bw_0_and_1 > 2/3 && random_number_bw_0_and_1 <= 5/6
        direction = "positive_z"
    else
        direction = "negative_z"
    end

    return direction
end

function stepvector(direction::String,steplength::Integer,dimensionality::Integer)
    
    step::Vector{Integer} = zeros(dimensionality)

    if direction == "positive_x"
        step = [steplength,0,0]
    elseif direction == "negative_x"
        step = [-steplength,0,0]
    elseif direction == "positive_y"
        step = [0,steplength,0]
    elseif direction == "negative_y"
        step = [0,-steplength,0]
    elseif direction == "positive_z"
        step = [0,0,steplength]
    else
        step = [0,0,-steplength]
    end

    return step
end

function takestep!(walkerposition::Vector{Integer},step::Vector{Integer})
    
    walkerposition += step;

    return walkerposition
end

end