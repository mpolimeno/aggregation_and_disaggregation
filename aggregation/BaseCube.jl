module BaseCube

export positionoffaces_of_basecube
export orientationoffaces_ofbasecube

function positionoffaces_of_basecube(numberoffaces_inacube::Integer,dimensionality::Integer,centerof_baseface::Integer)
    if centerof_baseface != 1
        error("centerof_baseface MUST be 1");
    end
    # Initialize the position of the faces
    position_offaces::Matrix{Integer} = zeros(numberoffaces_inacube,dimensionality)

    top_index::Integer = 1
    bottom_index::Integer = numberoffaces_inacube
    while top_index < bottom_index
        position_offaces[top_index,top_index] = centerof_baseface
        position_offaces[bottom_index,top_index] = -1*centerof_baseface
        top_index += 1
        bottom_index -= 1
    end
    
    return position_offaces
end

function orientationoffaces_ofbasecube(faceindex::Integer)
    orientation::String = "any orientation"

    if faceindex == 1
        orientation = "positive_x"
    elseif faceindex == 2
        orientation = "positive_y"
    elseif faceindex == 3
        orientation = "positive_z"
    elseif faceindex == 4
        orientation = "negative_z"
    elseif faceindex == 5
        orientation = "negative_y"
    elseif faceindex == 6
        orientation = "negative_x"
    else
        error("faceindex MUST be integer between 1 and 6")
    end

    return orientation

end

end