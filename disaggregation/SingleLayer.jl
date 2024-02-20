module SingleLayer

using LinearAlgebra

export build_singlelayermatrix

function maporientation_tointeger(face_orientation::String)
    # convert orientation to integer to build matrix (keeping the format from MATLAB's legacy code)
    normal_direction::Integer = 0
    if face_orientation == "positive_x" || face_orientation == "negative_x"
        normal_direction = 1
        elseif face_orientation == "positive_y" || face_orientation == "negative_y"
        normal_direction = 2
        elseif face_orientation == "positive_z" || face_orientation == "negative_z"
            normal_direction = 3
        else
            error("NOT a valid orientation")
    end
    
    return normal_direction
end

function compute_offdiagonalconstantterms(x_s::Float64,y_s::Float64,z_s::Float64,z::Float64,cij::Float64)
    if round(z) == round(z_s) && round(abs(x_s)) == 1 && round(abs(y_s)) == 1
        cij = 4*asinh(1)
    else
        if round(abs(x_s)) != 1 && round(abs(y_s)) != 1
            if round(z) == round(z_s)
                p1 = (1-y_s) * log(sqrt((1-x_s)*(1-x_s) + (1-y_s)*(1-y_s))   + (1-x_s)) + 
                     (1-x_s) * log(sqrt((1-x_s)*(1-x_s) + (1-y_s)*(1-y_s))   + (1-y_s)) -  (1-y_s)

                p2 = (-1-y_s)* log(sqrt((1-x_s)*(1-x_s) + (-1-y_s)*(-1-y_s)) + (1-x_s)) +
                     (1-x_s) * log(sqrt((1-x_s)*(1-x_s) + (-1-y_s)*(-1-y_s)) + (-1-y_s)) - (-1-y_s)

                p3 = (1-y_s) * log(sqrt((-1-x_s)*(-1-x_s) + (1-y_s)*(1-y_s)) + (-1-x_s)) + 
                     (-1-x_s)* log(sqrt((-1-x_s)*(-1-x_s) + (1-y_s)*(1-y_s)) + (1-y_s)) - (1-y_s)

                p4 = (-1-y_s)* log(sqrt((-1-x_s)*(-1-x_s) + (-1-y_s)*(-1-y_s)) + (-1-x_s)) + 
                     (-1-x_s)* log(sqrt((-1-x_s)*(-1-x_s) + (-1-y_s)*(-1-y_s)) + (-1-y_s)) - (-1-y_s)
            else
                p1 = - (z-z_s) * atan(((1-x_s)*(1-y_s))/((z-z_s)*sqrt(((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s)) + (1-y_s)*(1-y_s)))) + 
                       (z-z_s) * atan((1-y_s)/(z-z_s)) + (1-y_s) * log((1-x_s)+ sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s)+(1-y_s)*(1-y_s))) + 
                       (1-x_s) * log(sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s) + (1-y_s)*(1-y_s)) + (1-y_s)) - (1-y_s)

                p2 = - (z-z_s) * atan(((1-x_s)*(-1-y_s))/((z-z_s)*sqrt(((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s)) + (-1-y_s)*(-1-y_s)))) + 
                       (z-z_s) * atan((-1-y_s)/(z-z_s)) + (-1-y_s)* log((1-x_s) + sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s) + (-1-y_s)*(-1-y_s))) + 
                       (1-x_s) * log(sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s) + (-1-y_s)*(-1-y_s)) + (-1-y_s)) - (-1-y_s)

                p3 = - (z-z_s) * atan(((-1-x_s)*(1-y_s))/((z-z_s)*sqrt(((z-z_s)*(z-z_s) + (-1-x_s)*(-1-x_s)) + (1-y_s)*(1-y_s)))) + 
                       (z-z_s) * atan((1-y_s)/(z-z_s)) + (1-y_s) * log((-1-x_s) + sqrt((z-z_s)*(z-z_s) + (-1-x_s)*(-1-x_s) + (1-y_s)*(1-y_s))) + 
                       (-1-x_s)* log(sqrt((z-z_s)*(z-z_s) + (-1-x_s)*(-1-x_s) + (1-y_s)*(1-y_s)) + (1-y_s)) - (1-y_s)

                p4 = - (z-z_s) * atan(((-1-x_s)*(-1-y_s))/((z-z_s)*sqrt(((z-z_s)*(z-z_s) + (-1-x_s)*(-1-x_s)) + (-1-y_s)*(-1-y_s)))) + 
                       (z-z_s) * atan((-1-y_s)/(z-z_s)) + (-1-y_s)* log((-1-x_s) + sqrt((z-z_s)*(z-z_s) + (-1-x_s)*(-1-x_s) + (-1-y_s)*(-1-y_s))) + 
                       (-1-x_s)* log(sqrt((z-z_s)*(z-z_s) + (-1-x_s)*(-1-x_s) + (-1-y_s)*(-1-y_s)) + (-1-y_s)) - (-1-y_s)        
            end
        else
            if round(z) != round(z_s)
                p1 = -(z-z_s) *atan(((1-x_s)*(1-y_s))/((z-z_s)*sqrt(((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s)) + (1-y_s)*(1-y_s)))) + 
                       (z-z_s) *atan((1-y_s)/(z-z_s)) + (1-y_s)*log((1-x_s) + sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s) + (1-y_s)*(1-y_s))) + 
                       (1-x_s) *log(sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s) + (1-y_s)*(1-y_s)) + (1-y_s)) - (1-y_s)

                p2 = -(z-z_s) *atan(((1-x_s)*(-1-y_s))/((z-z_s)*sqrt(((z-z_s)*(z-z_s)  + (1-x_s)*(1-x_s)) + (-1-y_s)*(-1-y_s)))) + 
                       (z-z_s) *atan((-1-y_s)/(z-z_s)) + (-1-y_s)*log((1-x_s) + sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s) + (-1-y_s)*(-1-y_s))) + 
                       (1-x_s) *log(sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s) + (-1-y_s)*(-1-y_s)) + (-1-y_s)) - (-1-y_s)

                p3 = -(z-z_s) *atan(((-1-x_s)*(1-y_s))/((z-z_s)*sqrt(((z-z_s)*(z-z_s)  + (-1-x_s)*(-1-x_s)) + (1-y_s)*(1-y_s)))) + 
                       (z-z_s) *atan((1-y_s)/(z-z_s)) + (1-y_s)*log((-1-x_s) + sqrt((z-z_s)*(z-z_s) + (-1-x_s)*(-1-x_s) + (1-y_s)*(1-y_s))) + 
                       (-1-x_s)*log(sqrt((z-z_s)*(z-z_s) + (-1-x_s)*(-1-x_s) + (1-y_s)*(1-y_s)) + (1-y_s)) - (1-y_s)

                p4 = -(z-z_s) *atan(((-1-x_s)*(-1-y_s))/((z-z_s)*sqrt(((z-z_s)*(z-z_s) + (-1-x_s)*(-1-x_s)) + (-1-y_s)*(-1-y_s)))) + 
                       (z-z_s) *atan((-1-y_s)/(z-z_s)) + (-1-y_s)*log((-1-x_s) + sqrt((z-z_s)*(z-z_s) + (-1-x_s)*(-1-x_s) + (-1-y_s)*(-1-y_s))) + 
                       (-1-x_s)*log(sqrt((z-z_s)*(z-z_s) + (-1-x_s)*(-1-x_s) + (-1-y_s)*(-1-y_s)) + (-1-y_s)) - (-1-y_s)
            else # z = z_s
                if round(abs(x_s)) == 1 && round(abs(y_s)) != 1 
                    p1 = -((1/2)*(1-y_s)*(log((1-y_s)*(1-y_s)) - 2))
                    p2 = -((1-y_s)*(log(sqrt((1-y_s)*(1-y_s) + 4) + 2) - 1) + 2*asinh((1-y_s)/2))
                    p3 = -((1/2) *(-1-y_s)*(log((-1-y_s)*(-1-y_s)) - 2))
                    p4 = -((-1-y_s)*(log(sqrt((-1-y_s)*(-1-y_s) + 4) + 2) - 1) + 2*asinh((-1-y_s)/2))
                end
                if round(abs(x_s)) != 1 && round(abs(y_s)) == 1
                    r_s::Float64 = x_s; # swap them
                    x_s = y_s;
                    y_s = r_s;
                    
                    p1 = -((1/2)*(1-y_s)*(log((1-y_s)*(1-y_s)) - 2));
                    p2 = -((1-y_s)*(log(sqrt((1-y_s)*(1-y_s) + 4) + 2) - 1) + 2*asinh((1-y_s)/2));
                    p3 = -((1/2) *(-1-y_s)*(log((-1-y_s)*(-1-y_s)) - 2));
                    p4 = -((-1-y_s)*(log(sqrt((-1-y_s)*(-1-y_s) + 4) + 2) - 1) + 2*asinh((-1-y_s)/2));
                    
                    r_s = x_s; # unswap them
                    x_s = y_s;
                    y_s = r_s;
                end                         
            end   
        end
        cij = p1 - p2 - p3 + p4
    end

    return cij
end

function compute_offdiagonalxterms(x_s::Float64,y_s::Float64,z_s::Float64,z::Float64,offdiagonal_xij::Float64)
    if  round(z) != round(z_s) 
        if round(x_s) != 1 
            px1 = (1-x_s)*(atan(((1-y_s)*(1-x_s))/((z-z_s)*sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s) + (1-y_s)*(1-y_s)))))/((z-z_s)*(1-x_s))
            px2 = (1-x_s)*(atan(((1+y_s)*(1-x_s))/((z-z_s)*sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s) + (1+y_s)*(1+y_s)))))/((z-z_s)*(1-x_s))    
        else 
            px1 = 0
            px2 = 0
        end
        if round(x_s) != (-1)
            px3 = (1+x_s)*(atan(((1-y_s)*(1+x_s))/((z-z_s)*sqrt((z-z_s)*(z-z_s) + (1+x_s)*(1+x_s) + (1-y_s)*(1-y_s)))))/((z-z_s)*(1+x_s))
            px4 = (1+x_s)*(atan(((1+y_s)*(1+x_s))/((z-z_s)*sqrt((z-z_s)*(z-z_s) + (1+x_s)*(1+x_s) + (1+y_s)*(1+y_s)))))/((z-z_s)*(1+x_s))
        else
            px3 = 0
            px4 = 0
        end
    else
        px1 = 0
        px2 = 0
        px3 = 0
        px4 = 0
    end
    offdiagonal_xij = (z-z_s)*(z-z_s)*(px1+px2+px3+px4)

    return offdiagonal_xij
end

# Note to self: this function is way to big and untestable. Should be broken down into multiple functions ASAP

function build_singlelayermatrix(integration_point::Vector{Integer},face_center::Vector{Integer},face_orientation::String,dimensionality::Integer)
    position::Vector{Integer} = integration_point - face_center
    normal_direction::Integer = maporientation_tointeger(face_orientation)

    # main loop to build the entries of the matrix
    constant_terms::Matrix{Float64} = zeros(dimensionality,dimensionality)
    x_terms::Matrix{Float64} = zeros(dimensionality,dimensionality)
    for ii = 1:dimensionality
        for jj = 1:dimensionality
            if round(norm(position)) == 0 # at the singularity
                constant_terms[ii,jj] = (ii==jj) ? 8*asinh(1) : 0 # diagonal terms : non_diagonal terms
                x_terms[ii,jj] = (ii==jj && ii!=normal_direction) ? 4*asinh(1) : 0
            else # not at the singularity
                if ii == jj # on the diagonal
                    xs_index::Integer = mod(normal_direction+1,3)
                    ys_index::Integer = mod(normal_direction+2,3)
                    xs_index = (xs_index==0) ? 3 : xs_index
                    ys_index = (ys_index==0) ? 3 : ys_index

                    x_s::Float64 = position[xs_index,1]
                    y_s::Float64 = position[ys_index,1]
                    z_s::Float64 = position[normal_direction,1]
                    z::Float64  = 0
                    # assign the values to the constant terms first
                    offdiagonal_cij::Float64 = 0;
                    constant_terms[ii,jj] = compute_offdiagonalconstantterms(x_s,y_s,z_s,z,offdiagonal_cij)
                    #######################################################
                    # Now we assign the analytical values to the xx terms
                    if ii == normal_direction
                        offdiagonal_xij::Float64 = 0;
                        x_terms[ii,jj] = compute_offdiagonalxterms(x_s,y_s,z_s,z,offdiagonal_xij)
                    end
                else # ii!=jj, and so we are off the diagonal
                    constant_terms[ii,jj] = 0
                end
            end
        end
    end
  
    return x_terms
end

end