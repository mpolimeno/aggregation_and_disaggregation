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

function compute_constantterms(x_s::Float64,y_s::Float64,z_s::Float64,z::Float64,cij::Float64)
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

function computexterms_caseone(x_s::Float64,y_s::Float64,z_s::Float64,z::Float64,xij_caseone::Float64)
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
    xij_caseone = (z-z_s)*(z-z_s)*(px1+px2+px3+px4)

    return xij_caseone
end

function computexterms_casetwo(x_s::Float64,y_s::Float64,z_s::Float64,z::Float64,xij_casetwo::Float64)
    if round(z) != round(z_s)
        px1 = - (z-z_s)*atan(((1-x_s)*(1-y_s))/((z-z_s)*sqrt((z-z_s)*(z-z_s) +(1-x_s)*(1-x_s) + (1-y_s)*(1-y_s)))) + 
                (z-z_s)*atan((1-y_s)/(z-z_s)) + (1-y_s)*(log(sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s) + (1-y_s)*(1-y_s)) + (1-x_s)) - 1)

        px2 = - (z-z_s)*atan(((1-x_s)*(-1-y_s))/((z-z_s)*sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s) + (-1-y_s)*(-1-y_s)))) + 
                (z-z_s)*atan((-1-y_s)/(z-z_s)) + (-1-y_s)*(log(sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s) + (-1-y_s)*(-1-y_s)) + (1-x_s)) - 1)

        px3 = - (z-z_s)*atan(((-1-x_s)*(1-y_s))/((z-z_s)*sqrt((z-z_s)*(z-z_s) + (-1-x_s)*(-1-x_s) + (1-y_s)*(1-y_s)))) + 
                (z-z_s)*atan((1-y_s)/(z-z_s)) + (1-y_s)*(log(sqrt((z-z_s)*(z-z_s) + (-1-x_s)*(-1-x_s) + (1-y_s)*(1-y_s)) + (-1-x_s)) - 1)

        px4 = - (z-z_s)*atan(((-1-x_s)*(-1-y_s))/((z-z_s)*sqrt((z-z_s)*(z-z_s) + (-1-x_s)*(-1-x_s) + (-1-y_s)*(-1-y_s)))) + 
                (z-z_s)*atan((-1-y_s)/(z-z_s)) + (-1-y_s)*(log(sqrt((z-z_s)*(z-z_s) + (-1-x_s)*(-1-x_s) + (-1-y_s)*(-1-y_s)) + (-1-x_s)) - 1)
    else
        if round(x_s) ==1
            if round(y_s) != 1
                px1 = 0.5*(1-y_s)*(log((1-y_s)*(1-y_s)) - 2)
            else
                px1 = 0
            end
            if round(y_s) != (-1)
                px2 = 0.5*(-1-y_s)*(log((-1-y_s)*(-1-y_s)) - 2)
            else
                px2 = 0
            end
        else
            if round(y_s) != 1
                px1 = (1-y_s)*(log(sqrt( (1-x_s)*(1-x_s) + (1-y_s)*(1-y_s)) + (1-x_s)) - 1)
            else
                px1 = 0
            end
            if round(y_s) != (-1)
                px2 =  (-1-y_s)*(log(sqrt((1-x_s)*(1-x_s) + (-1-y_s)*(-1-y_s)) + (1-x_s)) - 1)
            else
                px2 = 0
            end
        end
        if round(x_s) == (-1)
            if round(y_s) != 1
                px3 = (1-y_s)*(log(sqrt((-1-x_s)*(-1-x_s) + (1-y_s)*(1-y_s)) + (-1-x_s)) - 1)
            else
                px3 = 0
            end
            if round(y_s) != (-1)
                px4 = (-1-y_s)*(log(sqrt((-1-x_s)*(-1-x_s) + (-1-y_s)*(-1-y_s)) + (-1-x_s)) - 1)
            else
                px4 = 0
            end
        else
            if round(y_s) != 1
                px3 = (1-y_s)*(log(sqrt((-1-x_s)*(-1-x_s) + (1-y_s)*(1-y_s)) + (-1-x_s)) - 1)
            else
                px3 = 0
            end
            if round(y_s) != (-1)
                px4 = (-1-y_s)*(log(sqrt((-1-x_s)*(-1-x_s) + (-1-y_s)*(-1-y_s)) + (-1-x_s)) - 1)
            else
                px4 = 0
            end
        end
    end
    xij_casetwo = px1 - px2 - px3 + px4

    return xij_casetwo
end

function computexterms_casethree(x_s::Float64,y_s::Float64,z_s::Float64,z::Float64,xij_casethree::Float64)
    if round(z) != round(z_s)
        px1 =  -log(sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s)  + (1-y_s)*(1-y_s)) + (1-y_s))
        px2 =  log(sqrt((z-z_s)*(z-z_s) + (1+x_s)*(1+x_s)  + (1-y_s)*(1-y_s)) + (1-y_s))
        px3 =  log(sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s)  + (-1-y_s)*(-1-y_s)) + (-1-y_s))
        px4 =  -log(sqrt((z-z_s)*(z-z_s) + (1+x_s)*(1+x_s)  + (-1-y_s)*(-1-y_s)) + (-1-y_s))

    else
        px1 = 0
        px2 = 0
        px3 = 0
        px4 = 0
    end
    xij_casethree = (z-z_s)*(px1 + px2 + px3 + px4)

    return xij_casethree
end

function computexterms_casefour(x_s::Float64,y_s::Float64,z_s::Float64,z::Float64,xij_casefour::Float64)
    if round(z) != round(z_s)
        px1 =  -log(sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s)  + (1-y_s)*(1-y_s)) + (1-y_s))
        px2 =  log(sqrt((z-z_s)*(z-z_s) + (1+x_s)*(1+x_s)  + (1-y_s)*(1-y_s)) + (1-y_s))
        px3 =  log(sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s)  + (-1-y_s)*(-1-y_s)) + (-1-y_s))
        px4 =  -log(sqrt((z-z_s)*(z-z_s) + (1+x_s)*(1+x_s)  + (-1-y_s)*(-1-y_s)) + (-1-y_s))
    else
        px1 = 0
        px2 = 0
        px3 = 0
        px4 = 0
    end
    xij_casefour= (z-z_s)*(px1 + px2 + px3 + px4)

    return xij_casefour
end

function computexterms_casefive(x_s::Float64,y_s::Float64,z_s::Float64,z::Float64,xij_casefive::Float64)
    px1 = -sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s) + (1-y_s)*(1-y_s))
    px2 = sqrt((z-z_s)*(z-z_s) + (1-x_s)*(1-x_s) + (1+y_s)*(1+y_s))
    px3 = sqrt((z-z_s)*(z-z_s) + (1+x_s)*(1+x_s) + (1-y_s)*(1-y_s))
    px4 = -sqrt((z-z_s)*(z-z_s) + (1+x_s)*(1+x_s) + (1+y_s)*(1+y_s))
    
    xij_casefive = px1 + px2 + px3 + px4
    return xij_casefive
end

function build_singlelayermatrix(integration_point::Vector{Integer},face_center::Vector{Integer},face_orientation::String,dimensionality::Integer)
    position::Vector{Integer} = integration_point - face_center
    normal_direction::Integer = maporientation_tointeger(face_orientation)
    numberoffaces_inacube::Integer = 6;

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
                    first_cij::Float64 = 0;
                    constant_terms[ii,jj] = compute_constantterms(x_s,y_s,z_s,z,first_cij)
                    #######################################################
                    # Now we assign the analytical values to the xx terms
                    if ii == normal_direction
                        # CASE 1
                        xij_caseone::Float64 = 0;
                        x_terms[ii,jj] = computexterms_caseone(x_s,y_s,z_s,z,xij_caseone)
                    else
                        # CASE 2
                        x_s2::Float64 = position[ii,1]
                        y_s2::Float64 = position[numberoffaces_inacube-normal_direction-ii,1]
                        z_s2::Float64 = position[normal_direction,1]
                        z2::Float64  = 0
                        xij_casetwo::Float64 = 0;
                        x_terms[ii,jj] = computexterms_casetwo(x_s2,y_s2,z_s2,z2,xij_casetwo)
                    end

                else # ii!=jj, and so we are off the diagonal
                    constant_terms[ii,jj] = 0
                    #######################################################
                    # Now we assign the analytical values to the xx terms
                    if ii != normal_direction && jj == normal_direction
                        # CASE 3
                        x_s3::Float64 = position[ii,1]
                        y_s3::Float64 = position[numberoffaces_inacube-normal_direction-ii,1]
                        z_s3::Float64 = position[normal_direction,1]
                        z3::Float64  = 0
                        xij_casethree::Float64 = 0;
                        x_terms[ii,jj] = computexterms_casethree(x_s3,y_s3,z_s3,z3,xij_casethree)
                    elseif ii == normal_direction && jj != normal_direction
                        # CASE 4
                        x_s4::Float64 = position[jj,1]
                        y_s4::Float64 = position[numberoffaces_inacube-normal_direction-jj,1]
                        z_s4::Float64 = position[normal_direction,1]
                        z4::Float64  = 0
                        xij_casefour::Float64 = 0;
                        x_terms[ii,jj] = computexterms_casefour(x_s4,y_s4,z_s4,z4,xij_casefour)
                    else
                        # CASE 5
                        x_s5::Float64 = position[ii,1]
                        y_s5::Float64 = position[jj,1]
                        z_s5::Float64 = position[normal_direction,1]
                        z5::Float64  = 0
                        xij_casefive::Float64 = 0;
                        x_terms[ii,jj] = computexterms_casefive(x_s5,y_s5,z_s5,z5,xij_casefive)
                    end
                end
            end
        end
    end
  
    return constant_terms, x_terms
end

end