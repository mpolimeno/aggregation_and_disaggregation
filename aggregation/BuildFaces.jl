module BuildFaces

include("BaseCube.jl")
using .BaseCube

export build_externalfacesofaggregate
export getorientation_externalfacesofaggregate

function build_externalfacesofaggregate(numberoffaces_inacube::Integer,dimensionality::Integer,centerof_baseface::Integer,position_of_cubes::Matrix{Integer})
    # get the position of the faces in the base cube to build all the faces of the aggregate
    facesof_basecube::Matrix{Integer} = positionoffaces_of_basecube(numberoffaces_inacube,dimensionality,centerof_baseface)
    number_of_cubes::Integer = size(position_of_cubes,1);
    
    allfaces = [];
    faces_toremove = [];
    for cube in axes(position_of_cubes,1)
        # assign position to all the faces in all the cubes in the aggregate
        # this will have internal and external faces
        facesof_currentcube::Matrix{Integer} = zeros(numberoffaces_inacube,dimensionality)
        for face=1:numberoffaces_inacube
            facesof_currentcube[face,:] = vec(position_of_cubes[cube,:] + facesof_basecube[face,:])
            push!(allfaces,facesof_currentcube[face,:]);
        end
        # check for adjacent cubes, so that to remove faces later
        for jj=1:(cube-1)
            for kk=1:numberoffaces_inacube
                if isequal(vec(position_of_cubes[cube,:] + 2*facesof_basecube[kk,:]),position_of_cubes[jj,:])
                    push!(faces_toremove,[cube kk])
                    push!(faces_toremove,[jj 7-kk])
                end
            end
        end
    end
    # now we remove the useless faces
    numberoffaces_toremove::Integer = size(faces_toremove,1)
    facestoremove_array::Matrix{Integer} = zeros(size(faces_toremove,1),2)
    for ii in axes(facestoremove_array,1)
        facestoremove_array[ii,:] = (faces_toremove[ii]);
    end
    faces_tokeep = [];
    for jj=1:(number_of_cubes*numberoffaces_inacube)
        removeface::Bool = false;
        for ii=1:numberoffaces_toremove
            current_face::Integer = numberoffaces_inacube*(facestoremove_array[ii,1]-1) + facestoremove_array[ii,2]
            if jj == current_face
                removeface = true
            end
        end
        if removeface == false
            push!(faces_tokeep,jj)
        end
    end

    externalfaces::Matrix{Integer} = zeros(size(faces_tokeep,1),dimensionality)
    for ii in axes(externalfaces,1)
        externalfaces[ii,:] = allfaces[faces_tokeep[ii]]
    end

    return externalfaces
end


function getorientation_externalfacesofaggregate(numberoffaces_inacube::Integer,dimensionality::Integer,centerof_baseface::Integer,position_of_cubes::Matrix{Integer})
    number_of_cubes::Integer = size(position_of_cubes,1);

    # get the position of the faces in the base cube to build all the faces of the aggregate
    facesof_basecube::Matrix{Integer} = positionoffaces_of_basecube(numberoffaces_inacube,dimensionality,centerof_baseface)
    
    allorientations = [];
    faces_toremove = [];
    for cube in axes(position_of_cubes,1)
        for face=1:numberoffaces_inacube
            orientationof_currentface::String = orientationoffaces_ofbasecube(face)
            push!(allorientations,orientationof_currentface);
        end
        # check for adjacent cubes, so that to remove faces later
        for jj=1:(cube-1)
            for kk=1:numberoffaces_inacube
                if isequal(vec(position_of_cubes[cube,:] + 2*facesof_basecube[kk,:]),position_of_cubes[jj,:])
                    push!(faces_toremove,[cube kk])
                    push!(faces_toremove,[jj 7-kk])
                end
            end
        end
    end
    # now we remove the useless faces
    numberoffaces_toremove::Integer = size(faces_toremove,1)
    facestoremove_array::Matrix{Integer} = zeros(size(faces_toremove,1),2)
    for ii in axes(facestoremove_array,1)
        facestoremove_array[ii,:] = (faces_toremove[ii]);
    end
    faces_tokeep = [];
    for jj=1:(number_of_cubes*numberoffaces_inacube)
        removeface::Bool = false;
        for ii=1:numberoffaces_toremove
            current_face::Integer = numberoffaces_inacube*(facestoremove_array[ii,1]-1) + facestoremove_array[ii,2]
            if jj == current_face
                removeface = true
            end
        end
        if removeface == false
            push!(faces_tokeep,jj)
        end
    end
    
    orientationof_externalfaces::Vector{String} = [];
    for ii=1:size(faces_tokeep,1)
        orientationof_externalfaces = push!(orientationof_externalfaces,allorientations[faces_tokeep[ii]])
    end

    return orientationof_externalfaces
end

end