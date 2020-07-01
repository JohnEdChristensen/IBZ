export unfold, BZMaps,comparePolygon

function unfold(IBZ,spaceGroup)
    ops = [spaceGroup[i,:,:] for i=1:size(spaceGroup)[1]]
    points = IBZ.points
    for m in ops
        #println(m)
        newPoints = Array{Float64,2}(undef,0,3)

        for i in 1:size(points)[1]
            #println(points[i,:])
            newVert =  transpose(points[i,:]) * m
            #println(newVert)
            newPoints = vcat(newPoints,newVert)


        end
        points = vcat(points, newPoints)
        points = reducePoints(points)

    end
    return chull(points)
end
function BZMaps(BZ,spaceGroup)
    ops = [spaceGroup[i,:,:] for i=1:size(spaceGroup)[1]]
    for m in ops
        for i in 1:size(BZ.points)[1]
            #newVert =  transpose(transpose(BZ.points[i,:]) * m)
            newVert =  m*BZ.points[i,:]
            #check if the newVert maps to another point
            point_maps = false
            for j in 1:size(BZ.points)[1]
                if newVert ≈ BZ.points[j,:]
                    point_maps = true
                end
            end
            if !point_maps
                return false
            end
        end
    end
    return true
end
function comparePolygon(bz1,bz2)
    num_points_in_2 = 0
    for i in bz1.vertices
        point = bz1.points[i,:]
        point_in_2 = false
        for j in bz2.vertices
            compare_point = bz2.points[j,:]
            if compare_point ≈ point
                point_in_2 = true
                num_points_in_2 +=1
            end
        end
        if !point_in_2
            return false
        end
    end
    return num_points_in_2 == size(bz2.vertices)[1]
end
