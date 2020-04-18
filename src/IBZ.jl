module IBZ

using LinearAlgebra,PyCall,QHull,Polyhedra,CDDLib, Combinatorics,Distances, Distributions

export reduce_bz,reducePoints
export make_bz, make_bz_2d
export reduce_bz_old
export rand_sc,rand_fcc,rand_bcc,rand_hex,rand_rhom_a,rand_rhom_b,rand_st, rand_bct_a,rand_bct_b,rand_so
export lattices,lattice_labels,expectedOrder

@doc """
function projectToPlane(plane,point)

returns the projection of point onto plane
"""
function projectToPlane(plane,point)
    plane = normalize(plane)
    t = (-plane[1]*point[1] - plane[2]*point[2] -plane[3]*point[3])/
    (plane[1]^2+plane[2]^2+plane[3]^2)
    return [point[1] + t*plane[1] 
            point[2] + t*plane[2]
            point[3] + t*plane[3]]
end


@doc """
getIntercept(p1,p2,n)
returns a point between p1 and p2 if the plane normal to n intersects
return Nothing if it does not
"""
function getIntercept(p1,p2,n,pPlane=[0.0,0.0,0.0])
    l = p1-p2
    if dot(l,n) ≈ 0
        #lines are parallel do nothing
        return Nothing
    end
    d = dot((pPlane-p1),n)/dot(l,n)
    pIntersect = d*l + p1
    #check if pIntersect is inbetween p1 and p2
    dP1P2 = norm(p1-p2)
    dP1PIntersect = norm(p1-pIntersect)
    dP2PIntersect = norm(p2-pIntersect)
    if dP1P2 ≈ dP1PIntersect + dP2PIntersect
        return pIntersect
    end
    # doesn't intersect the line segment
    return Nothing
end

@doc """
returns the ch points plus the points of intersection from the plane
and edges of the convex hull
"""
function getIntersectPoints(plane,ch)

    oldPoints = ch.points
    newPoints = Array{Float64,2}(undef,0,3)

    #create vertices at plane intersect
    for s ∈ ch.simplices
        #test the 3 lines making up each simplex
        i1 = getIntercept(oldPoints[s[1],:],oldPoints[s[2],:],plane)
        #println(i1)
        if i1 != Nothing
            newPoints = vcat(newPoints,i1')
        end
        i2 = getIntercept(oldPoints[s[2],:],oldPoints[s[3],:],plane)
        if i2 != Nothing
            newPoints = vcat(newPoints,i2')
        end
        i3 = getIntercept(oldPoints[s[3],:],oldPoints[s[1],:],plane)
        if i3 != Nothing
            newPoints = vcat(newPoints,i3')
        end
    end
    return newPoints

end

@doc """
function reducePoints(points)
gets rid of duplicate points and points that lie along an already
defined edge. also rounds points near zero to zero
"""
function reducePoints(points)
    ch = chull(points)
    reducedPoints = Array{Float64,2}(undef,0,3)
    for i ∈ ch.vertices
        if isapprox(ch.points[i,:],[0;0;0],atol=eps(Float64))
            reducedPoints = vcat(reducedPoints,[0 0 0])
        else
            for j ∈ 1:3
                if isapprox(ch.points[i,j],0,atol=eps(Float64))
                    ch.points[i,j] = 0
                end
            end
            reducedPoints = vcat(reducedPoints,ch.points[i,:]')
        end
    end
    return reducedPoints
end



@doc """
reflectionReduce(rPlane,ch)

removes all points from the convex hull ch that are on the opposite side of
the plane defined by the vector rPlane. Creates additional verticies at 
plane intersections of the simplicies of the convex hull

"""
function reflectionReduce(rPlane,ch)
    #println("rPlane")
    #println(rPlane)

    oldPoints = ch.points

    newPoints = getIntersectPoints(rPlane,ch)
    #remove all points that have a negative dot product
    for i ∈ ch.vertices
        p = oldPoints[i,:]
        if isapprox(dot(p,rPlane), 0,atol=eps(Float64)) || dot(p,rPlane) > 0
            #println("is above plane")
            #println(p)
            newPoints = vcat(newPoints,p')
        end
    end
    #println(newPoints)
    try
        rPoints = reducePoints(newPoints)
        return chull(rPoints)
    catch
        println("attempted to reduced already reduced volume")
        return ch
    end
end

@doc """
reduce_bz(lat, at, atom_pos)

returns the reduced Brillouin zone

# Arguments
- `lat:Array{Float64,3,3}`: The lattice vectors
- `at:Array{Int64,n}`: The atom types
- `atom_pos:Array{Int64,n,3}`: position of the atoms

"""
function reduce_bz(lat, at, atom_pos)

    bz = make_bz(lat,false)

    verts =  collect(points(polyhedron(bz,CDDLib.Library())))

    obz = copy(hcat(verts...)')
    obz = chull(obz)

    symmetry = pyimport("phenum.symmetry")
    #transpose to use phenumlib
    spaceGroup = symmetry.get_spaceGroup(transpose(lat),at,atom_pos)[1]
    ops = [spaceGroup[i,:,:] for i=1:size(spaceGroup)[1]]
    #keep track of the indices of the stabilzer operators ( the operators that did not move a given point)
    stabilizer_index = []

    index = 1
    while size(ops)[1] != 1
        for i=1:(size(ops)[1])
            m = ops[i]

            vert = obz.points[obz.vertices[index],:]

            rotated_vert = m*vert

            # check if m is a stabilizer, if it is, mark it to be kept
            if (vert ≈ rotated_vert)
                append!(stabilizer_index,i)
            else
                cut_plane = (vert - rotated_vert)

                bz = bz ∩ HalfSpace(cut_plane,0)
            end
        end

        #if we only have stabilizers move to a new point
        if size(stabilizer_index)[1] == size(ops)[1]
            index = (index + 1)%size(obz.vertices)[1]
        else
            index = 1
        end
        #remove all but the stabilizers
        ops = [ops[i,:,:][1] for i in stabilizer_index]
        stabilizer_index= []
    end

    #convert to verts
    verts =  collect(points(polyhedron(bz,CDDLib.Library())))
    # convert to 2d array
    return chull(copy(hcat(verts...)'))

end #function

function reduce_bz_2d(lat, spaceGroup)

    bz = make_bz_2d(lat,false)

    verts =  collect(points(polyhedron(bz,CDDLib.Library())))

    obz = copy(hcat(verts...)')
    obz = chull(obz)

    #symmetry = pyimport("phenum.symmetry")
    #transpose to use phenumlib
    #spaceGroup = symmetry.get_spaceGroup(transpose(lat),at,atom_pos)[1]
    #ops = [spaceGroup[i,:,:] for i=1:size(spaceGroup)[1]]
    #keep track of the indices of the stabilzer operators ( the operators that did not move a given point)
    stabilizer_index = []

    index = 1
    while size(ops)[1] != 1
        for i=1:(size(ops)[1])
            m = ops[i]

            vert = obz.points[obz.vertices[index],:]

            rotated_vert = m*vert

            # check if m is a stabilizer, if it is, mark it to be kept
            if (vert ≈ rotated_vert)
                append!(stabilizer_index,i)
            else
                cut_plane = (vert - rotated_vert)

                bz = bz ∩ HalfSpace(cut_plane,0)
            end
        end

        #if we only have stabilizers move to a new point
        if size(stabilizer_index)[1] == size(ops)[1]
            index = (index + 1)%size(obz.vertices)[1]
        else
            index = 1
        end
        #remove all but the stabilizers
        ops = [ops[i,:,:][1] for i in stabilizer_index]
        stabilizer_index= []
    end

    #convert to verts
    verts =  collect(points(polyhedron(bz,CDDLib.Library())))
    # convert to 2d array
    return chull(copy(hcat(verts...)'))

end #function

@doc """
returns the brillouin zone for the given lattice,

by default returns the verticies the of the bz, 
can also return the Half space representation

!!! does not currently work for skewed cells, need to implement minkowski reduction !!!!!
"""
function make_bz(lat,vertsOrHrep = true)
    #@show vertsOrHrep   
    #enumarate all lattice points 2 on each side
    lattice_points = collect(Iterators.product(-2:2,-2:2,-2:2))
    #make a array of vectors
    lattice_points = [collect(i) for i in vec(lattice_points)]
    #calculate the Cartesian position of each lattice point
    cart_point = [lat*i for i in lattice_points]
    #sort based on distance from origin
    distance = [norm(i) for i in cart_point]
    cart_point = cart_point[sortperm(distance)]
    #create the first half space
    p1 =cart_point[1]/2
    d = p1[1]^2+p1[2]^2+p1[3]^2
    bz = HalfSpace(p1,d)
    for lp in cart_point[2:end]
        p = lp/2
        d = p[1]^2+p[2]^2+p[3]^2
        bz = bz ∩ HalfSpace(p,d)
    end
    if vertsOrHrep
        #return polyhedron(bz,CDDLib.Library())
        verts =  collect(points(polyhedron(bz,CDDLib.Library())))
        #2d array
        bz = copy(hcat(verts...)')
    end
    return bz
end

function make_bz_2d(lat,vertsOrHrep = true)
    #@show vertsOrHrep   
    #enumarate all lattice points 2 on each side
    lattice_points = collect(Iterators.product(-2:2,-2:2))
    #make a array of vectors
    lattice_points = [collect(i) for i in vec(lattice_points)]
    #calculate the Cartesian position of each lattice point
    cart_point = [lat*i for i in lattice_points]
    #sort based on distance from origin
    distance = [norm(i) for i in cart_point]
    cart_point = cart_point[sortperm(distance)]
    #create the first half space
    p1 =cart_point[1]/2
    d = p1[1]^2+p1[2]^2
    bz = HalfSpace(p1,d)
    for lp in cart_point[2:end]
        p = lp/2
        d = p[1]^2+p[2]^2
        bz = bz ∩ HalfSpace(p,d)
    end
    if vertsOrHrep
        #return polyhedron(bz,CDDLib.Library())
        verts =  collect(points(polyhedron(bz,CDDLib.Library())))
        #2d array
        bz = copy(hcat(verts...)')
        order = detriangulate(bz,collect(1:size(bz)[1]))
        print(order,:)
    end
    return bz[order,:]
end

@doc """
reduce_bz_old(lat, at, atom_pos)

returns the reduced Brillouin zone, uses a more a cutting plane to reduce the bz instead of 
halfspaces, used only for refrence



!!! this method fails some test cases, I think it is some sort of finite preciesion error in 
reflection reduce !!!!

# Arguments
- `lat:Array{Float64,3,3}`: The lattice vectors
- `at:Array{Int64,n}`: The atom types
- `atom_pos:Array{Int64,n,3}`: position of the atoms

"""
function reduce_bz_old(lat, at, atom_pos)

    bz = make_bz(lat)
    bz = chull(bz)
    obz = bz

    symmetry = pyimport("phenum.symmetry")
    #transpose to use phenumlib
    spaceGroup = symmetry.get_spaceGroup(transpose(lat),at,atom_pos)[1]
    ops = [spaceGroup[i,:,:] for i=1:size(spaceGroup)[1]]
    stabilizer_index = []

    index = 1
    while size(ops)[1] != 1
        for i=1:(size(ops)[1])

            m = ops[i]

            #create point, rotated point pair to determine cutting plane
            vert = obz.points[obz.vertices[index],:]
            rotated_vert = m*vert

            # check if m is a stabilizer, if it is, mark it to be kept
            if (vert ≈ rotated_vert)
                append!(stabilizer_index,i)
            else
                cut_plane = (vert - rotated_vert)
                bz = reflectionReduce(cut_plane,bz)
            end
        end
        #if we only have stabilzers, move to a new point
        if size(stabilizer_index)[1] == size(ops)[1]
            index = (index + 1)%size(obz.vertices)[1]
        else
            index = 1
        end

        #remove all but the stabilizers
        ops = [ops[i,:,:][1] for i in stabilizer_index]
        stabilizer_index= []
    end

    return bz

end #function

function write_obj(verts,faces,path)
    io = open(path,"w")
    for v in eachrow(verts)
        write(io, "v ") 
        for p in v
            write(io,string(p, " " ))
        end
        write(io,"\n")
    end
    for f in faces
        write(io, "f ") 
        for i in f
            write(io,string(i, " " ))
        end
        write(io,"\n")
    end

    close(io)
end

function simplify_faces(chull)
    vertices = chull.vertices
    points = chull.points
    faces = chull.simplices
    polygons = get_polygons(faces,points)
    for i in 1:size(polygons)[1]
        polygons[i] = detriangulate(points,polygons[i])
    end
    return polygons
end
function detriangulate(points,polygon)
    indices = unique!(vcat(polygon...))
    possible_orderings = collect(permutations(indices))
    siz = size(possible_orderings[1])[1]
    min_distance  = Inf
    min_order = 0
    for (order_num,order) in enumerate(possible_orderings)
        #@show order
        total = 0
        for (i,index) in enumerate(order)
            #wrap to the beginning
            next_point = i == siz ? order[1] : order[i+1]
            #distance
            total += euclidean(points[index,:],points[next_point,:])
        end
        if total < min_distance
            @show min_distance
            @show total
            min_distance = total
            min_order = order_num
            if total == 0
                @show points
                @show order
            end
        end

    end 
    return possible_orderings[min_order]
end
function get_polygons(faces,points)
    #each element  contains the indices of the vertices that make up a polygon
    polygons = []
    processed_faces = []
    for (findex, face) in enumerate(faces)
        #@show findex
        #@show processed_faces
        #check if the face has already been processed
        if indexin(findex,processed_faces)[1] === nothing
            push!(processed_faces,findex)
            polygon = []
            p0 = points[face[1],:] 
            v1 = points[face[1],:] - points[face[2],:]
            v2 = points[face[1],:] - points[face[3],:]
            normal = normalize(cross(v1,v2))
            #equation of a plane ax+by+cz = d
            a = normal[1] 
            b = normal[2] 
            c = normal[3] 
            d = a *p0[1] + b*p0[2] +c*p0[3]
            #find all faces whose points are all in this plane
            for (findex2,face2) in enumerate(faces)
                in_plane = true
                #@show face2
                for i in face2
                    p = points[i,:]
                    if !(a*p[1]+b*p[2]+c*p[3] ≈  d)
                        in_plane=false
                        break
                    end
                end
                if in_plane
                    push!(polygon,face2)
                    push!(processed_faces,findex2)
                end
            end
            push!(polygons,polygon)
        end
    end
    return polygons
end
#lattices used for testing
function order_params(a,b,c)
    if a > b
        a,b = b,a
    end
    if b > c
        b,c= c,b
    end
    if a > b
        a,b = b,a
    end
    return a,b,c
end
    
function rand_sc(a=rand(Uniform(.1,5)))
    return [a 0 0;
            0 a 0;
            0 0 a]
end
function rand_fcc(a=rand(Uniform(.1,5)))
    return [0 a/2 a/2;
            a/2 0 a/2;
            a/2 a/2 0]
end
function rand_bcc(a=rand(Uniform(.1,5)))
    return [-a/2 a/2 a/2;
            a/2 -a/2 a/2;
            a/2 a/2 -a/2]
end
function rand_hex(a=rand(Uniform(.1,5)),c=rand(Uniform(.1,5)))
    #I don't know why this is necessary, but phenum lib hangs without this.
    #TODO fix
    if a < c
        temp = a
        a = c
        c = temp
    end
    return transpose([a/2 -(a*sqrt(3))/2 0;
                      a/2 (a*sqrt(3))/2 0;
                      0 0 c])
end
function rand_rhom_a(a = rand(Uniform(.1,3)),α = rand(Uniform(π/10,π/2)))
    return transpose([a*cos(α/2) -a*sin(α/2) 0;
                      a*cos(α/2) a*sin(α/2) 0;
                      a*cos(α)/cos(α/2) 0 a*sqrt(1-(cos(α)^2)/(cos(α/2)^2))])
end
function rand_rhom_b(a = rand(Uniform(.1,5)),α = rand(Uniform(π/2,2)))
    # goes to 2 so that I don't have a negative root
    return transpose([a*cos(α/2) -a*sin(α/2) 0;
                      a*cos(α/2) a*sin(α/2) 0;
                      a*cos(α)/cos(α/2) 0 a*sqrt(1-(cos(α)^2)/(cos(α/2)^2))])
end
function rand_st(a = rand(Uniform(.1,5)),c = rand(Uniform(.1,5)))
    return transpose([a 0 0;
                      0 a 0;
                      0 0 c])
end

function rand_bct_a(a = rand(Uniform(.1,5)),c = rand(Uniform(.1,5)))
    if c > a
        #swap
        a,c = c,a
    end
    return transpose([-a/2 a/2 c/2;
                      a/2 -a/2 c/2;
                      a/2 a/2 -c/2])
end

function rand_bct_b(a = rand(Uniform(.1,5)),c = rand(Uniform(.1,5)))
    if c < a
        #swap
        a,c = c,a
    end
    return transpose([-a/2 a/2 c/2;
                      a/2 -a/2 c/2;
                      a/2 a/2 -c/2])
end

function rand_so(a = rand(Uniform(.1,5)),b = rand(Uniform(.1,5)),c = rand(Uniform(.1,5)))
    a,b,c = order_params(a,b,c)
    if a > b || b > c
        print("error in so conditions")
    end
    return transpose([a 0 0;
                      0 b 0;
                      0 0 c])
end

function rand_baseco(a = rand(Uniform(.1,5)),b = rand(Uniform(.1,5)),c = rand(Uniform(.1,5)))
    a,b,c = order_params(a,b,c)
    if a > b || b > c
        print("error in caseco so conditions")
    end
    return transpose([a/2 -b/2 0;
                      a/2 b/2 0;
                      0   0 c])
end

function rand_bco(a = rand(Uniform(.1,5)),b = rand(Uniform(.1,5)),c = rand(Uniform(.1,5)))
    return transpose([-a/2 b/2 c/2;
                      a/2 -b/2 c/2;
                      a/2 b/2 -c/2])
end

#1/a^2 > 1/b^2 + 1/c^2
#TODO figure out why we need certain ordering? (a < b < c)
function rand_fco_a(a = rand(Uniform(.1,5)),b = rand(Uniform(.1,5)),c = rand(Uniform(.1,5)))
    if a > b
        a,b = b,a
    end
    if 1/a^2 <= 1/b^2 + 1/c^2
        range = sqrt((1/a^2 - 1/b^2)^(-1))
        c = rand(Uniform(range,2*range))
        if 1/a^2 <= 1/b^2 + 1/c^2
            print("incorrect conditions for rand_fco_a")
        end
    end
    return transpose([0 b/2 c/2;
                      a/2 0 c/2;
                      a/2 b/2 0])
end

#1/a^2 < 1/b^2 + 1/c^2
function rand_fco_b(a = rand(Uniform(.1,5)),b = rand(Uniform(.1,5)),c = rand(Uniform(.1,5)))
    if a > b
        a,b = b,a
    end
    if 1/a^2 >= 1/b^2 + 1/c^2
        range = sqrt((1/a^2 - 1/b^2)^(-1))
        c = rand(Uniform(0,range))
        if 1/a^2 <= 1/b^2 + 1/c^2
            print("incorrect conditions for rand_fco_a")
        end
    end
    return transpose([0 b/2 c/2;
                      a/2 0 c/2;
                      a/2 b/2 0])
end

#1/a^2 = 1/b^2 + 1/c^2
function rand_fco_c(a = rand(Uniform(.1,5)),b = rand(Uniform(.1,5)),c = rand(Uniform(.1,5)))
    if a > b
        a,b = b,a
    end
    if 1/a^2 != 1/b^2 + 1/c^2
        c = sqrt((1/a^2 - 1/b^2)^(-1))
        if 1/a^2 != 1/b^2 + 1/c^2
            print("incorrect conditions for rand_fco_c")
        end
    end
    return transpose([0 b/2 c/2;
                      a/2 0 c/2;
                      a/2 b/2 0])
end

function rand_sm(a = rand(Uniform(.1,5)),b = rand(Uniform(.1,5)),c = rand(Uniform(.1,5)),α  = rand(Uniform(0,π/2)))
    return transpose([a 0 0;
                      0 b 0;
                      0 c*cos(α) c*sin(α)])
end
#need help here, what is  k_gamma?
function rand_basecm_a(a = rand(Uniform(.1,5)),b = rand(Uniform(.1,5)),c = rand(Uniform(.1,5)),α  = rand(Uniform(0,π/2)))
    return transpose([a/2 b/2 0;
                      -a/2 b/2 0;
                      0 c*cos(α) c*sin(α)])
end


lattice_labels = ["sc" "fcc" "bcc" "hex" "rhoma" "rhomb" "st" "bcta" "bctb" "so" "baseco" "bco" "fcoa" "fcob" "fcoc" "sm" "basecm" "tric"]
#labels = ["sc"]
sym = pyimport("bzip.symmetry")
sc = rand_sc(1)
fcc = rand_fcc(1)
bcc = rand_bcc(1)
hex =  rand_hex()
rhoma = rand_rhom_a()
rhomb = rand_rhom_b()
st = rand_st() 
bcta =  rand_bct_a()
bctb =  rand_bct_b()
so =   rand_so()
baseco = rand_baseco()
bco =  rand_bco()
fcoa = rand_fco_a()
fcob = rand_fco_b()
fcoc = rand_fco_c()
sm =  rand_sm()
basecm =  rand_basecm_a()
tric =  sym.make_lattice_vectors("triclinic",[1, 2, 3],[pi/13, pi/7, pi/5])

lattices =      [sc, fcc, bcc, hex, rhoma,rhomb, st, bcta, bctb, so, baseco, bco, fcoa, fcob, fcoc, sm, basecm, tric]
expectedOrder = [48, 48,  48,  24,  12,   12,    16, 16,   16,   8,  8,      8,   8,    8,    8,    4,  4,      2]
    end #module


