using PEBSI,Test,LinearAlgebra,PyCall,QHull,Distributions

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

#function rand_sc()
#    a=rand(Uniform(.1,5))
#    return [a 0 0;
#            0 a 0;
#            0 0 a]
#end
#function rand_fcc()
#    a=rand(Uniform(.1,5))
#    return [0 a/2 a/2;
#            a/2 0 a/2;
#            a/2 a/2 0]
#end
#function rand_bcc()
#    a=rand(Uniform(.1,5))
#    return [-a/2 a/2 a/2;
#            a/2 -a/2 a/2;
#            a/2 a/2 -a/2]
#end
#function rand_hex()
#    a=rand(Uniform(.1,5))
#    c=rand(Uniform(.1,5))
#    #I don't know why this is necessary, but phenum lib hangs without this.
#    if a < c
#        temp = a
#        a = c
#        c = temp
#    end
#        
#    return transpose([a/2 -(a*sqrt(3))/2 0;
#            a/2 (a*sqrt(3))/2 0;
#            0 0 c])
#end
#function rand_rhom_a()
#    a = rand(Uniform(.1,3))
#    α = rand(Uniform(π/10,π/2))
#    return transpose([a*cos(α/2) -a*sin(α/2) 0;
#            a*cos(α/2) a*sin(α/2) 0;
#            a*cos(α)/cos(α/2) 0 a*sqrt(1-(cos(α)^2)/(cos(α/2)^2))])
#end
#function rand_rhom_b()
#    a = rand(Uniform(.1,5))
#    # goes to 2 so that I don't have a negative root
#    α = rand(Uniform(π/2,2))
#    return transpose([a*cos(α/2) -a*sin(α/2) 0;
#            a*cos(α/2) a*sin(α/2) 0;
#            a*cos(α)/cos(α/2) 0 a*sqrt(1-(cos(α)^2)/(cos(α/2)^2))])
#end
#function rand_st()
#    α = rand(Uniform(π/2,2))
#    return transpose([a*cos(α/2) -a*sin(α/2) 0;
#            a*cos(α/2) a*sin(α/2) 0;
#            a*cos(α)/cos(α/2) 0 a*sqrt(1-(cos(α)^2)/(cos(α/2)^2))])
#end
#labels = ["sc" "fcc" "bcc" "hex" "rhom" "st" "bct" "so" "baseco" "bco" "fco" "sm" "basecm" "tric"]
##labels = ["sc"]
#sym = pyimport("bzip.symmetry")
#sc = Matrix{Float64}(I,3,3)
#fcc = [0 .5 .5;
#      .5 0 .5;
#      .5 .5 0]
#bcc = [-.5 .5 .5;
#       .5 -.5 .5;
#       .5 .5 -.5]
#hex =  sym.make_lattice_vectors("hexagonal",[1, 1, 2],[pi/2, pi/2, 2*pi/3])
#rhom = sym.make_lattice_vectors("rhombohedral",[1, 1, 1],[pi/4, pi/4, pi/4])
#st =   sym.make_lattice_vectors("tetragonal",[1, 1, 2],[pi/2, pi/2, pi/2])
#bct =  sym.make_lattice_vectors("body-centered tetragonal",[1, 1, 2],[pi/2, pi/2, pi/2])
#so =   sym.make_lattice_vectors("orthorhombic",[1, 2, 3],[pi/2, pi/2, pi/2])
#baseco = sym.make_lattice_vectors("base-centered orthorhombic",[1, 2, 3],[pi/2, pi/2, pi/2])
#bco =  sym.make_lattice_vectors("body-centered orthorhombic",[1, 2, 3],[pi/2, pi/2, pi/2])
#fco =  sym.make_lattice_vectors("face-centered orthorhombic",[1, 2, 3],[pi/2, pi/2, pi/2])
#sm =   sym.make_lattice_vectors("monoclinic",[1, 2, 3],[pi/4, pi/2, pi/2])
#basecm = sym.make_lattice_vectors("base-centered monoclinic",[1, 2, 3],[pi/4, pi/2, pi/2])
#tric =  sym.make_lattice_vectors("triclinic",[1, 2, 3],[pi/13, pi/7, pi/5])
#
#
#
#lats = [sc, fcc, bcc, hex, rhom, st, bct, so, baseco, bco, fco, sm, basecm, tric]
#lats = [sc]
lats = lattices
labels =lattice_labels

#commented only for testing! TODO
#@testset "make_bz_test" begin
#    i=1
#    @testset "$l" for l in labels
#        lat =  lats[i]
#        at = [1]
#        atom_pos = [0 0 0]
#        make_IBZ = pyimport("bzip.make_IBZ")
#        ch = make_IBZ.find_bz(lat)
#        #convert to julia
#        ch = chull(ch.points)
#        chjulia = chull(make_bz(lat))
#        @test ch.volume ≈  chjulia.volume
#        @test comparePolygon(ch, chjulia)
#        i+=1
#    end
#end


@testset "bz reduce fixed" begin
    i = 1
    @testset "$l" for l in labels
        println(labels[i])
        lat =  lats[i]
        println(lat)
        at = [1]
        atom_pos = [0 0 0]
        symmetry = pyimport("phenum.symmetry")

        ch = make_bz(lat)
        #convert to format
        ch = chull(ch)
       # origionalPoints = ch.points
       #need to transpose for phenum lib to give the correct space group
        spaceGroupTranspose = symmetry.get_spaceGroup(transpose(lat),at,atom_pos)[1]
        #spaceGroup = symmetry.get_spaceGroup(lat,at,atom_pos)[1]
        rch = reduce_bz(lat,at,atom_pos)
        rch_old = reduce_bz_old(lat,at,atom_pos)
        @test expectedOrder[i] ≈ size(spaceGroupTranspose)[1]
        @test rch.volume/ch.volume ≈ 1/size(spaceGroupTranspose)[1]
        @test comparePolygon(ch, unfold(rch,spaceGroupTranspose))
        @test BZMaps(ch,spaceGroupTranspose)

        @test rch_old.volume/ch.volume ≈ 1/size(spaceGroupTranspose)[1]
        @test comparePolygon(ch, unfold(rch_old,spaceGroupTranspose))
        i += 1
    end
end

@testset "bz reduce rand" begin
    lats = Array[]
    test_size = 2
    for i in 1:test_size push!(lats,rand_sc()) end
    for i in 1:test_size push!(lats,rand_fcc()) end
    for i in 1:test_size push!(lats,rand_bcc()) end
    for i in 1:test_size push!(lats,rand_hex()) end
    for i in 1:test_size push!(lats,rand_rhom_a()) end
    for i in 1:test_size push!(lats,rand_rhom_b()) end
    for i in 1:test_size push!(lats,rand_st()) end
    for i in 1:test_size push!(lats,rand_bct_a()) end
    for i in 1:test_size push!(lats,rand_bct_b()) end
    @testset "$i" for (i,lat) in enumerate(lats)
        println(lat)
        at = [1]
        atom_pos = [0 0 0]
        symmetry = pyimport("phenum.symmetry")

        ch = make_bz(lat)
        #convert to format
        #ch = copy(hcat(ch...)')
        ch = chull(ch)

       # origionalPoints = ch.points
       #need to transpose for phenum lib to give the correct space group
        spaceGroupTranspose = symmetry.get_spaceGroup(transpose(lat),at,atom_pos)[1]
        #println(size(spaceGroupTranspose)[1])
        #println(size(spaceGroupTranspose)[1])
        rch = reduce_bz(lat,at,atom_pos)
        #rch_old = reduce_bz_old(lat,at,atom_pos)
        @test rch.volume/ch.volume ≈ 1/size(spaceGroupTranspose)[1]
        @test comparePolygon(ch, unfold(rch,spaceGroupTranspose))
        @test BZMaps(ch,spaceGroupTranspose)

        #@test rch_old.volume/ch.volume ≈ 1/size(spaceGroupTranspose)[1]
        #@test comparePolygon(ch, unfold(rch_old,spaceGroupTranspose))
    end
end
