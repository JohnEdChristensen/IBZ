using Test,IBZ,PyCall,QHull
include("testHelperFunctions.jl")
lats = lattices
labels =lattice_labels

@testset "bz reduce fixed" begin
    i = 0
    @testset "$l" for l in labels
        symmetry = pyimport("phenum.symmetry")
        at = [1]
        atom_pos = [0 0 0]
        convention = "ordinary"
        i+=1
        println(labels[i])
        lat =  lats[i]
        println(lat)
        println("Getting BZ...")
        ch = make_bz(lat,convention)
        #convert to format
        ch = chull(ch)
        # origionalPoints = ch.points
        #need to transpose for phenum lib to give the correct space group
        println("Getting symmetry...")
        spaceGroupTranspose = symmetry.get_spaceGroup(transpose(lat),at,atom_pos)[1]
        #spaceGroup = symmetry.get_spaceGroup(lat,at,atom_pos)[1]

        println("Getting IBZ...")
        rch = reduce_bz(lat,at,atom_pos,convention)
        #rch_old = reduce_bz_old(lat,at,atom_pos)

        @test expectedOrder[i] ≈ size(spaceGroupTranspose,1)

        @test ch.volume/rch.volume ≈ size(spaceGroupTranspose,1)
        println("Testing unflold...")
        @test comparePolygon(ch, unfold(rch,spaceGroupTranspose))
        println("Testing bz mapping...")
        @test BZMaps(ch,spaceGroupTranspose)

        #@test rch_old.volume/ch.volume ≈ 1/size(spaceGroupTranspose)[1]
        #@test comparePolygon(ch, unfold(rch_old,spaceGroupTranspose))
    end
end
