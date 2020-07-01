using Test,IBZ,PyCall
include("testHelperFunctions.jl")
lats = lattices
labels =lattice_labels

@testset "bz reduce rand" begin
    convention = "ordinary"
    lats = Array[]
    test_size =20
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

        ch = make_bz(lat,convention)
        #convert to format
        #ch = copy(hcat(ch...)')
        ch = chull(ch)

       # origionalPoints = ch.points
       #need to transpose for phenum lib to give the correct space group
        spaceGroupTranspose = symmetry.get_spaceGroup(transpose(lat),at,atom_pos)[1]
        #println(size(spaceGroupTranspose)[1])
        #println(size(spaceGroupTranspose)[1])

        rch = reduce_bz(lat,at,atom_pos,convention)
        #rch_old = reduce_bz_old(lat,at,atom_pos)
        @test rch.volume/ch.volume ≈ 1/size(spaceGroupTranspose)[1]
        @test comparePolygon(ch, unfold(rch,spaceGroupTranspose))
        @test BZMaps(ch,spaceGroupTranspose)

        #@test rch_old.volume/ch.volume ≈ 1/size(spaceGroupTranspose)[1]
        #@test comparePolygon(ch, unfold(rch_old,spaceGroupTranspose))
    end
end
