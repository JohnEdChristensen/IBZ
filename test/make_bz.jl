using Test,IBZ,PyCall
include("testHelperFunctions.jl")
lats = lattices
labels =lattice_labels

@testset "make_bz_test" begin
    symmetry = pyimport("phenum.symmetry")
    i=1
    at = [1]
    atom_pos = [0 0 0]
    convention = "ordinary"
    @testset "$l" for l in labels
        lat =  lats[i]
        rlat = make_recip_latvecs(lat,convention)
        bz = make_bz(lat,convention)
        ch = chull(bz)
        @test ch.volume ≈ det(rlat)

        spaceGroup = symmetry.get_spaceGroup(lats[i]',at,atom_pos)[1]
        ops = [spaceGroup[j,:,:] for j=1:size(spaceGroup)[1]]
        unfoldedpts=unique([op*bz[j,:] for op in ops, j in 1:size(bz,1)])
        @test all([any([unfoldedpts[j] ≈ bz[k,:] for j in
            1:length(unfoldedpts)]) for k in 1:size(bz,1)])
        i+=1
    end
end
