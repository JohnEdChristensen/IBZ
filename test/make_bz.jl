using Test,IBZ,PyCall
include("testHelperFunctions.jl")
lats = lattices
labels =lattice_labels

@testset "make_bz_test" begin
    i=1
    @testset "$l" for l in labels
        lat =  lats[i]
        at = [1]
        atom_pos = [0 0 0]
        make_IBZ = pyimport("bzip.make_IBZ")
        ch = make_IBZ.find_bz(lat)
        #convert to julia
        ch = chull(ch.points)
        chjulia = chull(make_bz(lat))
        @test ch.volume ≈  chjulia.volume
        @test comparePolygon(ch, chjulia)
        i+=1
    end
end

