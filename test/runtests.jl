using ChemFormulae, ChemElementsBB
using ChemFormulae: ElemInCompound

using Test

@testset "ChemFormulae.jl" begin
    B2 = ElemInCompound(chem_els.B, 2)
    @test B2.n == 2
    @test B2.weight == 21.62
    @test isnan(B2.mass_share) 

end
