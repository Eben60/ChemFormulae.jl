using ChemFormulae, ChemElementsBB
using ChemFormulae: ElemInCompound

using Test

@testset "ChemFormulae.jl" begin


abellait = ChemFormula("NaPb₂(CO3)₂(OH)");
@test abellait.brutto_string == "HC2O7NaPb2"
@test abellait.weight ≈ 574.41276928
@test abellait.atoms_total == 13

x = abellait[1];
@test x.elem.name == "Hydrogen"
@test x.mass_share ≈ 0.0017548356406900244

pb = abellait[:Pb];
@test pb.weight == 414.4
@test pb.mass_share ≈ 0.7214324300614544

@test haskey(abellait, 3)
@test haskey(abellait, 5)
@test ! haskey(abellait, 6)
@test haskey(abellait, :C)
@test ! haskey(abellait, :Li)

@test length(abellait) == 5
@test eltype(abellait) == typeof(pb)
@test firstindex(abellait) == 1
@test lastindex(abellait) == 5

@test abellait[1:2] == abellait[[1,2]] == abellait[[:H, :C]]

s = 0
for e in abellait
    s += e.n
end
@test s == 13

end
;