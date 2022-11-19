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

h2o_1 = ChemFormula("H2O")
h2o_2 = ChemFormula([:H=>2, :O=>1], "H2O")
h2o_3 = ChemFormula([1=>2, 8=>1], "H2O")
h2o_4 = ChemFormula(Dict(:H=>2, :O=>1), "H2O")
h2o_5 = ChemFormula((;O=1, H=2))

@test h2o_1.cc_string == h2o_2.cc_string

@test h2o_1.brutto_string == h2o_2.brutto_string
@test h2o_1.atoms_total == h2o_2.atoms_total
@test h2o_1.weight == h2o_2.weight
@test h2o_1.atoms == h2o_2.atoms
@test h2o_1.bysymbol == h2o_2.bysymbol

@test isequal(h2o_1, h2o_2)
@test isequal(h2o_1, h2o_3)
@test isequal(h2o_1, h2o_4)
@test isequal(h2o_1, h2o_5)

@test_broken h2o_1 == h2o_2

ybco = ChemFormula([:Y=>1, :Ba=>2, :Cu=>3, :O=>6.94])
@test ybco.brutto_string == "O6.94Cu3YBa2"
@test ybco.atoms_total isa Float64
@test ybco.atoms_total ≈ 12.94

end
;