module ChemFormulaTests

using ChemFormulae, ChemElementsBB
using ChemFormulae: ElemInCompound
using Test

@testset "ChemFormula" begin

groups = (;Me = "CH3", Et = "C2H5", Ph = "C6H5")

@test ChemFormula("EtOH"; groups) == ChemFormula("C2H5OH")

end # testset
end # module
;