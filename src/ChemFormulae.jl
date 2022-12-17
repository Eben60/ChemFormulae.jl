module ChemFormulae

# using ChemElementsBB # , ChemEquations
using Mendeleev, Unitful

include("compound.jl")
include("string_conversions.jl")
include("ChemFormula_def.jl")
include("overloads.jl")

export ChemFormula, @cf_str

end
