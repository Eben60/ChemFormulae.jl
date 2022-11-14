
mutable struct ElemInCompound
    const elem::ChemElemBB
    const n::Float64
    const weight::Float64
    mass_share::Float64
end

struct ChemFormula
    cc_string::String
    atoms::Vector{ElemInCompound}
    weight::Float64
end

ElemInCompound(elem::ChemElemBB, n) = ElemInCompound(elem::ChemElemBB, n, n*elem, NaN)

function ChemFormula(f::AbstractString)
    c = Compound(f)
    atoms = ElemInCompound[]
    for a in c.tuples
        sym = Symbol(a[1])
        n = a[2]
        eic = ElemInCompound(chem_els[sym], n)
        push!(atoms, eic)
    end
    weight = sum(x -> x.elem, atoms)

    for a in atoms
        a.mass_share = a.weight / weight
    end
    return ChemFormula(f, atoms, weight)
end
