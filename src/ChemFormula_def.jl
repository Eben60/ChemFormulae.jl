
mutable struct ElemInCompound
    const elem::ChemElemBB
    const n::Float64
    const weight::Float64
    mass_share::Float64
end

ElemInCompound(elem::ChemElemBB, n) = ElemInCompound(elem::ChemElemBB, n, n*elem, NaN)

struct ChemFormula
    cc_string::String
    brutto_string::String
    atoms::Vector{ElemInCompound}
    weight::Float64
end

el_in_comp_substr(e::ElemInCompound) = "$(e.elem.symbol)$(n2s(e.n))"

function ChemFormula(f::AbstractString)
    f = de_subscr(f)
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
    sort!(atoms; by = x -> x.elem.number)
    els = [el_in_comp_substr(x) for x in atoms]
    brutto_string = join(els, "")
    return ChemFormula(f, brutto_string, atoms, weight)   
end
