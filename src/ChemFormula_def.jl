
struct ElemInCompound
    elem::ChemElemBB
    n::Float64
    weight::Float64
    mass_share::Float64
end

ElemInCompound(elem::ChemElemBB, n) = ElemInCompound(elem::ChemElemBB, n, n*elem, NaN)

struct ChemFormula
    cc_string::String
    brutto_string::String
    weight::Float64
    atoms::Vector{ElemInCompound}
    bysymbol::Dict{Symbol, Int}
end

el_in_comp_substr(e::ElemInCompound) = "$(e.elem.symbol)$(n2s(e.n))"

function ChemFormula(f::AbstractString)
    f0 = de_subscr(f)
    c = Compound(f0)

    ats = [(;elem=chem_els[Symbol(a[1])], n = a[2]) for a in c.tuples]
    weight = sum(x -> x.elem*x.n, ats)
    atoms = [ElemInCompound(x.elem, x.n, x.elem.weight*x.n, x.elem.weight*x.n/weight) for x in ats]

    sort!(atoms; by = x -> x.elem.number)
    els = [el_in_comp_substr(x) for x in atoms]
    brutto_string = join(els, "")

    bysymbol = Dict{Symbol, Int}(atoms[i].elem.symbol=>i for i in eachindex(atoms))
    return ChemFormula(f, brutto_string, weight, atoms, bysymbol)   
end

Base.getindex(e::ChemFormula, i::Integer) = e.atoms[i]
Base.getindex(e::ChemFormula, i::Symbol) = e.atoms[e.bysymbol[i]]


# TODO indexing, iterations