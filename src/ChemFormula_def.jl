
struct ElemInCompound
    elem::ChemElemBB
    n::Float64
    weight::Float64
    mass_share::Float64
end

# ElemInCompound(elem::ChemElemBB, n) = ElemInCompound(elem::ChemElemBB, n, n*elem, NaN)

struct ChemFormula
    cc_string::String
    brutto_string::String
    atoms_total::Int
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

    atoms = [ElemInCompound(x.elem, x.n, x.elem.weight*x.n, x.elem.weight*x.n/weight) 
                            for x in ats]
    sort!(atoms; by = x -> x.elem.number)
    atoms_total = sum(x -> x.n, atoms)

    brutto_string = join([el_in_comp_substr(x) for x in atoms], "")

    bysymbol = Dict{Symbol, Int}(atoms[i].elem.symbol=>i for i in eachindex(atoms))
    return ChemFormula(f, brutto_string, atoms_total, weight, atoms, bysymbol)   
end

Base.getindex(e::ChemFormula, i::Integer) = e.atoms[i]
Base.getindex(e::ChemFormula, i::Symbol) = e.atoms[e.bysymbol[i]]

Base.getindex(e::ChemFormula, v::AbstractVector{I}) where I <: Integer = e.atoms[v]
Base.getindex(e::ChemFormula, v::AbstractVector{Symbol}) = 
        e.atoms[[e.bysymbol[i] for i in v]]

Base.haskey(e::ChemFormula, i::Integer) = firstindex(e.atoms) <= i <= lastindex(e.atoms)
Base.haskey(e::ChemFormula, i::Symbol) = haskey(e.bysymbol, i)

Base.firstindex(e::ChemFormula) = firstindex(e.atoms)
Base.lastindex(e::ChemFormula) = lastindex(e.atoms)
Base.eachindex(e::ChemFormula) = eachindex(e.atoms)

Base.eltype(e::ChemFormula) = ElemInCompound
Base.length(e::ChemFormula) = length(e.atoms)
Base.iterate(e::ChemFormula, state...) = iterate(e.atoms, state...)

Base.isequal(e1::ChemFormula, e2::ChemFormula) = e1.brutto_string == e2.brutto_string
Base.hash(e::ChemFormula, h::UInt) = hash(e.brutto_string, h)

# to be used for sorting only
Base.isless(e1::ChemFormula, e2::ChemFormula) = e1.brutto_string < e2.brutto_string
