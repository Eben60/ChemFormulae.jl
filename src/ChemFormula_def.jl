
struct ElemInCompound
    elem::ChemElemBB
    n::Float64
    weight::Float64
    mass_share::Float64
end

# ElemInCompound(elem::ChemElemBB, n) = ElemInCompound(elem::ChemElemBB, n, n*elem, NaN)

struct ChemFormula{T <: Union{Float64, Int}}
    cc_string::Union{String, Missing}
    brutto_string::String
    atoms_total::T
    weight::Float64
    atoms::Vector{ElemInCompound}
    bysymbol::Dict{Symbol, Int}
end

el_in_comp_substr(e::ElemInCompound) = "$(e.elem.symbol)$(n2s(e.n))"

function parse_chemformula(f::AbstractString)
    f0 = de_subscr(f)
    c = Compound(f0).tuples
    return [Symbol(a[1]) => a[2] for a in c]
end


function ChemFormula(f::AbstractString)
    ps = parse_chemformula(f)
    return ChemFormula(ps, f)
end


function ChemFormula(ps::Vector{Pair{S, I}}, f = missing) where I <: Real where S <: Union{Symbol, Integer}
    ats = [(;elem=chem_els[a[1]], n = a[2]) for a in ps]
    weight = sum(x -> x.elem*x.n, ats)
    atoms = [ElemInCompound(x.elem, x.n, x.elem.weight*x.n, x.elem.weight*x.n/weight) 
                            for x in ats]
    sort!(atoms; by = x -> x.elem.number)
    atoms_total = sum(x -> x.n, atoms)
    brutto_string = join([el_in_comp_substr(x) for x in atoms], "")
    bysymbol = Dict{Symbol, Int}(atoms[i].elem.symbol=>i for i in eachindex(atoms))
    return ChemFormula(f, brutto_string, atoms_total, weight, atoms, bysymbol)   
end

ChemFormula(d::Dict, f = missing) = ChemFormula(collect(d), f)
ChemFormula(p::NamedTuple, f = missing) = ChemFormula(collect(pairs(p)), f)