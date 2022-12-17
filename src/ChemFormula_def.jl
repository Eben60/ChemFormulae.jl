
struct ElemInCompound
    elem::ChemElem
    n::Float64
    weight::Float64
    mass_share::Float64
end

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
    return Vector{Pair{Symbol, Union{Float64, Int}}}([Symbol(a[1]) => a[2] for a in c])
end

function ChemFormula(f::AbstractString; groups = nothing)
    if  ! isnothing(groups) 
        groups = collect(pairs(groups))
        groups = [string(p.first) => "($(p.second))" for p in groups ]
        f1 = add1(f)
        f = replace(f1, groups...)
    end
    ps = parse_chemformula(f)
    return ChemFormula(ps, f)
end

weight_uless(elem) = elem.atomic_weight |> ustrip

function string_or_missing(s)
    ismissing(s) && return s
    return String(s)
end

function ChemFormula(ps::Vector{Pair{S, I}}, f = missing) where I <: Real where S <: Union{Symbol, Integer}
    ats = [(;elem=chem_elements[a[1]], n = a[2]) for a in ps]
    weight = sum(x -> weight_uless(x.elem)*x.n, ats)
    atoms = [ElemInCompound(x.elem, x.n, weight_uless(x.elem)*x.n, weight_uless(x.elem)*x.n/weight) 
                            for x in ats]
    sort!(atoms; by = x -> x.elem.atomic_number)
    atoms_total = sum(x -> x.n, atoms)
    brutto_string = join([el_in_comp_substr(x) for x in atoms], "")
    bysymbol = Dict{Symbol, Int}(atoms[i].elem.symbol=>i for i in eachindex(atoms))
    return ChemFormula(string_or_missing(f), brutto_string, atoms_total, weight, atoms, bysymbol)   
end

ChemFormula(p::Union{Dict, NamedTuple}, f = missing) = ChemFormula(collect(pairs(p)), f)

macro cf_str(str) ChemFormula(str) end