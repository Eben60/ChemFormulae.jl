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

Base.hash(e::ChemFormula, h::UInt) = hash(e.brutto_string, h)
Base.isequal(e1::ChemFormula, e2::ChemFormula) = e1.atoms == e2.atoms
Base.:(==)(e1::ChemFormula, e2::ChemFormula) = isequal(e1, e2)

# to be used for sorting only
Base.isless(e1::ChemFormula, e2::ChemFormula) = e1.brutto_string < e2.brutto_string
