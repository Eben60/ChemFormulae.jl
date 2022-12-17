# adapted from https://github.com/zlatanvasovic/ChemEquations.jl by Zlatan Vasović 


"Type stored in `Compound.tuples`."
const ElementTuple = Tuple{String, Union{Int, Float64}}

"Regex to match `{...}` charge string."
const CHARGEREGEX = r"{(.*)}"

"""
Stores chemical compound's elements and charge in a structured way.

!!! info
    Electron is stored as `"e"`.
"""
struct Compound
    tuples::Vector{ElementTuple}
    charge::Int
end

"""
Constructs a compound from `str`.

An element begins with an uppercase unicode letter and
ends with a lowercase unicode letter or a unicode symbol.

!!! info
    An element can also begin with a symbol if
    the symbol is the first character (e.g. `"⬡H"`).

Parsing is insensitive to whitespace and underscores (`_`),
but also to state symbols (`(s)`, `(l)`, `(g)`, `(aq)`).
Special parsing is implemented for:
- parens (e.g. `"(CH3COO)2Mg"`)
- compounds with `"*"` (e.g. `"CuSO4 * 5H2O"`)
- electrons (`"e"`)
Charge is in the form `"{±n}"` or `"{n±}"`.
It is automatically deduced for electron (`"e"`).

# Examples
```jldoctest
julia> Compound("H2O(s)")
H2O

julia> Compound("H3O{+}")
H3O{+}

julia> Compound("(CH3COO)2Mg")
C4H6O4Mg

julia> Compound("CuSO4 * 5H2O")
CuSO9H10

julia> Compound("⬡Cl")
⬡Cl
```
"""
function Compound(str::AbstractString)
    str = replace(str, [' ', '_'] => "")
    str = replace(str, r"\((s|l|g|aq)\)$" => "") # remove state symbols
    Compound(_elementtuples(str), _charge(str))
end

"Add 1 to elements and parens without a coefficient"
add1(str) = replace(str,  r"(?<x>\p{L}|\p{S}|\))(?=(\p{Lu}|\(|\)|\*|$))" => s"\g<x>1")

"Extracts element tuples from compound's string."
function _elementtuples(str::AbstractString)
    str = replace(str, CHARGEREGEX => "")
    if str ∈ ("", "e")
        return [("e", 1)]
    end

    tuples = ElementTuple[]

    # Add 1 to elements and parens without a coefficient
    str = add1(str)
#   str = replace(str,  r"(?<x>\p{L}|\p{S}|\))(?=(\p{Lu}|\(|\)|\*|$))" => s"\g<x>1")

    # Expand parens
    for substr ∈ eachmatch(r"\(((\p{L}|\d)+)\)(\d+)", str)
        capture = replace(
            substr.captures[1],
            isdigit => x -> string(parse(Int, substr.captures[3]) * parse(Int, x))
        )
        str = replace(str, substr.match => capture)
    end

    # Expand compounds with '*', e.g. CuSO4*5H2O

    if occursin('*', str)
        str = split(str, '*')
        for (i, substr) ∈ enumerate(str)
            if isdigit(substr[1])
                k, substr = match(r"(^\d+)(.+)", substr).captures
                substr = replace(substr, isdigit => x -> string(parse(Int, k) * parse(Int, x)))
                str[i] = substr
            end
        end
        str = join(str)
    end

    str = split(str, r"(?=\p{Lu})")
    for element ∈ str
        element, k = split(element, r"(?=[\d\.])(?<![\d\.])")

        index = findfirst(x -> x[1] == element, tuples)
        n = parse(Float64, k)
        isinteger(n) && (n = Int(n))

        if isnothing(index)
            push!(tuples, (element, n))
        else
            tuples[index] = (element, tuples[index][2] + n)
        end
    end

    return tuples
end

"Extracts charge from compound's string into a number of specified type."
function _charge(str::AbstractString)
    if str == "e"
        return -1
    end

    strmatch = match(CHARGEREGEX, str)
    if isnothing(strmatch)
        return 0
    else
        str = strmatch.captures[1]
        if str ∈ ("-", "+")
            str *= "1"
        elseif str[end] ∈ ('-', '+')
            str = str[end] * str[1:end-1]
        end
        return parse(Int, str)
    end
end

"""
Constructs a compound with `cc"str"` syntax, instead of `Compound(str)`.

# Examples
```jldoctest
julia> cc"H3O{+1}"
H3O{+}
```
"""
macro cc_str(str) Compound(str) end

"""
Checks whether two compounds are chemically equal.

# Examples
```jldoctest
julia> cc"MgOHOH" == cc"Mg(OH)2"
true
```
"""
function Base.:(==)(compound_1::Compound, compound_2::Compound)
    return sort(compound_1.tuples) == sort(compound_2.tuples) &&
        compound_1.charge == compound_2.charge
end

"""
Creates a string to represent the compound.

All elements are displayed only once (e.g. `"H2O"` and not `"HOH"`),
in the order in which they were originally given (e.g. `"MgO2H2"` from `cc"Mg(OH)2"`),
with coefficients equal to 1 not displayed (e.g. `"H"` and not `"H1"`).

# Examples
```jldoctest
julia> string(cc"CuSO4 * 5 H2O")
"CuSO9H10"
```
"""
function Base.string(compound::Compound)
    if compound == cc"e"
        return "e"
    end

    str = ""
    for (element, k) ∈ compound.tuples
        str *= element
        if k ≠ 1
            str *= string(k)
        end
    end
    if hascharge(compound)
        str *= "{"
        if compound.charge > 0
            str *= "+"
        end
        if compound.charge == -1
            str *= "-"
        elseif compound.charge ≠ 1
            str *= string(compound.charge)
        end
        str *= "}"
    end
    return str
end

"Displays the compound using [`Base.string(::Compound)`](@ref)."
function Base.show(io::IO, compound::Compound)
    print(io, string(compound))
end

"""
Returns compound's elements as strings.

# Examples
```jldoctest
julia> elements(cc"CH3COOH")
3-element Array{String,1}:
 "C"
 "H"
 "O"
```
"""
function elements(compound::Compound)
    return first.(compound.tuples)
end

"True if the compound's charge is nonzero."
function hascharge(compound::Compound)
    return compound.charge ≠ 0
end
