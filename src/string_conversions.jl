
const to_subscr = Dict(
    '0' => '₀',
    '1' => '₁',
    '2' => '₂',
    '3' => '₃',
    '4' => '₄',
    '5' => '₅',
    '6' => '₆',
    '7' => '₇',
    '8' => '₈',
    '9' => '₉',
)

const from_subscr = [v=>k for (k, v) in to_subscr]

function n2s(x)
    x = round(x; digits=3)
    if isinteger(x)
        x = Int(x)
    end
    x == 1 && return ""
    return string(x)
end

de_subscr(s) = replace(s, from_subscr...)

