
struct OurRational{T<:Integer} <: Real
    num::T
    den::T
    function OurRational{T}(num::T, den::T) where T<:Integer
        if num == 0 && den == 0
             error("invalid rational: 0//0")
        end
        num = flipsign(num, den)
        den = flipsign(den, den)
        g = gcd(num, den)
        num = div(num, g)
        den = div(den, g)
        new(num, den)
    end
end

OurRational(n::T, d::T) where {T<:Integer} = OurRational{T}(n,d)
