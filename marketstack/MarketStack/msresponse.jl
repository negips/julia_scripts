import Tables

struct MarketStackResponse
    data::Vector{AbstractVector}
    names::Vector{AbstractString}

    MarketStackResponse(data::Vector{AbstractVector}, names::Vector{<:AbstractString}) = begin
        l1 = length(data[1])
        @assert all(t -> length(t) == l1, data)
        @assert length(data) == length(names)
        new(data, names)
    end
end

function narrow_types(v::AbstractVector)
    T = mapreduce(typeof, promote_type, v)
    convert(AbstractVector{T}, v)
end

MarketStackResponse(data::AbstractVector{<:AbstractVector{T}} where T, names::AbstractVector{<:AbstractString}) = begin
    MarketStackResponse(narrow_types(data), names)
end

MarketStackResponse(data::AbstractMatrix{T} where T, names::AbstractMatrix{<:AbstractString}) = begin
    v = AbstractVector[narrow_types(c) for c in eachcol(data)]
    n = vec(names)
    MarketStackResponse(v, n)
end

MarketStackResponse(raw::Tuple{AbstractMatrix{Any}, AbstractMatrix{<:AbstractString}}) = begin
    MarketStackResponse(raw[1], raw[2])
end

MarketStackResponse(d::Dict) = Tables.table(reshape(collect(values(d)), (1,:)), header=keys(d))

Tables.istable(::MarketStackResponse) = true
Tables.rowaccess(::MarketStackResponse) = false
Tables.columnaccess(::MarketStackResponse) = true
Tables.columns(t::MarketStackResponse) = t
Tables.getcolumn(t::MarketStackResponse, i::Int) = t.data[i]
Tables.getcolumn(t::MarketStackResponse, nm::Symbol) = begin
    ind = findfirst(==(nm), Symbol.(t.names))
    t.data[ind]
end
Tables.columnnames(t::MarketStackResponse) = Symbol.(t.names)
