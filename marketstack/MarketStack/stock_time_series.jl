function MS_stock_quote(symbol::String, endpoint::String; client = MSGlobal[], outputsize::String="compact", datatype::Union{String, Nothing}=nothing, parser = "default", exchange::String="")
#    @argcheck in(interval, ["1min", "5min", "15min", "30min", "60min"])
#    @argcheck in(endpoint, ["eod", "intraday", "splits", "dividends", "tickers", "exchanges", "currencies", "timezones"])   
    @argcheck in(endpoint,MSEndPoints)
    if !isempty(exchange)
      @argcheck in(exchange,ExchangeCodes)
    end  
    @argcheck in(outputsize, ["compact", "full"])
    @argcheck in(datatype, ["json", "csv", nothing])

    params = Dict(
        "access_key"=>key(client),
        "symbols"=>symbol
    )

    if !isempty(exchange)
      ep = "exchange/$exchange/$endpoint"
    else
      ep = endpoint
    end

    uri = _build_uri(client.scheme, client.host, ep, params)
    println(uri)
    data = retry(_get_request, delays=Base.ExponentialBackOff(n=3, first_delay=5, max_delay=30))(uri)
    p = _parser(parser, datatype)
    return p(data)
end


function MS_time_series_intraday_extended(symbol::String, interval::String="60min", slice::String="year1month1"; client = MSGLOBAL[], parser = "default")
    @argcheck in(interval, ["1min", "5min", "15min", "30min", "60min"])
    sliceMatch = match(r"year(?<year>\d+)month(?<month>\d+)", slice)
    @argcheck !Compat.isnothing(sliceMatch)
    @argcheck parse(Int, sliceMatch["year"]) > 0
    @argcheck parse(Int, sliceMatch["year"]) < 3
    @argcheck parse(Int, sliceMatch["month"]) > 0
    @argcheck parse(Int, sliceMatch["month"]) < 13
    params = Dict(
        "function"=>"TIME_SERIES_INTRADAY_EXTENDED",
        "symbol"=>symbol,
        "interval"=>interval,
        "slice"=>slice,
        "apikey"=>key(client)
    )
    uri = _build_uri(client.scheme, client.host, "query", params)
    return uri
#    data = retry(_get_request, delays=Base.ExponentialBackOff(n=3, first_delay=5, max_delay=1000))(uri)
#    p = _parser(parser, "csv")
#    return p(data)
end

function MS_time_series_intraday(symbol::String, interval::String="1min"; client = MSGLOBAL[], outputsize::String="compact", datatype::Union{String, Nothing}=nothing, parser = "default")
    @argcheck in(interval, ["1min", "5min", "15min", "30min", "60min"])
    @argcheck in(outputsize, ["compact", "full"])
    @argcheck in(datatype, ["json", "csv", nothing])
    params = Dict(
        "function"=>"TIME_SERIES_INTRADAY",
        "symbol"=>symbol,
        "interval"=>interval,
        "outputsize"=>outputsize,
        "datatype"=>datatype,
        "apikey"=>key(client)
    )
    uri = _build_uri(client.scheme, client.host, "query", params)
    return uri
#    data = retry(_get_request, delays=Base.ExponentialBackOff(n=3, first_delay=5, max_delay=1000))(uri)
#    p = _parser(parser, datatype)
#    return p(data)
end

for func in (:daily, :daily_adjusted, :weekly, :weekly_adjusted, :monthly, :monthly_adjusted)
    x = "MS_time_series_$(func)"
    fname = Symbol(x)
    @eval begin
        function ($fname)(symbol::String; client = MSGLOBAL[], outputsize::String="compact", datatype::Union{String, Nothing}=nothing, parser = "default")
            @argcheck in(outputsize, ["compact", "full"])
            @argcheck in(datatype, ["json", "csv", nothing])
            params = Dict(
                "function"=>uppercase($x),
                "symbol"=>symbol,
                "outputsize"=>outputsize,
                "datatype"=>isnothing(datatype) ? "csv" : datatype,
                "apikey"=>key(client)
            )
            uri = _build_uri(client.scheme, client.host, "query", params)
            data = retry(_get_request, delays=Base.ExponentialBackOff(n=3, first_delay=5, max_delay=1000))(uri)
            p = _parser(parser, datatype)
            return p(data)
        end

        export $fname
    end
end





