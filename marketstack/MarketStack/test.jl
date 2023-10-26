println("MarketStack Testing")

include("MarketStack.jl")

using .MarketStack
using DataFrames, PyPlot, Dates
using JSON,HTTP

close("all")

sym = "TITAN.XBOM"
endpoint = "eod"

dts = today() #Date("2023-09-06")
dte = Date("2023-09-05")
options = Dict(
                "exchange"    => "XBOM",
                "limit"       => "1",
                "offset"      => "0",
               ) 

stock_quote = MS_stock_quote(sym,endpoint,datatype="json",Options=options)

data  = DataFrame(stock_quote["data"])

dates = data.date;
open  = data.open
clos  = data.close

plot(dates,open)
plot(dates,clos)


nasdaq_key = get(ENV,"NASDAQ_API_KEY","")

println(nasdaq_key)

host = "https://data.nasdaq.com/api/v3/datasets/"
dataset="bse"
stockcode="BOM500114"
outfmt="data.json"

url = "$host/$dataset/$stockcode/$outfmt?api_key=$nasdaq_key" 

res    = HTTP.get(url)
body   = copy(res.body)  # TODO: re-write to avoid copying
bodyjs = JSON.parse(String(body))

data   = bodyjs["dataset_data"]["data"]
columns = bodyjs["dataset_data"]["column_names"]
dc     = Dict()

l = length(data)
dates       = Vector{String}(undef,l)
openval     = Vector{Float64}(undef,l)
highval     = Vector{Float64}(undef,l)
lowval      = Vector{Float64}(undef,l)
closeval    = Vector{Float64}(undef,l)
wap         = Vector{Float64}(undef,l)
nshares     = Vector{Float64}(undef,l)
ntrades     = Vector{Float64}(undef,l)
turnover    = Vector{Float64}(undef,l)
deliverable = Vector{Float64}(undef,l)
percentage  = Vector{Float64}(undef,l)
spreadhl    = Vector{Float64}(undef,l)
spreadco    = Vector{Float64}(undef,l)


for i in 1:l
  di              = data[i]
  dates[i]        = di[1]
  openval[i]      = di[2]
  highval[i]      = di[2]
  lowval[i]       = di[2]
  closeval[i]     = di[2]
  wap[i]          = di[2]
  nshares[i]      = di[2]
  ntrades[i]      = di[2]
  turnover[i]     = di[2]
  deliverable[i]  = di[2]
  percentage[i]   = di[2]
  spreadhl[i]     = di[2]
  spreadco[i]     = di[2]
end

dateval = Dates.value.(Date.(dates))

close("all")
plot(dateval .- minimum(dateval) ,highval)







