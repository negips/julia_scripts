println("MarketStack Testing")

include("MarketStack.jl")

using .MarketStack
using DataFrames, PyPlot, Dates
using JSON,HTTP
using CSV

close("all")

fname = "BSE.Securities"
securities = CSV.read(fname,DataFrame)

symbols = securities[!,"Security Id"]
names   = securities[!,"Security Name"]
nsec    = length(symbols)
volumes = zeros(Float64,nsec)

options = Dict(
                "exchange"    => "XNAS",
                "limit"       => "100",
                "offset"      => "0",
               ) 


suff = ".XBOM"
symbols2 = symbols[1:100]
s    = symbols[1]
syms = "$s$suff"

for i in 2:length(symbols2)
  global syms
  symbol = symbols2[i]
#  suff = ".XBOM"
  syms = "$syms,$symbol$suff"
  endpoint = "eod"

#  stock_quote = MS_stock_quote(sym,endpoint,datatype="json",Options=options)
#  global data  = DataFrame(stock_quote["data"])

#  volumes[i] = data.volume[1]
end  

syms = "AAPL"
endpoint = "eod"

stock_quote = MS_stock_quote(syms,endpoint,datatype="json",Options=options)
data  = DataFrame(stock_quote["data"])

dates = data.date

for i in 1:length(dates)
  d = dates[i][1:10]
  dates[i] = d
end  

dates = Date.(dates)


#  dts = today() #Date("2023-09-06")
#  dte = Date("2023-09-05")















