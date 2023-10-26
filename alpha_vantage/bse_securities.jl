println("AlphaVantage Testing")


using AlphaVantage
using DataFrames, PyPlot, Dates
using JSON,HTTP
using CSV

close("all")

fname = "BSE.Securities"
fname = "MarketWatch.csv"
securities = CSV.read(fname,DataFrame)

symbols = securities[!,"Security Name"]
#names   = securities[!,"Security Name"]
nsec    = length(symbols)
volumes = zeros(Float64,nsec)

options = Dict(
                "exchange"    => "XNAS",
                "limit"       => "100",
                "offset"      => "0",
               ) 


symbols2 = symbols[1:10]

NF = [""] #["DATAINFRA","JIOFIN"]

for i in 1:length(symbols2)
  symbol = symbols2[i]
  suff   = ".BSE"
  syms   = "$symbol$suff"

  try
    stock_quote         = time_series_daily(syms);
    data                = DataFrame(stock_quote)
    if columnindex(data, :volume)>0      
      vol                = data.volume[1]  
      volumes[i]         = vol
      println("$i $symbol: $vol")
    else
      println("$i $symbol: Volume data not Found")
      println(typeof(symbol))
      append!(NF,symbol)
    end  
  catch
    println("$symbol: Stock Data not found")
    append!(NF,symbol)
  end  
      
  sleep(15)

  
end  



#  dts = today() #Date("2023-09-06")
#  dte = Date("2023-09-05")















