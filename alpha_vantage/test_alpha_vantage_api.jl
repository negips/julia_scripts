println("Alpha Vantage Testing")

using AlphaVantage
using DataFrames, PyPlot, Dates

#gr(size=(800,470))
# Get daily S&P 500 data
spy = time_series_daily("SPY");
# Convert to a DataFrame
data = DataFrame(spy);
# Convert timestamp column to Date type
data[!, :timestamp] = Dates.Date.(data[!, :timestamp]);
data[!, :open] = Float64.(data[!, :open])
# Plot the timeseries
plot(data[!, :timestamp], data[!, :open], label=["Open"])
savefig("sp500.png")

# using JSON
# using HTTP
# 
# API_KEY = get(ENV,"ALPHA_VANTAGE_API_KEY","")
# 
# println(API_KEY)
# 
# url = "https://www.alphavantage.co"
# params = Dict("function" => "GLOBAL_QUOTE",
#               "symbol" => "IBM",
#               "datatype" => "JSON",
#               "apikey" => API_KEY)
# 
# js = JSON.json(params);
# 
# res = HTTP.request("POST",url,"Content-Type" => "application/json",js)













