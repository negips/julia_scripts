println("Alpha Vantage Testing")

using AlphaVantage
using DataFrames, PyPlot, Dates

#gr(size=(800,470))

# Get daily data

sym = "TITAN.BSE"

spy = time_series_daily(sym);

# Convert to a DataFrame
data = DataFrame(spy);

# Convert timestamp column to Date type
data[!, :timestamp] = Dates.Date.(data[!, :timestamp]);
data[!, :open] = Float64.(data[!, :open])

# Plot the timeseries
plot(data[!, :timestamp], data[!, :open], label=["Open"])

#savefig("sp500.png")














