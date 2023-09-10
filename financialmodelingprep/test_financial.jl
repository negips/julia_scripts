#!/usr/bin/julia

using FinancialModelingPrep
using DataFrames
using PyPlot
using Dates

# load your API key
# 80266911b9310375ee54298ce5dd5387
#FMP_API_KEY = ENV["FMP_API_KEY"]

# create a new FMP API instance
#fmp = FMP(apikey = "80266911b9310375ee54298ce5dd5387" )

close("all")

Stock_Name = "MSFT"

income_statement = income_statements(fmp, Stock_Name)
cash_flow_statement = cash_flow_statements(fmp, Stock_Name)
price = price_quote(fmp, Stock_Name)
hist_price = historical_price_quote(fmp, Stock_Name)

df1 = DataFrame(income_statement)
df2 = DataFrame(cash_flow_statement)
df3 = DataFrame(price)
df4 = DataFrame(hist_price)

dt  = DateTime.(df4.date);
stock = df4.close;

plot(dt,stock)







