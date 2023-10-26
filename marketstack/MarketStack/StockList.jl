#  Just the List of Stocks

using CSV
using DataFrames

bse = "BSE.Securities"

d = CSV.read(bse,DataFrame)


