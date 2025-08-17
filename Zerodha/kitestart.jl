# Start Kite Connect

using HTTP
using JSON
using SHA
using DotEnv
using Printf

include("KiteAPI/KiteAPI.jl")

DotEnv.load!()

key         = ENV["FIRST_API_KEY"]
secret      = ENV["FIRST_API_SECRET"]

# Creates tokens
kite        = KiteAPI.KiteConnection(key,secret)

@printf("Access Token: %s\n",kite.access_token)

sto   = "BEL"
Ex    = "NSE"

#ltp   = KiteAPI.last_trading_price(Ex,sto)
profile   = KiteAPI.get_user_profile(kite)
equity    = KiteAPI.get_user_margins(kite,"equity")
commodity = KiteAPI.get_user_margins(kite,"commodity")
holdings  = KiteAPI.get_holdings(kite)
positions = KiteAPI.get_positions(kite)

# @printf("LTP of Stock %s at Exchange %s = %s\n",sto,Ex,ltp)

# Close session.
res       = KiteAPI.kite_close_session(kite)

println("Done.")





