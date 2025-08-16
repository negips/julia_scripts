"""
    KiteAPI

Julia module to interface with Zerodha's KiteConnect API
"""
module KiteAPI

using HTTP
using JSON
using SHA
using Printf

const API_ENDPOINT = "https://api.kite.trade"
API_KEY = ""
API_SECRET = ""
ACCESS_TOKEN = ""

include("connect.jl")
include("quote.jl")
include("user.jl")
include("portfolio.jl")

# Connect
export      init, 
            gen_request_token,
            gen_access_token, 
            set_access_token

# Quote
export      last_trading_price

# User
export      get_user_profile,
            get_user_margins

# Portfolio
export      get_holdings,
            get_auctions,
            get_positions

end # module




