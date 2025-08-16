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

export      init, 
            gen_request_token,
            gen_access_token, 
            set_access_token

export      last_trading_price

end # module
