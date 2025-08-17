
"""
    last_trading_price(exchange::String,stock::String)

Get LTP of trading symbol

`instrument` should be in EXCHANGE:TRADINGSYMBOL format, eg. "NSE:INFY"
"""
function last_trading_price(exchange::String,stock::String)
  # parts = split(instrument, ":")
  # if (length(parts) != 2)
  #   throw(ArgumentError("instrument should be in EXCHANGE:TRADINGSYMBOL format"))
  # end

  instrument = exchange*":"*stock
  url_fragment = "/quote/ltp?i=$instrument"
  url = API_ENDPOINT*url_fragment 
  header  = [ "X-Kite-Version" => "3",
              "Content-Type"   => "application/x-www-form-urlencoded",
              "Authorization"  => "token $(API_KEY):$(ACCESS_TOKEN)"
            ]

  res = HTTP.get(url,header)

  last_price = res["data"][instrument]["last_price"]

  return last_price
end

