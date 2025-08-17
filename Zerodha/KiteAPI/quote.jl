"""
    last_trading_price(exchange::String,stock::String)

This API returns the complete market data snapshot of up to 500 instruments in one go. It includes the quantity, OHLC, and Open Interest fields, and the complete bid/ask market depth amongst others.

Instruments are identified by the exchange:tradingsymbol combination and are passed as values to the query parameter i which is repeated for every instrument. If there is no data available for a given key, the key will be absent from the response. 
Get LTP of trading symbol

`instrument` should be in EXCHANGE:TRADINGSYMBOL format, eg. "NSE:INFY"
"""
function last_trading_price(conn::KiteConnection,exchange::String,stock::String)

  instrument = exchange*":"*stock
  url_fragment = "/quote/ltp?i=$instrument"
  url = API_ENDPOINT*url_fragment 
  # header  = [ "X-Kite-Version" => "3",
  #             "Content-Type"   => "application/x-www-form-urlencoded",
  #             "Authorization"  => "token $(API_KEY):$(ACCESS_TOKEN)"
  #           ]
  header    = kite_std_header(conn)

  res = HTTP.get(url,header)

  last_price = res["data"][instrument]["last_price"]

  return last_price
end
#---------------------------------------------------------------------- 
