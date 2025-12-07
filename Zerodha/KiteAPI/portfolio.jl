"""
    get_holdings(conn::KiteConnection)

Holdings contain the user's portfolio of long term equity delivery stocks. An instrument in a holdings portfolio remains there indefinitely until its sold or is delisted or changed by the exchanges. Underneath it all, instruments in the holdings reside in the user's DEMAT account, as settled by exchanges and clearing institutions.

"""
function get_holdings(conn::KiteConnection)

  url_fragment    = "/portfolio/holdings"
  url             = API_ENDPOINT*url_fragment 
  header          = kite_std_header(conn)
  resp            = HTTP.get(url,header)
  jresp           = JSON.parse(String(resp.body))

  return jresp
end
#---------------------------------------------------------------------- 
"""
    get_auctions(conn::KiteConnection)

This API returns a list of auctions that are currently being held, along with details about each auction such as the auction number, the security being auctioned, the last price of the security, and the quantity of the security being offered. Only the stocks that you hold in your demat account will be shown in the auctions list.

"""
function get_auctions(conn::KiteConnection)

  url_fragment    = "/portfolio/holdings/auctions"
  url             = API_ENDPOINT*url_fragment 
  header          = kite_std_header(conn)
  resp            = HTTP.get(url,header)
  jresp           = JSON.parse(String(resp.body))

  return jresp
end
#----------------------------------------------------------------------
"""
    get_positions(conn::KiteConnection)

Positions contain the user's portfolio of short to medium term derivatives (futures and options contracts) and intraday equity stocks. Instruments in the positions portfolio remain there until they're sold, or until expiry, which, for derivatives, is typically three months. Equity positions carried overnight move to the holdings portfolio the next day.

The positions API returns two sets of positions, net and day. net is the actual, current net position portfolio, while day is a snapshot of the buying and selling activity for that particular day. This is useful for computing intraday profits and losses for trading strategies.

"""
function get_positions(conn::KiteConnection)

  url_fragment    = "/portfolio/positions"
  url             = API_ENDPOINT*url_fragment 
  header          = kite_std_header(conn)
  resp            = HTTP.get(url,header)
  jresp           = JSON.parse(String(resp.body))

  return jresp
end
#----------------------------------------------------------------------



