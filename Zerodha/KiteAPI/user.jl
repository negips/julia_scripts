"""
    get_user_profile(conn::KiteConnection)

While a successful token exchange returns the full user profile, it's possible to retrieve it any point of time with the /user/profile API. Do note that the profile API does not return any of the tokens.

"""
function get_user_profile(conn::KiteConnection)

  url_fragment    = "/user/profile"
  url             = API_ENDPOINT*url_fragment 
  header          = kite_std_header(conn)
  resp            = HTTP.get(url,header)
  jresp           = JSON.parse(String(resp.body))

  return jresp
end
#---------------------------------------------------------------------- 
"""
    get_user_margins(conn::KiteConnection,segment::String)

A GET request to /user/margins/:segment returns funds, cash, and margin information for the user. Segment in the URI can be either equity or commodity.

"""
function get_user_margins(conn::KiteConnection,segment::String)

  @assert (segment == "equity" || segment == "commodity") "segment must be either commodity or equity"

  url_fragment    = "/user/margins/"*segment
  url             = API_ENDPOINT*url_fragment 
  header          = kite_std_header(conn)
  resp            = HTTP.get(url,header)
  jresp           = JSON.parse(String(resp.body))

  return jresp
end
#---------------------------------------------------------------------- 





