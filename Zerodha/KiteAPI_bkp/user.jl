"""
    get_user_profile()

Get the user profile

"""
function get_user_profile()

  url_fragment = "/user/profile"
  url = API_ENDPOINT*url_fragment 
  header  = [ "X-Kite-Version" => "3",
              "Content-Type"   => "application/x-www-form-urlencoded",
              "Authorization"  => "token $(API_KEY):$(ACCESS_TOKEN)"
            ]

  resp      = HTTP.get(url,header)
  jresp     = JSON.parse(String(resp.body))

  # last_price = res["data"][instrument]["last_price"]

  return jresp
end
#---------------------------------------------------------------------- 
"""
    get_user_margins(segment::String)

Get the user funds and margins

"""
function get_user_margins(segment::String)

  @assert (segment == "equity" || segment == "commodity") "segment must be either commodity or equity"   

  url_fragment = "/user/margins/"*segment
  url = API_ENDPOINT*url_fragment 
  header  = [ "X-Kite-Version" => "3",
              "Content-Type"   => "application/x-www-form-urlencoded",
              "Authorization"  => "token $(API_KEY):$(ACCESS_TOKEN)"
            ]

  resp      = HTTP.get(url,header)
  jresp     = JSON.parse(String(resp.body))

  return jresp
end

