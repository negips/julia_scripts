"""
  'kite_get_tokens(api_key::String, api_secret::String)`

Setup your KiteAPI session by providing your API key and
API secret which you get from Zerodha
"""
function kite_get_tokens(api_key::String,api_secret::String)

  # Request token
  rtoken = gen_request_token(api_key)
  # Access token
  atoken = gen_access_token(api_key,api_secret,rtoken)

  return rtoken,atoken
end

"""
  'gen_request_token()`

Generate the request token
"""
function gen_request_token(key::String)

   url = "https://kite.trade/connect/login?api_key="*key
   # println(url)
   @printf("Login at: %s\n",url)
   
   @printf("Enter request_token:\n")
   request_token = readline()
   @printf("Request Token: %s\n",request_token)

   return request_token
end  
#---------------------------------------------------------------------- 
"""
  `gen_access_token(request_token::String)`

Generate the access token by passing in yout request token
"""
function gen_access_token(key::String,secret::String,request_token::String)
  
  checksum  = bytes2hex(sha256(key * request_token * secret))
  url       = "$API_ENDPOINT/session/token"
  header    = [ "X-Kite-Version" => X_KITE_VER,
                "Content-Type"   => "application/x-www-form-urlencoded" ]
  body      = "api_key=$(key)&request_token=$(request_token)&checksum=$(checksum)"

  res       = HTTP.post(url, header, body)
  jres      = JSON.parse(String(res.body))
  access_token = jres["data"]["access_token"]

  return access_token
end
#---------------------------------------------------------------------- 
"""
  `kite_std_header(conn::KiteConnection)`

Create the standard header for HTTP GET/PUT/POST
"""
function kite_std_header(conn)

  key       = conn.api_key
  atoken    = conn.access_token
  header  = [ "X-Kite-Version" => X_KITE_VER,
              "Content-Type"   => "application/x-www-form-urlencoded",
              "Authorization"  => "token $(key):$(atoken)"
            ]

  return header
end
#---------------------------------------------------------------------- 
"""
  `kite_close_session(conn::KiteConnection)`

This call invalidates the access_token and destroys the API session. After this, the user should be sent through a new login flow before further interactions. This does not log the user out of the official Kite web or mobile applications.
"""
function kite_close_session(conn)

  key       = conn.api_key
  atoken    = conn.access_token
  header    = [ "X-Kite-Version" => "3" ]
  url       = "$API_ENDPOINT/session/token?api_key=$(key)&access_token=$(atoken)"
  resp      = HTTP.request("DELETE", url, header)
  jresp     = JSON.parse(String(resp.body))
  
  return jresp
end
#---------------------------------------------------------------------- 







