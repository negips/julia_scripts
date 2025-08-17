"""
  'init(api_key::String, api_secret::String)`

Setup your KiteAPI session by providing your API key and
API secret which you get from Zerodha
"""
function init(api_key::String, api_secret::String)
  global API_KEY = api_key
  global API_SECRET = api_secret
end

"""
  'gen_request_token()`

Setup your KiteAPI session by providing your API key and
API secret which you get from Zerodha
"""
function gen_request_token()

   url = "https://kite.trade/connect/login?api_key="*API_KEY
   # println(url)
   @printf("Login at: %s\n",url)
   
   @printf("Enter request_token:\n")
   request_token = readline()
   @printf("Request Token: %s\n",request_token)

   return request_token
end  
#---------------------------------------------------------------------- 


function get_http_headers()
  hdr = [ "X-Kite-Version" => "3",
          "Authorization" => "token $API_KEY:$ACCESS_TOKEN"
        ]
  return hdr
end

function http_get(url_fragment::String)

  url = "$API_ENDPOINT/$url_fragment"
  r = HTTP.request("GET", url, get_http_headers())

  return JSON.parse(String(r.body))
end

"""
  `gen_access_token(request_token::String)`

Generate the access token by passing in yout request token
"""
function gen_access_token(request_token::String)
  checksum  = bytes2hex(sha256(API_KEY * request_token * API_SECRET))
  url       = "$API_ENDPOINT/session/token"
  header    = [ "X-Kite-Version" => "3",
                "Content-Type"   => "application/x-www-form-urlencoded" ]
  body      = "api_key=$API_KEY&request_token=$request_token&checksum=$checksum"

  res       = HTTP.post(url, header, body)
  r         = JSON.parse(String(res.body))
  global ACCESS_TOKEN = r["data"]["access_token"]

  return nothing
end

"""
  `set_access_token(access_token::String)`

Set the access token as a global in the Module
"""
function set_access_token(access_token::String)

  global ACCESS_TOKEN = access_token
end



