# Start Kite Connect

using HTTP
using JSON
using SHA
using DotEnv
using Printf

include("KiteAPI/KiteAPI.jl")

function kite_gen_request_token(key)

   url = "https://kite.trade/connect/login?api_key="*key
   println(url)

   return url
end  
#---------------------------------------------------------------------- 
function kite_gen_access_token(rtoken,key,secret,endpoint)

  str       = key * rtoken * secret
  # println(str)
  checksum  = bytes2hex(sha256(str))

  url       = "$endpoint/session/token"

  header    = [ "X-Kite-Version" => "3",
                 "Content-Type"   => "application/x-www-form-urlencoded" ]
               
  body      = "api_key=$(key)&request_token=$(rtoken)&checksum=$(checksum)"
  
  resp      = HTTP.post(url,header,body)

  return resp
end  
#----------------------------------------------------------------------

DotEnv.load!()

key         = ENV["KITE_API_KEY"]
secret      = ENV["KITE_API_SECRET"]

KiteAPI.init(key,secret)

#token_url   = kite_gen_request_token(KiteAPI.API_KEY)

# Get Request Token
request_token = KiteAPI.gen_request_token()

# Generate Access Token
KiteAPI.gen_access_token(request_token)

@printf("Access Token: %s\n",KiteAPI.ACCESS_TOKEN)

sto   = "BEL"
Ex    = "NSE"

ltp   = KiteAPI.last_trading_price(Ex,sto)

@printf("LTP of Stock %s at Exchange %s = %s\n",sto,Ex,ltp)

println("Done.")





