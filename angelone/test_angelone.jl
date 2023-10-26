println("Testing PyCall for AngelOne")

using PyCall
using JSON
using HTTP

#angel = pyimport("SmartApi")

marketapi  = get(ENV,"ANGELONE_MARKET_API_KEY","")
tradingapi = get(ENV,"ANGELONE_TRADING_API_KEY","")

payload = Dict("mode" => "FULL",
               "exchangeTokens" => Dict("NSE" => "3045") 
		  )

headers = Dict(
    "Content-Type" => "application/json",
    "Accept" => "application/json",
    "X-UserType" => "USER",
    "X-SourceID" => "WEB",
    "X-ClientLocalIP" => "CLIENT_LOCAL_IP",
    "X-ClientPublicIP" => "CLIENT_PUBLIC_IP",
    "X-MACAddress" => "MAC_ADDRESS",
    "X-PrivateKey" => "API_KEY"
   )

host = "https://smartapi.angelbroking.com/publisher-login?api_key=$marketapi"

res = HTTP.request("POST",host,body=payload,headers=headers)

js  = JSON.parse(String(res.body))

#conn = http.client.HTTPSConnection("apiconnect.angelbroking.com")
#payload = {
#    "mode": "FULL",
#    "exchangeTokens": {
#        "NSE": ["3045"]
#    }
#}
#headers = {
#  'X-PrivateKey': 'API_KEY',
#  'Accept': 'application/json',
#  'X-SourceID': 'WEB',
#  'X-ClientLocalIP': 'CLIENT_LOCAL_IP',
#  'X-ClientPublicIP': 'CLIENT_PUBLIC_IP',
#  'X-MACAddress': 'MAC_ADDRESS',
#  'X-UserType': 'USER',
#  'Authorization': 'Bearer AUTHORIZATION_TOKEN',
#  'Accept': 'application/json',
#  'X-SourceID': 'WEB',
#  'Content-Type': 'application/json'
#}
#conn.request("POST", "rest/secure/angelbroking/market/v1/quote/", payload, headers)
#res = conn.getresponse()
