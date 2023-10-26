import HTTP


host  = "https://apiconnect.angelbroking.com"
apikey = get(ENV, "ANGELONE_MARKET_API_KEY", "")

headers = Dict("X-PrivateKey" => apikey, 
               "Accept" => "application/json", 
               "X-SourceID" => "WEB", 
               "X-ClientLocalIP" => "CLIENT_LOCAL_IP", 
               "X-ClientPublicIP" => "CLIENT_PUBLIC_IP", 
               "X-MACAddress" => "MAC_ADDRESS", 
               "X-UserType" => "USER", 
               "Authorization" => "Bearer AUTHORIZATION_TOKEN", 
               "Accept" => "application/json", 
               "X-SourceID" => "WEB", 
               "Content-Type" => "application/json"
              )

payload = Dict("mode" => "FULL", 
               "exchangeTokens" => Dict("NSE" => ["3045"])
              )


endpoint = "rest/secure/angelbroking/market/v1/quote/"

url = "$host/$endpoint"

res = HTTP.get(url, headers, payload);

#res = getresponse(conn)

#data = read(res)

#println(join([data.decode("utf-8")], " "));
