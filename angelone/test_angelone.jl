println("Testing PyCall for AngelOne")

using PyCall
using JSON
using HTTP

#angel = pyimport("SmartApi")

payload = Dict("clientcode" => "CLIENT_ID",
               "password" => "CLIENT_PIN",
		   "totp" => "TOTP_CODE")

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



