#     Add Definition of Constructors here
#---------------------------------------------------------------------- 
"""
      function KiteConnection(key::String,secret::String)
        
      Create a new connection from the API key and API secret.

"""
function KiteConnection(key::String,secret::String)

  # Get new tokens
  rtoken,atoken = kite_get_tokens(key,secret)

  # Create Connection object
  return KiteConnection(key,secret,rtoken,atoken)
end        

#---------------------------------------------------------------------- 

