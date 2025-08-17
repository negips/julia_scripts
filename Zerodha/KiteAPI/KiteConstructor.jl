#     Add Definition of Constructors here
#---------------------------------------------------------------------- 
"""
      function KiteConnection(field::Array)
        
      Constructor for the KiteConnection

"""
function KiteConnection(key::String,secret::String)


  rtoken = ""
  atoken = ""
  return KiteConnection(key,secret,rtoken,atoken)
end        

#---------------------------------------------------------------------- 

