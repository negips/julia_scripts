mutable struct MarketStackClient
    scheme::String
    key::String
    host::String
end

MSClient(; scheme = "http", key = "", host = marketstack_api) = MarketStackClient(scheme, key, host)
#
const MSGlobal = Ref(MSClient(key = get(ENV, "MARKETSTACK_API_KEY", ""), host = get(ENV, "MARKETSTACK_HOST", "api.marketstack.com/v1")))
#
function key(client::MarketStackClient)
    if isempty(client.key)
        @warn "No API key found"
    end
    return client.key
end
