VERSION >= v"1.6.0"

module MarketStack

#using Compat # For compatibility across Versions
using ArgCheck
using DelimitedFiles
using HTTP
using JSON


MSEndPoints = ["eod", "intraday", "splits", "dividends", "tickers", "exchanges", "currencies", "timezones"];
StockExchanges    = Dict( "XASE"    =>    "American Stock Exchange"           ,
                          "XASX"    =>    "Australian Stock Exchange"         ,
                          "XLIM"    =>    "Bolsa de Valores de Lima"          ,
                          "BVFM"    =>    "B3 - Brasil Bolsa Balcao S.A"      ,
                          "XBAH"    =>    "Bahrein Bourse"                    ,
                          "XBEL"    =>    "Belgrade Stock Exchange"           ,
                          "XBOG"    =>    "Bolsa de Valores de Colombia"      ,
                          "BMEX"    =>    "Bolsas y Marcados Espanoles"       ,
                          "XBOM"    =>    "Bombay Stock Exchange"             ,
                          "XMIL"    =>    "Borsa Italiana"                    ,
                          "XBUE"    =>    "Buenos Aires Stock Exchange"       ,
                          "XTSU"    =>    "Borse Stuttgart"                   ,
                          "XCNQ"    =>    "Canadian Securities Exchange"      ,
                          "XCSE"    =>    "Copenhagen Stock Exchange"         ,
                          "XFRA"    =>    "Deutsche Borse"                    ,
                          "XETRA"   =>    "Deutsche Borse Xetra"              ,
                          "XDFM"    =>    "Dubai Financial Market"            ,
                          "XCAI"    =>    "Egyptian Exchange"                 ,
                          "XAMS"    =>    "Euronext Amesterdam"               ,
                          "XBRU"    =>    "Euronext Brussels"                 ,
                          "XLIS"    =>    "Euronext Lison"                    ,
                          "XPAR"    =>    "Euronext Paris"                    ,
                          "XFKA"    =>    "Fukuoka Stock Exchange"            ,
                          "XHEL"    =>    "Helsinki Stock Exchange"           ,
                          "XSTC"    =>    "Ho Chi Minh Stock Exchange"        ,
                          "XHKG"    =>    "Hong Kong Stock Exchange"          ,
                          "IEXG"    =>    "Investors Exchange"                ,
                          "XIST"    =>    "Istanbul Stock Exchange"           ,
                          "XIDX"    =>    "Jakarta Stock Exchange"            ,
                          "XJSE"    =>    "Johannesburg Stock Exchange"       ,
                          "XKRX"    =>    "Korean Stock Exchange"             ,
                          "XLON"    =>    "London Stock Exchange"             ,
                          "XKLS"    =>    "Malaysia Stock Exchange"           ,
                          "XMEX"    =>    "Mexican Stock Exchange"            ,
                          "MISX"    =>    "Moscow Stock Exchange"             ,
                          "XNAS"    =>    "NASDAQ Stock Exchange"             ,
                          "ARCX"    =>    "NYSE ARCA"                         ,
                          "XNGO"    =>    "Nagoya Stock Exchange"             ,
                          "XICE"    =>    "Nasdaq Island"                     ,
                          "XRIS"    =>    "Nasdaq Riga"                       ,
                          "XLIT"    =>    "Nasdaq Vilnius"                    ,
                          "XNSE"    =>    "National Stock Exchange India"     ,
                          "XNYS"    =>    "New York Stock Exchange"           ,
                          "XNZE"    =>    "New Zealand Stock Exchange"        ,
                          "XNSA"    =>    "Nigerian Stock Exchange"           ,
                          "OOTC"    =>    "OTC Bulletin Board"                ,
                          "PSGM"    =>    "OTC Grey Market"                   ,
                          "OTCM"    =>    "OTC Markets"                       ,
                          "PINC"    =>    "OTC PINK Current"                  ,
                          "OTCQB"   =>    "OTCQB Marketplace"                 ,
                          "OTCQX"   =>    "OTCQX Marketplace"                 ,
                          "XOSL"    =>    "Oslo Stock Exchange"               ,
                          "DSMD"    =>    "Qatar Stock Exchange"              ,
                          "XSWX"    =>    "SIX Swiss Exchange"                ,
                          "XSGO"    =>    "Santiago Stock Exchange"           ,
                          "XSAP"    =>    "Sapporo Stock Exchange"            ,
                          "XSAU"    =>    "Saudi Stock Exchange"              ,
                          "XSHG"    =>    "Shanghai Stock Exchange"           ,
                          "XSHE"    =>    "Shenzhen Stock Exchange"           ,
                          "XSES"    =>    "Singapore Stock Exchange"          ,
                          "XBKK"    =>    "Stock Exchange of Thailand"        ,
                          "XSTO"    =>    "Stockholm Stock Exchange"          ,
                          "XTSX"    =>    "TSX Venture Exchange"              ,
                          "XTAI"    =>    "Taiwan Stock Exchange"             ,
                          "XTAL"    =>    "Tallinn Stock Exchange"            ,
                          "XTAE"    =>    "Tel Aviv Stock Exchange"           ,
                          "XTKS"    =>    "Tokyo Stock Exchange"              ,
                          "XTSE"    =>    "Toronto Stock Exchange"            ,
                          "NMFQS"   =>    "US Mutual Funds"                   ,
                          "XWAR"    =>    "Warsaw Stock Exchange"             ,
                          )

ExchangeCodes = String.(keys(StockExchanges))
ExchangeNames = String.(values(StockExchanges))

include("msclient.jl")
include("msresponse.jl")
include("utils.jl")
include("stock_time_series.jl")
#include("digital_currency.jl")
#include("foreign_exchange_currency.jl")
#include("stock_technical_indicators.jl")
#include("sector_performance.jl")
#include("fundamental_values.jl")
#include("fundamentals.jl")
#include("economic_indicators.jl")

# Nothing to Export Yet
export MS_stock_quote
export StockExchanges

## avclient
#export key, AlphaVantageClient, AlphaVantageResponse
#
## stock_time_series
#export
#    time_series_intraday,
#    time_series_intraday_extended
## `time_series_daily` etc are exported in macro
#
## digital_currency
#export crypto_rating, digital_currency_intraday
## `digital_currency_daily` etc are exported in macro
#
## foreign_exchange_currency
#export
#    currency_exchange_rate,
#    fx_intraday,
#    fx_daily
## `fx_weekly` etc are exported in macro
#
## stock_technical_indicators
## `VWAP` etc are exported in macro
#
## sector_performance
#export sector_performance
#
## fundamentals
#export
#    company_overview,
#    income_statement,
#    balance_sheet,
#    cash_flow,
#    listing_status,
#    earnings,
#    earnings_calendar,
#    ipo_calendar
#
## economic_indicators
## others exported in macro
#export
#    real_gdp,
#    treasury_yield,
#    federal_fund_rate,
#    cpi

end # module
