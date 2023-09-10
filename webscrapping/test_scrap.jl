!#/bin/julia

println("Testing Website Scrapping")

using HTTP
using Gumbo
using AbstractTrees
using Cascadia
using TableScraper
using DataFrames, DataFramesMeta

url = "https://www.moneycontrol.com/mutual-funds/performance-tracker/returns/small-cap-fund.html"

mfpage            = HTTP.get(url)
mfpage_parsed     = parsehtml(String(mfpage.body))
root              = mfpage_parsed.root

tbls = scrape_tables(url,identity)
df   = DataFrame(tbls[2])

r,c = size(df)
Links = Array{String}(undef,r)

for i in 1:r
  ch = df[i,1].children
  Links[i] = ch[1].attributes["href"]
end  

#display(Links)
l1 = Links[1]
l2 = replace(l1,"nav/"=>"")
il = findlast("/",l2)
i1 = il.start
ll = length(l2)
id = l2[i1+1:ll]

#portfolio_1 = l2[1:i1]*"portfolio-overview/"*id
portfolio_1 = l2[1:i1]*"portfolio-holdings/"*id

tbls_2 = scrape_tables(portfolio_1)
df2    = DataFrame(tbls_2[5])

# Seems to work



