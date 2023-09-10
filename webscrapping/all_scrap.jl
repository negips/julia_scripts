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

r,c         = size(df)
Links       = Array{String}(undef,r)
Names       = Array{String}(undef,r)
Plan_Type   = Array{String}(undef,r)
Fund_Type   = Array{String}(undef,r)
Crisil_Rank = Vector{Float64}(undef,r)
AuM         = Vector{Float64}(undef,r)
Gains_1w    = Vector{Float64}(undef,r)
Gains_1m    = Vector{Float64}(undef,r)
Gains_3m    = Vector{Float64}(undef,r)
Gains_6m    = Vector{Float64}(undef,r)
Gains_ytd   = Vector{Float64}(undef,r)
Gains_1y    = Vector{Float64}(undef,r)
Gains_2y    = Vector{Float64}(undef,r)
Gains_3y    = Vector{Float64}(undef,r)
Gains_5y    = Vector{Float64}(undef,r)
Gains_10y   = Vector{Float64}(undef,r)

for i in 1:r
  j = 1
  ch = df[i,j].children
  Links[i] = ch[1].attributes["href"]
  Names[i] = text(ch[1])
  #
  j = j + 1
  ch = df[i,j].children
  Plan_Type[i] = text(ch[1])
  j = j + 1
  ch = df[i,j].children
  Fund_Type[i] = text(ch[1])
  #
  j = j + 1
  ch = df[i,j].children
  cr = 0.0
  try 
    cr             = parse(Float64,text(ch[1]))
  catch e
    cr             = 0.0
  end
  Crisil_Rank[i] = cr
  #
  j = j + 1
  ch = df[i,j].children
  aum1 = text(ch[1])
  aum2 = replace(aum1,","=>"")
  AuM[i] = parse(Float64,aum2)
  #
  j = j + 1
  ch = df[i,j].children
  gains1 = text(ch[1])
  gains2 = replace(gains1,"%"=>"")
  g      = 0.0
  try 
    g             = parse(Float64,gains2)
  catch e
    g             = 0.0
  end
  Gains_1w[i] = g
  #
  j = j + 1
  ch = df[i,j].children
  gains1 = text(ch[1])
  gains2 = replace(gains1,"%"=>"")
  g      = 0.0
  try 
    g             = parse(Float64,gains2)
  catch e
    g             = 0.0
  end
  Gains_1m[i] = g
  #
  j = j + 1
  ch = df[i,j].children
  gains1 = text(ch[1])
  gains2 = replace(gains1,"%"=>"")
  g      = 0.0
  try 
    g             = parse(Float64,gains2)
  catch e
    g             = 0.0
  end
  Gains_3m[i] = g
  #
  j = j + 1
  ch = df[i,j].children
  gains1 = text(ch[1])
  gains2 = replace(gains1,"%"=>"")
  g      = 0.0
  try 
    g             = parse(Float64,gains2)
  catch e
    g             = 0.0
  end
  Gains_6m[i] = g
  #
  j = j + 1
  ch = df[i,j].children
  gains1 = text(ch[1])
  gains2 = replace(gains1,"%"=>"")
  g      = 0.0
  try 
    g             = parse(Float64,gains2)
  catch e
    g             = 0.0
  end
  Gains_1y[i] = g
  #
  j = j + 1
  ch = df[i,j].children
  gains1 = text(ch[1])
  gains2 = replace(gains1,"%"=>"")
  g      = 0.0
  try 
    g             = parse(Float64,gains2)
  catch e
    g             = 0.0
  end
  Gains_2y[i] = g
  #
  j = j + 1
  ch = df[i,j].children
  gains1 = text(ch[1])
  gains2 = replace(gains1,"%"=>"")
  g      = 0.0
  try 
    g             = parse(Float64,gains2)
  catch e
    g             = 0.0
  end
  Gains_3y[i] = g
  #
  j = j + 1
  ch = df[i,j].children
  gains1 = text(ch[1])
  gains2 = replace(gains1,"%"=>"")
  g      = 0.0
  try 
    g             = parse(Float64,gains2)
  catch e
    g             = 0.0
  end
  Gains_5y[i] = g
  #
  j = j + 1
  ch = df[i,j].children
  gains1 = text(ch[1])
  gains2 = replace(gains1,"%"=>"")
  g      = 0.0
  try 
    g             = parse(Float64,gains2)
  catch e
    g             = 0.0
  end
  Gains_10y[i] = g


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

tbls2_0 = scrape_tables(portfolio_1)
df2_0  = DataFrame(tbls2_0[5])

tbls2  = scrape_tables(portfolio_1,identity)
df2    = DataFrame(tbls2[5])

r2,c2             = size(df2)
Company_Link      = Array{String}(undef,r2)
Company_Name      = Array{String}(undef,r2)
#Plan_Type         = Array{String}(undef,r)
#Fund_Type         = Array{String}(undef,r)
#Crisil_Rank       = Vector{Float64}(undef,r)
#AuM               = Vector{Float64}(undef,r)
#Gains_1w          = Vector{Float64}(undef,r)
#Gains_1m          = Vector{Float64}(undef,r)
#Gains_3m          = Vector{Float64}(undef,r)
#Gains_6m          = Vector{Float64}(undef,r)
#Gains_ytd         = Vector{Float64}(undef,r)
#Gains_1y          = Vector{Float64}(undef,r)
#Gains_2y          = Vector{Float64}(undef,r)
#Gains_3y          = Vector{Float64}(undef,r)
#Gains_5y          = Vector{Float64}(undef,r)
#Gains_10y         = Vector{Float64}(undef,r)

for i in 1:r2
  j = 1
  l = length(df2[i,j].children)
  ch = df2[i,j].children[l].children
  Company_Link[i] = ch[1].attributes["href"]
  Company_Name[i]  = text(ch[1])
  #
  println(i)
end  

# Seems to work






















