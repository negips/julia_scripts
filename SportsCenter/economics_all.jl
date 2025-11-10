# Economics of Sports Center.

using Printf
using PyPlot

include("structures.jl")
#---------------------------------------------------------------------- 

const OneLakh           = 1.0e5
const OneCrore          = 1.0e7
ifverbose               = true
lafs                    = 12
lgfs                    = 12

Land_Investment         = 10.0e7
Construction_Investment = 10.0e7 
Total_Investment        = Land_Investment + Construction_Investment
Tax_Rate                = 0.34
Interest_Rate           = 0.05


Start_Year              = 2
End_Year                = 25
Total_Years             = End_Year - Start_Year + 1
Years                   = LinRange(Start_Year,End_Year,Total_Years)

Req_Yearly_Profit       = zeros(Float64,Total_Years)
Req_Monthly_Profit      = zeros(Float64,Total_Years)

for i in 1:Total_Years
  
  global Req_Yearly_Profit
  global Req_Monthly_Profit
  yr  = Years[i]
  Req_Yearly_Profit[i]   = Total_Investment*(1.0 + Interest_Rate)/yr/OneCrore
  Req_Monthly_Profit[i]  = Total_Investment*(1.0 + Interest_Rate)/yr/12.0/OneLakh

end  

close("all")
cm  = get_cmap("tab10")

h1  = figure(num=1,figsize=[10.0,7.0])
ax1 = gca()
ax1.plot(Years,Req_Monthly_Profit,linestyle="-",marker="o",markersize=6,color=cm(0))
xticks = Start_Year:2:End_Year
ax1.set_xticks(xticks)
ax1.grid()
ax1.set_ylabel("Required Monthly Profit (Lakhs)")
ax1.set_xlabel("Years")

# ax2 = ax1.twinx()
# ax2.plot(Years,Req_Yearly_Profit,linestyle="--",marker="s",markersize=4,color=cm(1))
# ax2.set_ylabel("Yearly Profit (Cr.)")

#---------------------------------------------------------------------- 

# Initialize sports,prices,slots, etc.
include("configuration_init.jl")

Depreciation                  = 0.125*Land_Investment
Monthly_Operating_Expenses    = 8.0*OneLakh

# Quater-capacity
#---------------------------------------------------------------------- 
include("configuration_quarter_capacity.jl")
QuaterCapEconomics = GetCenterEconomics(All_sports,Tax_Rate,Depreciation,Monthly_Operating_Expenses)
qcap_pr            = fill(QuaterCapEconomics.MonthlyProfit,Total_Years)./OneLakh
ax1.plot(Years,qcap_pr,linestyle="--",color=cm(2))

# Half-capacity
#---------------------------------------------------------------------- 
include("configuration_half_capacity.jl")
HalfCapEconomics  = GetCenterEconomics(All_sports,Tax_Rate,Depreciation,Monthly_Operating_Expenses)
hcap_pr           = fill(HalfCapEconomics.MonthlyProfit,Total_Years)./OneLakh
ax1.plot(Years,hcap_pr,linestyle="--",color=cm(3))

# ThreeFourths-capacity
#---------------------------------------------------------------------- 
include("configuration_threefourths_capacity.jl")
TFCapEconomics    = GetCenterEconomics(All_sports,Tax_Rate,Depreciation,Monthly_Operating_Expenses)
tfcap_pr          = fill(TFCapEconomics.MonthlyProfit,Total_Years)./OneLakh
ax1.plot(Years,tfcap_pr,linestyle="--",color=cm(4))

# Full-capacity
#---------------------------------------------------------------------- 
include("configuration_full_capacity.jl")
FCapEconomics     = GetCenterEconomics(All_sports,Tax_Rate,Depreciation,Monthly_Operating_Expenses)
fcap_pr           = fill(FCapEconomics.MonthlyProfit,Total_Years)./OneLakh
ax1.plot(Years,fcap_pr,linestyle="--",color=cm(5))
#---------------------------------------------------------------------- 


# Parametric with occupancy

Occup              = 0.2:0.05:1.0
nOccup             = length(Occup)
OccupancyEconomics = Vector{CenterEconomics}(undef,nOccup)
PaybackPeriod      = zeros(Float64,nOccup)

for i in 1:nOccup
  occ = Occup[i]
  SetAllOccupancy!(All_sports,occ)
  # Remove Basketball and Volleyball
  SetOccupancy!(All_sports[6],0.0)
  SetOccupancy!(All_sports[7],0.0)

  OccupancyEconomics[i] = GetCenterEconomics(All_sports,Tax_Rate,Depreciation,Monthly_Operating_Expenses)
  PaybackPeriod[i]      = Total_Investment*(1.0 + Interest_Rate)/OccupancyEconomics[i].YearlyProfit
end  

h2  = figure(num=2,figsize=[10.0,7.0])
ax3 = gca()
ax3.plot(Occup,PaybackPeriod,linestyle="-",marker="o",markersize=6,color=cm(0))
ax3.grid()
ax3.set_ylabel("Time-Adjusted Payback Period (Years)",fontsize=lafs)
ax3.set_xlabel("Occupancy",fontsize=lafs)
heading = @sprintf("Total Investment: %.1f Cr.", Total_Investment/OneCrore)
ax3.set_title(heading)

h2.savefig("PaybackPeriod_Occupancy.png")

#---------------------------------------------------------------------- 
println("Done.")













