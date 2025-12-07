# Economics of Sports Center.

using Printf
using PyPlot

const OneLakh           = 1.0e5
const OneCrore          = 1.0e7

Land_Investment         = 5.0e7
Construction_Investment = 10.0e7 
Total_Investment        = Land_Investment + Construction_Investment

Tax_Rate                = 0.34
Interest_Rate           = 0.05
Break_Even              = 20.0            # years

Investment_Value        = Total_Investment*(1.0 + Interest_Rate)^Break_Even
Required_Profit         = Total_Investment*(1.0 + Interest_Rate)/Break_Even

@printf "Doing a Simple Calculation on Total Investment %.4f Cr:\n" Total_Investment/OneCrore
@printf "Break Even %i years.\n" Break_Even
@printf "Assuming Interest Rate of %f %%\n" Interest_Rate
@printf "Year 1: Annual Profit: %.3f Cr.\n" Required_Profit/OneCrore
@printf "Year 1: Monthly Profit: %.3f Lakh.\n" Required_Profit/12/OneLakh
@printf "----------------------------------------\n" 


P1                      = (Land_Investment*(1.0 + Interest_Rate)^4)/Break_Even 
P2                      = (0.5*Construction_Investment*(1.0 + Interest_Rate)^2)/Break_Even
P3                      = (0.5*Construction_Investment*(1.0 + Interest_Rate))/Break_Even

Required_Profit2        = P1 + P2 + P3

@printf "Splitting Investment over 3 Years:\n"
@printf "Year 1: Annual Profit: %.3f Cr.\n" Required_Profit2/OneCrore
@printf "Year 1: Monthly Profit: %.3f Lakh.\n" Required_Profit2/12/OneLakh
@printf "----------------------------------------\n\n" 


Total_Years = 20
Req_Yearly_Profit  = zeros(Float64,Total_Years)
Req_Monthly_Profit = zeros(Float64,Total_Years)

for i in 1:Total_Years
  
  global Req_Yearly_Profit
  global Req_Monthly_Profit
  yr  = i
  Req_Yearly_Profit[i]   = Total_Investment*(1.0 + Interest_Rate)/yr/OneCrore
  Req_Monthly_Profit[i]  = Total_Investment*(1.0 + Interest_Rate)/yr/12.0/OneLakh

end  

close("all")
cm  = get_cmap("tab10")

h1  = figure(num=1,figsize=[9.0,7.0])
ax1 = gca()
ax1.plot(LinRange(1,Total_Years,Total_Years),Req_Monthly_Profit,linestyle="-",marker="o",markersize=6,color=cm(0))
xticks = Vector(0:2:20)
ax1.set_xticks(xticks)
ax1.grid()
ax1.set_ylabel("Yearly Profit (Lakhs)")
ax1.set_xlabel("Years")

ax2 = ax1.twinx()
ax2.plot(LinRange(1,Total_Years,Total_Years),Req_Yearly_Profit,linestyle="--",marker="s",markersize=4,color=cm(1))
ax2.set_ylabel("Monthly Profit (Cr.)")
#ax2.grid()

#@printf "Monthly Revenue %.3f Lakh.\n" Total_Monthly_Revenue/OneLakh
#@printf "Monthly Profit %.3f Lakh.\n" Monthly_Profit/OneLakh

#---------------------------------------------------------------------- 
println("Done.")













