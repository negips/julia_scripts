# Economics of Sports Center.

using Printf


mutable struct Sport

  Name::String
  Courts::Int
  Slots::Int
  Occupancy::Float64
  Price::Float64
  Daily_Collection::Float64
  Monthly_Collection::Float64

end   
#---------------------------------------------------------------------- 
function Sport(name::String,courts::Int,slots::Int,occupancy::Float64,price::Float64)

  daily_court_hours     = courts*slots
  daily_collection      = daily_court_hours*occupancy*price
  monthly_collection    = daily_collection*30

  sport = Sport(name,courts,slots,occupancy,price,daily_collection,monthly_collection)

  return sport
end  
#---------------------------------------------------------------------- 

include("configuration_init.jl")
include("configuration2.jl")


Total_Monthly_Revenue = 0.0
Total_Daily_Revenue   = 0.0
for sport in All_sports
  global Total_Monthly_Revenue, Total_Daily_Revenue
  Total_Monthly_Revenue = Total_Monthly_Revenue + sport.Monthly_Collection
  Total_Daily_Revenue   = Total_Daily_Revenue + sport.Daily_Collection
end
Total_Yearly_Revenue    = Total_Monthly_Revenue*12


# 
Land_Investment         = 7.5e7
Construction_Investment = Land_Investment
Total_Investment        = Land_Investment + Construction_Investment

Depreciation            = 0.125*Land_Investment
Tax_rate                = 0.34
Operating_Expenses      = 0.5*Total_Yearly_Revenue

Taxable_Income          = Total_Yearly_Revenue - Operating_Expenses
Yearly_Tax              = Tax_rate*Taxable_Income

Yearly_Profit           = Total_Yearly_Revenue - Operating_Expenses - Yearly_Tax
Monthly_Profit          = Yearly_Profit/12.0 

const OneLakh           = 1.0e5
const OneCrore          = 1.0e7

@printf "Yearly Revenue %.3f Cr.\n" Total_Yearly_Revenue/OneCrore
@printf "Yearly Operating Expenses %.3f Cr.\n" Operating_Expenses/OneCrore
@printf "Yearly Tax %.3f Cr.\n" Yearly_Tax/OneCrore
@printf "Yearly Profit  %.3f Cr.\n" Yearly_Profit/OneCrore
@printf "----------------------------------------\n" 
@printf "Monthly Revenue %.3f Lakh.\n" Total_Monthly_Revenue/OneLakh
@printf "Monthly Profit %.3f Lakh.\n" Monthly_Profit/OneLakh

#---------------------------------------------------------------------- 
println("Done.")













