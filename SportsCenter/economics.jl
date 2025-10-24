# Economics of Sports Center.

using Printf


mutable struct DailySport

  Name::String
  Courts::Int
  Slots::Int
  Occupancy::Float64
  Price::Float64
  Daily_Collection::Float64
  MonthlyCollection::Float64

end   
#---------------------------------------------------------------------- 
mutable struct MemberSport

  Name::String
  Capacity::Int
  Slots::Int
  Occupancy::Float64
  Price::Float64
  Daily_Collection::Float64
  MonthlyCollection::Float64

end   
#---------------------------------------------------------------------- 
function DailySport(name::String,courts::Int,slots::Int,occupancy::Float64,price::Float64)

  daily_court_hours     = courts*slots
  daily_collection      = daily_court_hours*occupancy*price
  monthly_collection    = daily_collection*30

  sport = DailySport(name,courts,slots,occupancy,price,daily_collection,monthly_collection)

  return sport
end  
#---------------------------------------------------------------------- 
function SetPrice!(sport::DailySport,price::Float64)

  daily_court_hours           = sport.Courts*sport.Slots
  occupancy                   = sport.Occupancy
  daily_collection            = daily_court_hours*occupancy*price
  monthly_collection          = daily_collection*30

  sport.Price                 = price
  sport.Daily_Collection      = daily_collection
  sport.MonthlyCollection    = monthly_collection

  return nothing
end  
#---------------------------------------------------------------------- 
function SetOccupancy!(sport::DailySport,occupancy::Float64)

  price                       = sport.Price
  daily_court_hours           = sport.Courts*sport.Slots

  daily_collection            = daily_court_hours*occupancy*price
  monthly_collection          = daily_collection*30

  sport.Occupancy             = occupancy
  sport.Daily_Collection      = daily_collection
  sport.MonthlyCollection     = monthly_collection

  return nothing
end  
#---------------------------------------------------------------------- 

const OneLakh           = 1.0e5
const OneCrore          = 1.0e7
ifverbose               = true

include("configuration_init.jl")
#include("configuration4.jl")
#include("configuration_quarter_capacity.jl")
#include("configuration_half_capacity.jl")
include("configuration_threefourths_capacity.jl")
#include("configuration_full_capacity.jl")

Total_Monthly_Revenue   = 0.0
Total_Daily_Revenue     = 0.0
for sport in All_sports
  global Total_Monthly_Revenue, Total_Daily_Revenue
  Total_Monthly_Revenue = Total_Monthly_Revenue + sport.MonthlyCollection
  Total_Daily_Revenue   = Total_Daily_Revenue + sport.Daily_Collection

  if (ifverbose)
    revenue       = sport.MonthlyCollection
    occup         = sport.Occupancy
    @printf("%12s, at Occupancy: %.2f, Monthly Collection: %.3f Lakh\n", sport.Name, occup, revenue/OneLakh)
  end
end
Total_Yearly_Revenue          = Total_Monthly_Revenue*12

# 
Land_Investment               = 5.0*OneCrore
Construction_Investment       = 10.0*OneCrore
Total_Investment              = Land_Investment + Construction_Investment

Depreciation                  = 0.125*Land_Investment
Tax_rate                      = 0.34
Monthly_Operating_Expenses    = 8.0*OneLakh
Yearly_Operating_Expenses     = Monthly_Operating_Expenses*12

Taxable_Income          = Total_Yearly_Revenue - Yearly_Operating_Expenses - Depreciation
Yearly_Tax              = Tax_rate*Taxable_Income

Yearly_Profit           = Total_Yearly_Revenue - Yearly_Operating_Expenses - Yearly_Tax
Monthly_Profit          = Yearly_Profit/12.0 

@printf "----------------------------------------\n" 
@printf "Yearly Revenue %.3f Cr.\n" Total_Yearly_Revenue/OneCrore
@printf "Yearly Operating Expenses %.3f Cr.\n" Yearly_Operating_Expenses/OneCrore
@printf "Yearly Tax %.3f Cr.\n" Yearly_Tax/OneCrore
@printf "Yearly Profit  %.3f Cr.\n" Yearly_Profit/OneCrore
@printf "----------------------------------------\n" 
@printf "Monthly Revenue %.3f Lakh.\n" Total_Monthly_Revenue/OneLakh
@printf "Monthly Profit %.3f Lakh.\n\n" Monthly_Profit/OneLakh
#---------------------------------------------------------------------- 
println("Done.")













