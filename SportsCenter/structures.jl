# Economics of Sports Center.

mutable struct DailySport

  Name::String
  Courts::Int
  Slots::Int
  Occupancy::Float64
  Price::Float64
  DailyRevenue::Float64
  MonthlyRevenue::Float64

end   
#---------------------------------------------------------------------- 
mutable struct MemberSport

  Name::String
  Capacity::Int
  Slots::Int
  Occupancy::Float64
  Price::Float64
  Daily_Collection::Float64
  MonthlyRevenue::Float64

end   
#---------------------------------------------------------------------- 
mutable struct CenterEconomics

  Sports::Vector{DailySport}
  DailyRevenue::Float64
  MonthlyRevenue::Float64
  YearlyyRevenue::Float64
  DailyOpExpenses::Float64
  MonthlyOpExpenses::Float64
  YearlyOpExpenses::Float64
  Depreciation::Float64
  TaxRate::Float64
  MonthlyProfit::Float64
  YearlyProfit::Float64

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
  sport.DailyRevenue          = daily_collection
  sport.MonthlyRevenue        = monthly_collection

  return nothing
end  
#---------------------------------------------------------------------- 
function SetOccupancy!(sport::DailySport,occupancy::Float64)

  price                       = sport.Price
  daily_court_hours           = sport.Courts*sport.Slots

  daily_collection            = daily_court_hours*occupancy*price
  monthly_collection          = daily_collection*30

  sport.Occupancy             = occupancy
  sport.DailyRevenue          = daily_collection
  sport.MonthlyRevenue        = monthly_collection

  return nothing
end  
#---------------------------------------------------------------------- 
function SetAllOccupancy!(Sports::Vector{DailySport},occupancy::Float64)

  ns = length(Sports)
  for i in 1:ns
    SetOccupancy!(Sports[i],occupancy)
  end  

  return nothing
end  
#---------------------------------------------------------------------- 
function GetCenterEconomics(AllSports::Vector{DailySport},taxrate,depreciation,MonthlyOpExp)

      TotalMonthlyRevenue     = 0.0
      TotalDailyRevenue       = 0.0
      for sport in AllSports
        TotalMonthlyRevenue   = TotalMonthlyRevenue + sport.MonthlyRevenue
        TotalDailyRevenue     = TotalDailyRevenue + sport.DailyRevenue
      end
      TotalYearlyRevenue      = TotalMonthlyRevenue*12
      
      YearlyOpExp             = MonthlyOpExp*12
      DailyOpExp              = MonthlyOpExp/30
      TaxableIncome           = TotalYearlyRevenue - YearlyOpExp - depreciation
      YearlyTax               = taxrate*TaxableIncome
     
      YearlyProfit            = TotalYearlyRevenue - YearlyOpExp - YearlyTax
      MonthlyProfit           = YearlyProfit/12.0 

      Economics               = CenterEconomics(AllSports,TotalDailyRevenue,TotalMonthlyRevenue,TotalYearlyRevenue,DailyOpExp,MonthlyOpExp,YearlyOpExp,depreciation,taxrate,MonthlyProfit,YearlyProfit)

      return Economics

end
#---------------------------------------------------------------------- 
function PrintRevenue(Sports::Vector{DailySport})

  One_Lakh         = 1.0e5
  for sport in Sports
    revenue       = sport.MonthlyRevenue
    occup         = sport.Occupancy
    @printf("%12s, at Occupancy: %.2f, Monthly Collection: %.3f Lakh\n", sport.Name, occup, revenue/One_Lakh)
  end
end  
#---------------------------------------------------------------------- 




