# Economics of Sports Center.

using Printf

const OneLakh           = 1.0e5
const OneCrore          = 1.0e7

Land_Investment         = 4.0e7
Construction_Investment = 10.0e7 
Total_Investment        = Land_Investment + Construction_Investment

Depreciation            = 0.125*Land_Investment
Tax_Rate                = 0.34
Interest_Rate           = 0.05

Break_Even              = 5.0            # years

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
@printf "----------------------------------------\n" 

#@printf "Monthly Revenue %.3f Lakh.\n" Total_Monthly_Revenue/OneLakh
#@printf "Monthly Profit %.3f Lakh.\n" Monthly_Profit/OneLakh

#---------------------------------------------------------------------- 
println("Done.")













