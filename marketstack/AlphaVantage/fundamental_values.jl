COMPANY_OVERVIEW_KEYS =  ["SharesOutstanding",
"ExDividendDate",
"52WeekLow",
"ReturnOnEquityTTM",
"LatestQuarter",
"200DayMovingAverage",
"EVToEBITDA",
"RevenuePerShareTTM",
"Beta",
"Sector", 
"ForwardAnnualDividendYield",
"Exchange",
"PercentInsiders",
"QuarterlyEarningsGrowthYOY",
"Currency",
"EBITDA",
"ShortRatio",
"DividendYield",
"AnalystTargetPrice",
"DilutedEPSTTM",
"BookValue",
"LastSplitDate",
"SharesFloat",
"PriceToSalesRatioTTM",
"FullTimeEmployees",
"ShortPercentOutstanding",
"52WeekHigh",
"ReturnOnAssetsTTM",
"PayoutRatio",
"PriceToBookRatio",
"Symbol",
"QuarterlyRevenueGrowthYOY",
"DividendDate",
"ProfitMargin",
"50DayMovingAverage",
"LastSplitFactor",
"RevenueTTM",
"PEGRatio",
"OperatingMarginTTM",
"Industry",
"ForwardPE",
"EVToRevenue",
"ForwardAnnualDividendRate",
"ShortPercentFloat",
"SharesShortPriorMonth",
"TrailingPE",
"SharesShort",
"Country",
"Address",
"EPS",
"GrossProfitTTM",
"Name",
"MarketCapitalization",
"PERatio",
"DividendPerShare",
"Description",
"FiscalYearEnd",
"AssetType",
"PercentInstitutions"]

COMMON_KEYS = ["reportedCurrency", "netIncome"]

COMMON = [[("quarterlyReports", "income_statement", x, "fiscalDateEnding") for x in COMMON_KEYS];
            [("annualReports", "income_statement", x, "fiscalDateEnding") for x in COMMON_KEYS]]

INCOME_STATEMENT_KEYS = [
"incomeTaxExpense",
"otherNonOperatingIncome",
"minorityInterest",
"discontinuedOperations",
"incomeBeforeTax",
"totalOtherIncomeExpense",
"interestIncome",
"researchAndDevelopment",
"grossProfit",
"totalRevenue",
"otherOperatingExpense",
"taxProvision",
"extraordinaryItems",
"ebit",
"otherItems",
"netIncomeApplicableToCommonShares",
"totalOperatingExpense",
"costOfRevenue",
"interestExpense",
"sellingGeneralAdministrative",
"operatingIncome",
"netIncomeFromContinuingOperations",
"netInterestIncome",
"effectOfAccountingCharges",
"nonRecurring",
"preferredStockAndOtherAdjustments"]

INCOME_STATEMENT = [[("quarterlyReports", "income_statement", x, "fiscalDateEnding") for x in INCOME_STATEMENT_KEYS];
                    [("annualReports", "income_statement", x, "fiscalDateEnding") for x in INCOME_STATEMENT_KEYS]]  


BALANCE_SHEET_KEYS = ["totalPermanentEquity",
"warrants",
"negativeGoodwill",
"preferredStockTotalEquity",
"accumulatedAmortization",
"inventory",
"additionalPaidInCapital",
"commonStockTotalEquity",
"longTermInvestments",
"netTangibleAssets",
"cashAndShortTermInvestments",
"longTermDebt",
"otherShareholderEquity",
"totalCurrentAssets",
"treasuryStock",
"otherAssets",
"capitalSurplus",
"totalNonCurrentAssets",
"accountsPayable",
"totalNonCurrentLiabilities",
"otherLiabilities",
"totalShareholderEquity",
"liabilitiesAndShareholderEquity",
"otherCurrentAssets",
"totalCurrentLiabilities",
"otherNonCurrrentAssets",
"shortTermDebt",
"commonStockSharesOutstanding",
"capitalLeaseObligations",
"netReceivables",
"retainedEarningsTotalEquity",
"earningAssets",
"totalAssets",
"commonStock",
"cash",
"deferredLongTermLiabilities",
"totalLongTermDebt",
"retainedEarnings",
"shortTermInvestments",
"propertyPlantEquipment",
"goodwill",
"preferredStockRedeemable",
"totalLiabilities",
"otherNonCurrentLiabilities",
"currentLongTermDebt",
"intangibleAssets",
"accumulatedDepreciation",
"otherCurrentLiabilities",
"deferredLongTermAssetCharges"]

BALANCE_SHEET = [[("quarterlyReports", "balance_sheet", x, "fiscalDateEnding") for x in BALANCE_SHEET_KEYS];
                [("annualReports", "balance_sheet", x, "fiscalDateEnding") for x in BALANCE_SHEET_KEYS]]  

 CASH_FLOW_KEYS = ["cashflowFromInvestment",
 "changeInInventory",
 "changeInAccountReceivables",
 "changeInCashAndCashEquivalents",
 "otherOperatingCashflow",
 "dividendPayout",
 "changeInReceivables",
 "capitalExpenditures",
 "changeInExchangeRate",
 "operatingCashflow",
 "cashflowFromFinancing",
 "changeInLiabilities",
 "stockSaleAndPurchase",
 "otherCashflowFromFinancing",
 "changeInOperatingActivities",
 "depreciation",
 "changeInCash",
 "netBorrowings",
 "investments",
 "changeInNetIncome",
 "otherCashflowFromInvestment"]

 CASH_FLOW = [[("quarterlyReports", "cash_flow", x, "fiscalDateEnding") for x in CASH_FLOW_KEYS];
              [("annualReports", "cash_flow", x, "fiscalDateEnding") for x in CASH_FLOW_KEYS]]  


 EARNINGS_KEYS_Q = ["reportedDate",
 "estimatedEPS",
 "surprise",
 "surprisePercentage",
 "reportedEPS"]

 EARNINGS_KEYS_A = ["reportedEPS"]

 EARNINGS = [[("quarterlyEarnings", "earnings", x, "reportedDate") for x in EARNINGS_KEYS_Q];
             [("annualEarnings", "earnings", x, "fiscalDateEnding") for x in EARNINGS_KEYS_A]]  

 FUNDAMENTAL_VALUES = vcat(INCOME_STATEMENT, BALANCE_SHEET, CASH_FLOW, EARNINGS)