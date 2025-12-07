# Economics of Sports Center.

occupancy               = 0.40
for i in 1:nSports
  SetOccupancy!(All_sports[i],occupancy)
end  
# Remove Basketball and Volleyball
SetOccupancy!(All_sports[6],0.0)
SetOccupancy!(All_sports[7],0.0)


# i                       = 1         # 1
# name                    = "Badminton"
# occupancy               = 0.40
# SetOccupancy!(All_sports[i],occupancy)
# # All_sports[1]           = DailySport(name,courts,slots,occupancy,price) 
# 
# i                       = i+1       # 2
# name                    = "Tennis"
# occupancy               = 0.40
# SetOccupancy!(All_sports[i],occupancy)
# 
# i                       = i+1       # 3
# name                    = "TableTennis"
# occupancy               = 0.40
# SetOccupancy!(All_sports[i],occupancy)
# 
# i                       = i+1       # 4
# name                    = "Squash"
# occupancy               = 0.40
# SetOccupancy!(All_sports[i],occupancy)
# 
# i                       = i+1       # 5
# name                    = "Padel"
# occupancy               = 0.40
# SetOccupancy!(All_sports[i],occupancy)
# 
# i                       = i+1       # 6
# name                    = "Basketball"
# occupancy               = 0.00
# SetOccupancy!(All_sports[i],occupancy)
# 
# i                       = i+1       # 7
# name                    = "Volleyball"
# occupancy               = 0.00
# SetOccupancy!(All_sports[i],occupancy)
# 
# i                       = i+1       # 8
# name                    = "Skating"
# occupancy               = 0.40
# SetOccupancy!(All_sports[i],occupancy)
# 
# i                       = i+1       # 9
# name                    = "Studio"
# occupancy               = 0.40
# SetOccupancy!(All_sports[i],occupancy)
# 
# i                       = i+1       # 10
# name                    = "Gym"
# occupancy               = 0.40
# SetOccupancy!(All_sports[i],occupancy)
# 
# i                       = i+1       # 11
# name                    = "Futsal"
# occupancy               = 0.40
# SetOccupancy!(All_sports[i],occupancy)
# 
# i                       = i+1       # 12
# name                    = "Snooker"
# occupancy               = 0.40
# SetOccupancy!(All_sports[i],occupancy)













