#!/bin/julia
  
h2          = figure(num=2)
ax2         = h2.subplots()

ifall       = true      # Plot all

if ifall
  for k in 1:nv
    local rts         = uroots[:,k]
    local ufine       = LinRange(minimum(rts),maximum(rts),1000)
    local pot         = cubic_potential_func(ufine,rts)
    local pot_k       = ax2.plot(ufine,pot,linestyle="--",color=cm2(k-1),linewidth=2,label="v=$(vpar[k])")
    local F_k         = ax1.plot([ufine[1]; ufine[end]],[vpar[k]; vpar[k]],linestyle="--",color=cm2(k-1))
  end
  ax2.set_xlabel(L"u", fontsize=lafs)
  ax2.set_ylabel(L"ϕ(u)", fontsize=lafs)
  legend()
 
else
  k           = 3
  rts         = uroots[:,k]
  ufine       = LinRange(minimum(rts),maximum(rts),1000)
  pot         = cubic_potential_func(ufine,rts)
  pot_k       = ax2.plot(ufine,pot,linestyle="--",color=cm2(k-1),linewidth=2,label="v=$(vpar[k])")
  F_k         = ax1.plot([ufine[1]; ufine[end]],[vpar[k]; vpar[k]],linestyle="--",color=cm2(k-1))
  
  ax2.set_xlabel(L"u", fontsize=lafs)
  ax2.set_ylabel(L"ϕ(u)", fontsize=lafs)
end  

println("Done")









