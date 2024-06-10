println("Yinyun - Neuron Model")

using LinearAlgebra


# C         - constant capacitance
# V         - voltage
# g         - Channel conductance
# m         - dynamic variables (probability of opening and closing of Sodium channel) ϵ [0,1]
# h         - dynamic variables (probability of inactivation of sodium channel) ϵ [0,1]
# n         - dynamic variables (probability of opening and closing of Potassium channel) ϵ [0,1]
# α         - opening/closing voltage dependence
# GABA      - conductance
