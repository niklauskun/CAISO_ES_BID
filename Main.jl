using JuMP
using Gurobi
using MAT
using Printf
using DataFrames
using CSV
# using Plots
using JLD

# High-level Settings
Zone = "WALNUT" # price zone name

# read price
fileln = matopen(string("./Data/",Zone,".mat"))
RTP = read(fileln, "Q")
close(fileln)

# simulation setting
T = 288; # time step per day
M = 1/12; # duaration per step in hour
ID = 1;
ED = 365;
N_sim = ED-ID+1; # number of days
L = RTP[:,1];
Lp = L .> 0

# BES setting
E = [0.2 0.2 0.2 0.2 0.2]
Pmax = 0.25; # power rating MW
S = length(E)
# Nonlinear Setting
Pd = transpose(Pmax*[0.7 1.0 1.0 0.9 0.5])
Pc = transpose(Pmax*[0.5 0.8 1.0 1.0 1.0])
eta = [0.80 0.85 0.90 0.85 0.80]
MC = [25 20 15 20 25]; # marginal discharge cost

# Linear Setting
# Pd = transpose(Pmax*[1.0 1.0 1.0 1.0 1.0])
# Pc = transpose(Pmax*[1.0 1.0 1.0 1.0 1.0])
# eta = [0.90 0.90 0.90 0.90 0.90]
# MC = [20 20 20 20 20]; # marginal discharge cost

e0 = [0.0 0.0 0.0 0.0 0.0]
ef = e0


# initialize optimization model
model = Model(Gurobi.Optimizer)
set_silent(model) # no outputs
# discharge power
@variable(model, d[1:T,1:S], lower_bound = 0)
# charge power
@variable(model, c[1:T,1:S], lower_bound = 0)
# energy level
@variable(model, e[1:T,1:S], lower_bound = 0)
@variable(model, C[1:S]) # value of battery capacity at the end of operation
@variable(model, R[1:S]) # market revenue
# binary variables
@variable(model, u[1:T,1:S], Bin)

# arbitrage revenue
@constraint(model, ArbRev[s=1:S], R[s] == M*sum(L.*(d[:,s]-c[:,s])) )
# piece-wise linear degradation (marginal discharge cost)
@constraint(model, DegCost[s=1:S], C[s] == M*MC[s]*sum(d[:,s]) )
# initial SoC evolution
@constraint(model, SoCInit[s=1:S], e[1,s] - e0[s] == M*(c[1,s]*eta[s] - d[1,s]/eta[s]) )
# rest SoC evolution
@constraint(model, SoCCont[t=2:T,s=1:S], e[t,s] - e[t-1,s] == M*(c[t,s]*eta[s] - d[t,s]/eta[s]) )
# final energy level
@constraint(model, Enelast[s=1:S], e[T,s] >= ef[s] )
# max discharge power
@constraint(model, DisMax[t=1:T], sum(d[t,:]./Pd) <= Lp[t] )
# max charge power
@constraint(model, ChrMax[t=1:T], sum(c[t,:]./Pc) <= 1 )
# max energy level
@constraint(model, SoCMax1[t=1:T], e[t,1] <= E[1] )
@constraint(model, SoCMax2[t=1:T,s=2:S], e[t,s] <= E[s]*u[t,s-1] )
# min energy level
@constraint(model, SoCMin1[t=1:T,s=1:S-1], e[t,s] >= E[s]*u[t,s] )
@constraint(model, SoCMin2[t=1:T], e[t,S] >= 0 )

# maximize revenue plus degradation value
@objective(model, Max, sum(R-C))



# initialize
R_s = zeros(1, N_sim)
P_s = zeros(1, N_sim)
C_s = zeros(1, N_sim)
S_s = zeros(5,N_sim*288)

@time begin

@printf("Optimization starts...\n")
for n = ID:ED

# update prices
local L = RTP[:,n]
# local L = collect(1:288)
local Lp = L .> 0

# update prices in constraints
for t = 1:T
    set_normalized_rhs(DisMax[t], Lp[t])
    for s = 1:S
        set_normalized_coefficient(ArbRev[s], d[t,s], -M*L[t] )
        set_normalized_coefficient(ArbRev[s], c[t,s], M*L[t] )
    end
end


optimize!(model)

global R_s[n] = value(sum(R))# objective_value(model);
global C_s[n] = value(sum(C))
global P_s[n] = value(sum(R-C))
global S_s[:,(n-1)*T+1:n*T] = transpose(value.(e))


termination_status(model)
@printf("Finished Day %d, Cum Rev %d, Cum Profit %d, Cum Cost %d, OptStatus: %s \n", n, sum(R_s), sum(P_s), sum(C_s), termination_status(model))


end
end
