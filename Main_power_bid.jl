using JuMP
using Gurobi
using MAT
using Printf
using DataFrames
using CSV
# using Plots
using JLD
using DelimitedFiles


# High-level Settings
Zone = "WALNUT" # price zone name

# read price
fileln = matopen(string("./Data/",Zone,".mat"))
RTP = read(fileln, "Q")
close(fileln)
dbid = readdlm("./Data/WALNUT_2016_Discharge_One.csv", ',', Float64)
cbid = readdlm("./Data/WALNUT_2016_Charge_One.csv", ',', Float64)

# simulation setting
T = 288; # time step per day
M = 1/12; # duaration per step in hour
ID = 1;
ED = 365;
N_sim = ED-ID+1; # number of days
L = RTP[1,1];
Lp = L .> 0
DB = dbid[:,1]
CB = cbid[:,1]

# BES setting
E = [0.2 0.2 0.2 0.2 0.2]
Pmax = 0.25; # power rating MW
S = length(E)
d0 = zeros(S,1)
c0 = zeros(S,1)

# Linear Setting
Pd = transpose(Pmax*[1.0 1.0 1.0 1.0 1.0])
Pc = transpose(Pmax*[1.0 1.0 1.0 1.0 1.0])
eta = [0.90 0.90 0.90 0.90 0.90]
MC = [20 20 20 20 20]; # marginal discharge cost

# Nonlinear Setting
Pd2 = transpose(Pmax*[0.7 1.0 1.0 0.9 0.5])
Pc2 = transpose(Pmax*[0.5 0.8 1.0 1.0 1.0])
eta2 = [0.80 0.85 0.90 0.85 0.80]
MC2 = [25 20 15 20 25]; # marginal discharge cost

e0 = [0.0 0.0 0.0 0.0 0.0]
ef = e0


# initialize optimization model
model = Model(Gurobi.Optimizer)
set_silent(model) # no outputs
# discharge power
@variable(model, d[1:S], lower_bound = 0)
# charge power
@variable(model, c[1:S], lower_bound = 0)
# energy level
@variable(model, e[1:S], lower_bound = 0)
@variable(model, C[1:S]) # value of battery capacity at the end of operation
@variable(model, R[1:S]) # market revenue
@variable(model, V[1:S]) # opportunity value
# binary variables
@variable(model, u[1:S], Bin)

# arbitrage revenue
@constraint(model, ArbRev[s=1:S], R[s] == M*L*(d[s]-c[s]) )
# piece-wise linear degradation (marginal discharge cost)
@constraint(model, DegCost[s=1:S], C[s] == M*MC[s]*d[s] )
# piece-wise linear degradation opportunity value
@constraint(model, OppoVal[s=1:S], V[s] == M*(CB[s]*c[s]-DB[s]*d[s]) )
# SoC evolution
@constraint(model, SoCEvo[s=1:S], e[s] - e0[s] == M*(c[s]*eta[s] - d[s]/eta[s]) )
# max discharge power
@constraint(model, DisMax, sum(d[:]./Pd) <= Lp )
# max charge power
@constraint(model, ChrMax, sum(c[:]./Pc) <= 1 )
# max energy level
@constraint(model, SoCMax1, e[1] <= E[1] )
@constraint(model, SoCMax2[s=2:S], e[s] <= E[s]*u[s-1] )
# min energy level
@constraint(model, SoCMin1[s=1:S-1], e[s] >= E[s]*u[s] )
@constraint(model, SoCMin2, e[S] >= 0 )

# maximize revenue plus degradation value
@objective(model, Max, sum(R+V))


# initialize optimization model
model2 = Model(Gurobi.Optimizer)
set_silent(model2) # no outputs
# discharge power
@variable(model2, d2[1:S], lower_bound = 0)
# charge power
@variable(model2, c2[1:S], lower_bound = 0)
# energy level
@variable(model2, e2[1:S], lower_bound = 0)
@variable(model2, C2[1:S]) # value of battery capacity at the end of operation
@variable(model2, R2[1:S]) # market revenue
# binary variables
@variable(model2, u2[1:S], Bin)

# arbitrage revenue
@constraint(model2, ArbRev2[s=1:S], R2[s] == M*L*(d2[s]-c2[s]) )
# piece-wise linear degradation (marginal discharge cost)
@constraint(model2, DegCost2[s=1:S], C2[s] == M*MC2[s]*d2[s] )
# SoC evolution
@constraint(model2, SoCEvo2[s=1:S], e2[s] - e0[s] == M*(c2[s]*eta2[s] - d2[s]/eta2[s]) )
# max discharge power
@constraint(model2, DisMax2, sum(d2[:]./Pd2) <= Lp )
# max charge power
@constraint(model2, ChrMax2, sum(c2[:]./Pc2) <= 1 )
# max energy level
@constraint(model2, SoCMax12, e2[1] <= E[1] )
@constraint(model2, SoCMax22[s=2:S], e2[s] <= E[s]*u2[s-1] )
# min energy level
@constraint(model2, SoCMin12[s=1:S-1], e2[s] >= E[s]*u2[s] )
@constraint(model2, SoCMin22, e2[S] >= 0 )

# maximize revenue plus degradation value
@objective(model2, Min, sum((d2-d0).^2+(c2-c0).^2))


# initialize
R_s = zeros(T, N_sim)
P_s = zeros(T, N_sim)
C_s = zeros(T, N_sim)
R_s2 = zeros(T, N_sim)
P_s2 = zeros(T, N_sim)
C_s2 = zeros(T, N_sim)

@time begin

@printf("Optimization starts...\n")
for n = ID:ED
for t = 1:T
    H = floor(Int,((n-1)*T+t-1)*M)+1
    # update prices
    local L = RTP[t,n]
    # local L = collect(1:288)
    local Lp = L .> 0
    #updata bids
    local DB = dbid[:,H]
    local CB = cbid[:,H]
    set_normalized_rhs(DisMax, Lp)
    for s = 1:S
        set_normalized_coefficient(ArbRev[s], d[s], -M*L )
        set_normalized_coefficient(ArbRev[s], c[s], M*L )
        set_normalized_coefficient(OppoVal[s], d[s], M*DB[s] )
        set_normalized_coefficient(OppoVal[s], c[s], -M*CB[s] )
        set_normalized_rhs(SoCEvo[s], e0[s])
    end
    optimize!(model)

    local d0 = value.(d)
    local c0 = value.(c)
    new_objective = @expression(model2, sum((d2-d0).^2+(c2-c0).^2))
    set_objective_function(model2, new_objective)
    set_normalized_rhs(DisMax2, Lp)
    for s = 1:S
        set_normalized_coefficient(ArbRev2[s], d2[s], -M*L )
        set_normalized_coefficient(ArbRev2[s], c2[s], M*L )
        set_normalized_rhs(SoCEvo2[s], e0[s])
    end
    optimize!(model2)

    #update SoC
    for s = 1:S
        global e0[s] = round(value(e2[s]),digits=6)
    end

    global R_s[t,n] = value(sum(R))# objective_value(model);
    global C_s[t,n] = value(sum(C))
    global P_s[t,n] = value(sum(R-C))
    termination_status(model)
    global R_s2[t,n] = value(sum(R2))# objective_value(model);
    global C_s2[t,n] = value(sum(C2))
    global P_s2[t,n] = value(sum(R2-C2))
    termination_status(model2)
end

@printf("Finished Day %d, Cum Rev %d, Cum Profit %d, Cum Cost %d, OptStatus: %s \n", n, sum(R_s), sum(P_s), sum(C_s), termination_status(model))
@printf("Finished Day %d, Cum Rev %d, Cum Profit %d, Cum Cost %d, OptStatus: %s \n", n, sum(R_s2), sum(P_s2), sum(C_s2), termination_status(model2))

end
end
