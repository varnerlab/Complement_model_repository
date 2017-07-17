# Run the model for the ensemble of solutions -
include("SolveBalances.jl")
using PyPlot
using PyCall
@pyimport numpy as np


# Load the ensemble files from disk -
pc_array_full = readdlm("./data/pc_array_O1_O2.dat")
rank_array = readdlm("./data/rank_array_O1_O2.dat")

# Select the desired rank -
idx_rank = find(rank_array .<= 5.0)

# Setup time scale -
tStart = 0.0;
tStop = 25.0;
tStep = 0.1;
number_of_timesteps = 200
time_experimental = linspace(tStart,tStop,number_of_timesteps)

# initialize data array -
data_array = zeros(length(time_experimental),1)

data_dictionary = Dict()
initial_condition_array = zeros(19)
initial_condition_array[1] = 0.0      # 1 zymosan
initial_condition_array[2] = 1.9      # 2 C4
initial_condition_array[3] = 0.322    # 3 C2
initial_condition_array[8] = 7.57     # 4 C3
initial_condition_array[14] = 0.195   # 5 C5
initial_condition_array[17] = 2.23    # 6 Factor H
initial_condition_array[18] = 0.417   # 7 C4BP

initial_condition_array[19] = 0.6    # 19 C3a (ic initial)
initial_condition_array[15] = 0.00006   # 19 C5a (ic initial)

data_dictionary["INITIAL_CONDITION_ARRAY"] = initial_condition_array;
bad_set = Set()

# Run the simulation for this rank -
number_of_samples = 1
for (index,rank_index_value) in enumerate(idx_rank)

  if (mod(index,1)==0)

    # Grab the parameter set from the cache -
    parameter_array = pc_array_full[:,rank_index_value]
    data_dictionary["PARAMETER_ARRAY"] = parameter_array

    # Run the model -
    (t,x) = SolveBalances(tStart,tStop,tStep,data_dictionary)

    # @show x[1,1]

    # Need to interpolate the simulation onto the experimental time scale -
    IC3a = np.interp(time_experimental[:],t,x[:,19])
    IC5a = np.interp(time_experimental[:],t,x[:,15])

    # grab -
    data_array = [data_array IC3a]

    # How many samples?
    number_of_samples = number_of_samples + 1

    @show index
  end
end

# calculate mean, and std
mean_value = mean(data_array,2)
std_value = std(data_array,2)


# make best fit -
SF = 2.58
UB = mean_value + (SF)*std_value
LB = mean_value - (SF)*std_value
idx_z = find(LB.<0)
LB[idx_z] = 0.0
fill_between(vec(time_experimental),vec(LB),vec(UB),color=[0.8,0.8,0.8],lw=2)

# 99% of mean -
SF = (2.58/sqrt(number_of_samples))
UB = mean_value + (SF)*std_value
LB = mean_value - (SF)*std_value
idx_z = find(LB.<0)
LB[idx_z] = 0.0

# Make the plot -
plot(time_experimental,mean_value,"k-")
fill_between(vec(time_experimental),vec(LB),vec(UB),color="gray",lw=2)

# Load the data -
# 0.0 mg/ml -
MEASURED_ARRAY_C3a = readdlm("./data/Shaw2015_Fig2a_C3a.txt",',')
MEASURED_ARRAY_C5a = readdlm("./data/Shaw2015_Fig3ai_C5a_original.txt",',')

MC3a = 0.11002*MEASURED_ARRAY_C3a[:,2]
MC5a = 0.00009615*MEASURED_ARRAY_C5a[:,2]
plot(MEASURED_ARRAY_C3a[:,1],MC3a,color="black",marker="o",markersize=8.0)
