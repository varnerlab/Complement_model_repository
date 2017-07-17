# Run the model for the ensemble of solutions -
include("SolveBalances.jl")
include("Common.jl")
using PyPlot
using PyCall
@pyimport numpy as np

# Load the data -
# 1 mg/ml
# MEASURED_ARRAY_C3a = readdlm("./data/Shaw2015_Fig2e_C3a.txt",',')
# MEASURED_ARRAY_C5a = readdlm("./data/Shaw2015_Fig3c_C5a.txt",',')

# Control -
# MEASURED_ARRAY_C3a = readdlm("./data/Shaw2015_Fig2a_C3a.txt",',')
# MEASURED_ARRAY_C5a = readdlm("./data/Shaw2015_Fig3ai_C5a_original.txt",',')

# 0.1 mg/ml -
# MEASURED_ARRAY_C3a = readdlm("./data/Shaw2015_Fig2d_C3a.txt",',')
# MEASURED_ARRAY_C5a = readdlm("./data/Shaw2015_Fig3b_C5a_original.txt",',')

# 0.01 mg/ml -
# MEASURED_ARRAY_C3a = readdlm("./data/Shaw2015_Fig2c_C3a.txt",',')
# MEASURED_ARRAY_C5a = readdlm("./data/Shaw2015_Fig3aiii_C5a_original.txt",',')

# 0.001 mg/ml -
MEASURED_ARRAY_C3a = readdlm("./data/Shaw2015_Fig2b_C3a.txt",',')
MEASURED_ARRAY_C5a = readdlm("./data/Shaw2015_Fig3aii_C5a_original.txt",',')

# Load the ensemble files from disk -
pc_array_full = readdlm("./data/pc_array_O1_O2.dat")
ec_array_full = readdlm("./data/ec_array_O1_O2.dat")

# Find lowest error solution -
total_error = sum(ec_array_full[:,1:end],1)

# Which col is the min error?
max_index = indmax(total_error)

# min error parameter set -
max_error_parameter_set = pc_array_full[:,max_index]
number_of_parameters = length(max_error_parameter_set)

# Generate N = 100 new *random* parameter family -
number_of_parameter_sets = 100
coefficient_of_variation = 0.29
random_parameter_family = zeros(number_of_parameters,number_of_parameter_sets)
for parameter_set_index = 1:number_of_parameter_sets
  for parameter_index = 1:number_of_parameters
    mean_value = max_error_parameter_set[parameter_index]
    random_parameter_family[parameter_index,parameter_set_index] = rand_normal(mean_value,coefficient_of_variation*mean_value)
  end
end

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
initial_condition_array[1] = 1.0      # 1 zymosan
initial_condition_array[2] = 1.9      # 2 C4
initial_condition_array[3] = 0.322    # 3 C2
initial_condition_array[8] = 7.57     # 4 C3
initial_condition_array[14] = 0.195   # 5 C5
initial_condition_array[17] = 2.23    # 6 Factor H
initial_condition_array[18] = 0.417   # 7 C4BP

initial_condition_array[19] = 0.55     # 19 C3a (ic initial)
initial_condition_array[15] = 0.006   # 19 C3a (ic initial)
data_dictionary["INITIAL_CONDITION_ARRAY"] = initial_condition_array;

# Setup the data set -
aic_data_set = Set{Array{Float64,1}}()

# Run the simulation for this rank -
number_of_samples = 1
for index = 1:number_of_parameter_sets

  # Grab the parameter set from the cache -
  parameter_array = random_parameter_family[:,index]
  data_dictionary["PARAMETER_ARRAY"] = parameter_array

  # How many parameters?
  number_of_parameters = length(parameter_array)

  # Run the model -
  (T,X) = SolveBalances(tStart,tStop,tStep,data_dictionary)

  # Calculate the error -
  IC3a = np.interp(MEASURED_ARRAY_C3a[:,1],T,X[:,19])
  IC5a = np.interp(MEASURED_ARRAY_C5a[:,1],T,X[:,15])

  # Compute the error values -
  MC3a = 0.11002*MEASURED_ARRAY_C3a[:,2]
  MC5a = 0.00009615*MEASURED_ARRAY_C5a[:,2]

  # Error -
  error_c3a = (1/norm(MC3a))*sum((IC3a - MC3a).^2)
  error_c5a = (1/norm(MC5a))*sum((IC5a - MC5a).^2)

  # How many steps?
  number_of_steps = length(MC3a)

  # Compute the aic -
  error_array = zeros(2)
  error_array[1] = 2*number_of_parameters+number_of_steps*log(error_c3a)
  error_array[2] = 2*number_of_parameters+number_of_steps*log(error_c5a)

  # Add -
  push!(aic_data_set,error_array)
  @show (index,2*number_of_parameters+number_of_steps*log(error_c3a),2*number_of_parameters+number_of_steps*log(error_c5a))
end

# ok, so lets comput the mean and std -
mean_array_C3 = Float64[]
mean_array_C5 = Float64[]
while (length(aic_data_set)>0)

  data_array = pop!(aic_data_set)
  push!(mean_array_C3,data_array[1])
  push!(mean_array_C5,data_array[2])
end

mean_value_C3 = mean(mean_array_C3)
mean_value_C5 = mean(mean_array_C5)
std_value_C3 = std(mean_array_C3)
std_value_C5 = std(mean_array_C5)

results_array = [mean_value_C3 std_value_C3 ; mean_value_C5 std_value_C5]
