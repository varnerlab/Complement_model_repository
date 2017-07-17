# Run the model for the ensemble of solutions -
include("SolveBalances.jl")
include("Common.jl")
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
distance_array_set = Set{Array{Float64,2}}()

data_dictionary = Dict()
initial_condition_array = zeros(19)
initial_condition_array[1] = 0.0    # 1 zymosan
initial_condition_array[2] = 1.9      # 2 C4
initial_condition_array[3] = 0.322    # 3 C2
initial_condition_array[8] = 7.57     # 4 C3
initial_condition_array[14] = 0.195   # 5 C5
initial_condition_array[17] = 2.23    # 6 Factor H
initial_condition_array[18] = 0.417   # 7 C4BP

initial_condition_array[19] = 0.55     # 19 C3a (ic initial)
initial_condition_array[15] = 0.0006   # 15 C5a (ic initial)

data_dictionary["INITIAL_CONDITION_ARRAY"] = initial_condition_array;
number_parameters = 0

for (index,rank_index_value) in enumerate(idx_rank)

  if (mod(index,10)==0)

    # Grab the parameter set from the cache -
    parameter_array = pc_array_full[:,rank_index_value]
    data_dictionary["PARAMETER_ARRAY"] = parameter_array

    number_parameters = length(parameter_array)

    # Run the model -
    (t_base,x_base) = SolveBalances(tStart,tStop,tStep,data_dictionary)

    # Initialize the distance array -
    distance_array = zeros(length(parameter_array),length(parameter_array))

    # For each parameter combination, we need to measured the distance between the systems states -
    for (parameter_index_outer,parameter_value_outer) in enumerate(parameter_array)
      for (parameter_index_inner,parameter_value_inner) in enumerate(parameter_array)

        # Update the parameters -
        copy_of_parameter_array = copy(parameter_array)
        copy_of_parameter_array[parameter_index_outer] = parameter_value_outer*(1.10)
        copy_of_parameter_array[parameter_index_inner] = parameter_value_inner*(1.10)

        # Update the data dictionary -
        data_dictionary["PARAMETER_ARRAY"] = copy_of_parameter_array

        # Run the model -
        (t_perturbed,x_pertubed) = SolveBalances(tStart,tStop,tStep,data_dictionary)

        # measure the distance -
        local_distance = norm(norm(x_base - x_pertubed))

        # add this value to the distance_array -
        distance_array[parameter_index_outer,parameter_index_inner] = local_distance

        @show (parameter_index_outer,parameter_index_inner,index)
      end # end inner -
    end # end outer -

    # When I get here, I should have a full distance array -
    push!(distance_array_set,distance_array)

  end # end the mod -
end # end rank iteration -

# write -
write_array_set_to_disk(copy(distance_array_set),"./clustergram_data","DArray-control")

# calculate the array -
mean_array = calculate_mean_array_from_set(copy(distance_array_set),number_parameters,number_parameters)
