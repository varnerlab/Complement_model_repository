# Estimates model parameters for proof-of-concept biochemical model -
using POETs
include("complement_lib.jl")

function estimate_complement_ensemble()

  number_of_subdivisions = 20
  number_of_parameters = 28
  number_of_objectives = 2

  ec_array_full = readdlm("./data/ec_array_O1_O2_O3.dat")
  total_error = sum(ec_array_full[:,1:end],1)
  min_index = indmin(total_error)
  pc_array_full = readdlm("./data/pc_array_O1_O2_O3.dat")
  initial_parameter_array = pc_array_full[:,min_index]
  #initial_parameter_array = 1.0*ones(28)
  #initial_parameter_array = vec(readdlm("best_solution.txt"))

  ec_array = zeros(number_of_objectives)
  pc_array = zeros(number_of_parameters)
  for index in collect(1:number_of_subdivisions)

    # Run JuPOETs -
    (EC,PC,RA) = estimate_ensemble(objective_function,neighbor_function,acceptance_probability_function,cooling_function,initial_parameter_array;rank_cutoff=4,maximum_number_of_iterations=20,show_trace=true)

    # Package -
    ec_array = [ec_array EC]
    pc_array = [pc_array PC]

    # Take the *best* value from the current ensemble, refine it and go around again -
    total_error = sum(ec_array[:,2:end],1)

    # Which col is the min error?
    min_index = indmin(total_error)
    @show index,total_error[min_index]

    # Refine the best solution -
    initial_parameter_array = pc_array[:,min_index].*(1+0.15*randn(number_of_parameters))
    #initial_parameter_array = local_refienment_step(initial_parameter_array)
  end

  return (ec_array,pc_array)
end
