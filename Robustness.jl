# Julia script to compute robsutness coefficients -
include("SolveBalances.jl")
using PyPlot
using PyCall
@pyimport numpy as np

# main: calculates the ratio of the AUC in the perturbed state, to the base state
function main(readout_index,perturbed_index_array,delta_array)

  # Load the ensemble files from disk -
  pc_array_full = readdlm("./data/pc_array_O1_O2.dat")
  rank_array = readdlm("./data/rank_array_O1_O2.dat")

  # Select the desired rank -
  idx_rank = find(rank_array .<= 5.0)

  # initialized the robsutness array -
  alpha_array = Float64[]

  # Run the simulation for this rank -
  number_of_samples = 1
  for (index,rank_index_value) in enumerate(idx_rank)

    if (mod(index,1)==0)

      # Grab the parameter set from the cache -
      parameter_array = pc_array_full[:,rank_index_value]

      # Calculate the AUC for the base state -
      AUC_BASE = calculate_auc_base(parameter_array,readout_index)

      # Calculate the AUC for the perturbed case -
      AUC_PERTURBED = calculate_auc_perturbed(parameter_array,readout_index,perturbed_index_array,delta_array)

      @show AUC_BASE,AUC_PERTURBED

      alpha_value = log10(abs(AUC_PERTURBED/AUC_BASE))
      push!(alpha_array,alpha_value)

    end
  end

  return alpha_array
end

function calculate_auc_base(parameter_array,readout_index)

  tStart = 0.0;
  tStop = 25.0;
  tStep = 0.1;

  data_dictionary = Dict()
  initial_condition_array = zeros(19)
  initial_condition_array[1] = 1.0      # 1 zymosan
  initial_condition_array[2] = 1.9      # 2 C4
  initial_condition_array[3] = 0.322    # 3 C2
  initial_condition_array[8] = 7.57     # 4 C3
  initial_condition_array[14] = 0.195   # 5 C5
  initial_condition_array[17] = 2.23    # 6 Factor H
  initial_condition_array[18] = 0.417   # 7 C4BP
  data_dictionary["INITIAL_CONDITION_ARRAY"] = initial_condition_array;
  data_dictionary["PARAMETER_ARRAY"] = parameter_array

  # Run the model -
  (t,x) = SolveBalances(tStart,tStop,tStep,data_dictionary)

  # Calculate the AUC for readout_index -
  return trapz(t,x[:,readout_index])
end

function calculate_auc_perturbed(parameter_array,readout_index,perturbed_index_array,delta_array)

  tStart = 0.0;
  tStop = 25.0;
  tStep = 0.1;

  data_dictionary = Dict()
  initial_condition_array = zeros(19)
  initial_condition_array[1] = 1.0    # 1 zymosan
  initial_condition_array[2] = 1.9      # 2 C4
  initial_condition_array[3] = 0.322    # 3 C2
  initial_condition_array[8] = 7.57     # 4 C3
  initial_condition_array[14] = 0.195   # 5 C5
  initial_condition_array[17] = 2.23    # 6 Factor H
  initial_condition_array[18] = 0.417   # 7 C4BP

  # Issue the perturbation -
  for (index,perturbed_index) in enumerate(perturbed_index_array)
    delta = delta_array[index]
    initial_condition_array[perturbed_index] = delta*initial_condition_array[perturbed_index]
  end

  data_dictionary["INITIAL_CONDITION_ARRAY"] = initial_condition_array;
  data_dictionary["PARAMETER_ARRAY"] = parameter_array

  # Run the model -
  (t,x) = SolveBalances(tStart,tStop,tStep,data_dictionary)

  # Calculate the AUC for readout_index -
  return trapz(t,x[:,readout_index])
end

function trapz{Tx<:Number, Ty<:Number}(x::Vector{Tx}, y::Vector{Ty})
    # Trapezoidal integration rule
    local n = length(x)
    if (length(y) != n)
        error("Vectors 'x', 'y' must be of same length")
    end
    r = zero(zero(Tx) + zero(Ty))
    if n == 1; return r; end
    for i in 2:n
        r += (x[i] - x[i-1]) * (y[i] + y[i-1])
    end
    #= correction -h^2/12 * (f'(b) - f'(a))
    ha = x[2] - x[1]
    he = x[end] - x[end-1]
    ra = (y[2] - y[1]) / ha
    re = (y[end] - y[end-1]) / he
    r/2 - ha*he/12 * (re - ra)
    =#
    return r/2
end
