# Script to sample the sensitivity of the C3a residual for the Alternate pathway -

# Run the model for the ensemble of solutions -
include("SolveBalances.jl")

# some global parameters -
BIG = 1e10
SMALL = 1e-6

using PyCall
@pyimport numpy as np

function calculate_error_obj(T,X,MEASURED_ARRAY_C3a,MEASURED_ARRAY_C5a)

  obj_array = BIG*ones(2)

  # Need to interpolate the simulation onto the experimental time scale -
  number_of_measurements = 20
  IC3a = np.interp(MEASURED_ARRAY_C3a[:,1],T,X[:,19])
  IC5a = np.interp(MEASURED_ARRAY_C5a[:,1],T,X[:,15])

  # Compute the error values -
  MC3a = 0.11002*MEASURED_ARRAY_C3a[:,2]
  MC5a = 0.00009615*MEASURED_ARRAY_C5a[:,2]
  # error_C3a = sum((MC3a - IC3a).^2)
  # error_C5a = sum((MC5a - IC5a).^2)

  # Compute the scaled IC3a and IC5a -
  min_IC3a = minimum(IC3a)
  max_IC3a = maximum(IC3a)
  min_IC5a = minimum(IC5a)
  max_IC5a = maximum(IC5a)

  min_MC3a = minimum(MC3a)
  max_MC3a = maximum(MC3a)
  min_MC5a = minimum(MC5a)
  max_MC5a = maximum(MC5a)

  delta_C3a = max_IC3a - min_IC3a
  delta_C5a = max_IC5a - min_IC5a

  # if (delta_C3a<SMALL)
  #   delta_C3a = SMALL
  # end
  #
  # if (delta_C5a<SMALL)
  #   delta_C5a = SMALL
  # end

  shape_C3a = (IC3a - min_IC3a)/(delta_C3a)
  shape_C5a = (IC5a - min_IC5a)/(delta_C5a)


  shape_MC3a = (MC3a - min_MC3a)/(max_MC3a - min_MC3a)
  shape_MC5a = (MC5a - min_MC5a)/(max_MC5a - min_MC5a)

  # Compute the scale error -
  scale_error_C3a = (max_IC3a - max_MC3a)/(max_MC3a)
  scale_error_C5a = (max_IC5a - max_MC5a)/(max_MC5a)

  #error_C3a = sum((shape_C3a - shape_MC3a).^2)+(scale_error_C3a)^2
  #error_C5a = sum((shape_C5a - shape_MC5a).^2)+(scale_error_C5a)^2
  error_C3a = sum((IC3a - MC3a).^2)
  error_C5a = sum((IC5a - MC5a).^2)

  # Compute the objective array -
  obj_array[1] = error_C3a
  obj_array[2] = error_C5a

  value = sum(obj_array);
  if (value>10)
    value = 10
  end

  return (value)
end


function main()

  # Load the data -
  MEASURED_ARRAY_C3a = readdlm("./data/Shaw2015_Fig2e_C3a.txt",',')
  MEASURED_ARRAY_C5a = readdlm("./data/Shaw2015_Fig3c_C5a.txt",',')

  # Load the paremeter samples -
  pc_array_full = readdlm("./sens/Samples.txt")
  (nrows,ncols) = size(pc_array_full)

  # for each sample, grab the parameters, solve the model and compute the error -
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

  tStart = 0.0;
  tStop = 25.0;
  tStep = 0.1;

  # initialize the error_array -
  error_array = zeros(nrows)

  for sample_index = 1:nrows

    # Set the parameters -
    parameter_array = pc_array_full[sample_index,:]
    data_dictionary["PARAMETER_ARRAY"] = parameter_array

    # Run the model -
    (T,X) = SolveBalances(tStart,tStop,tStep,data_dictionary)

    # Calculate the error -
    error_array[sample_index] = calculate_error_obj(T,X,MEASURED_ARRAY_C3a,MEASURED_ARRAY_C5a)

    @show sample_index
  end

  return error_array
end
