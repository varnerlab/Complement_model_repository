# include model data -
include("SolveBalances.jl")

using PyCall
@pyimport numpy as np

using Debug

# some global parameters -
BIG = 1e10
SMALL = 1e-6

# Globally load data so we don't have load on each iteration -
MEASURED_ARRAY_C3a_O1 = readdlm("./data/Shaw2015_Fig2a_C3a.txt",',')
MEASURED_ARRAY_C5a_O1 = readdlm("./data/Shaw2015_Fig3ai_C5a.txt",',')

MEASURED_ARRAY_C3a_O2 = readdlm("./data/Shaw2015_Fig2e_C3a.txt",',')
MEASURED_ARRAY_C5a_O2 = readdlm("./data/Shaw2015_Fig3c_C5a_original.txt",',')

MEASURED_ARRAY_C3a_O3 = readdlm("./data/Shaw2015_Fig2c_C3a.txt",',')
MEASURED_ARRAY_C5a_O3 = readdlm("./data/Shaw2015_Fig3aiii_C5a_original.txt",',')

function calculate_error_obj_1(T,X)

  obj_array = BIG*ones(2)

  # Need to interpolate the simulation onto the experimental time scale -
  number_of_measurements = 20
  IC3a = np.interp(MEASURED_ARRAY_C3a_O1[:,1],T,X[:,19])
  IC5a = np.interp(MEASURED_ARRAY_C5a_O1[:,1],T,X[:,15])

  # Compute the error values -
  MC3a = 0.11002*MEASURED_ARRAY_C3a_O1[:,2]
  MC5a = 0.00009615*MEASURED_ARRAY_C5a_O1[:,2]

  # Compute the scaled IC3a and IC5a -
  min_IC3a = minimum(IC3a)
  max_IC3a = maximum(IC3a)
  min_IC5a = minimum(IC5a)
  max_IC5a = maximum(IC5a)

  min_MC3a = minimum(MC3a)
  max_MC3a = maximum(MC3a)
  min_MC5a = minimum(MC5a)
  max_MC5a = maximum(MC5a)

  shape_C3a = (IC3a - min_IC3a)/(max_IC3a - min_IC3a)
  shape_C5a = (IC5a - min_IC5a)/(max_IC5a - min_IC5a)
  shape_MC3a = (MC3a - min_MC3a)/(max_MC3a - min_MC3a)
  shape_MC5a = (MC5a - min_MC5a)/(max_MC5a - min_MC5a)

  # Compute the scale error -
  scale_error_C3a = (max_IC3a - max_MC3a)/(max_MC3a)
  scale_error_C5a = (max_IC5a - max_MC5a)/(max_MC5a)

  error_C3a = sum((shape_C3a - shape_MC3a).^2) + 10*(scale_error_C3a)^2
  error_C5a = sum((shape_C5a - shape_MC5a).^2) + 10*(scale_error_C5a)^2

  # Compute the objective array -
  obj_array[1] = error_C3a
  obj_array[2] = error_C5a

  return (sum(obj_array))
end

function calculate_error_obj_2(T,X)

  obj_array = BIG*ones(2)

  # Need to interpolate the simulation onto the experimental time scale -
  number_of_measurements = 20
  IC3a = np.interp(MEASURED_ARRAY_C3a_O2[:,1],T,X[:,19])
  IC5a = np.interp(MEASURED_ARRAY_C5a_O2[:,1],T,X[:,15])

  # Compute the error values -
  MC3a = 0.11002*MEASURED_ARRAY_C3a_O2[:,2]
  MC5a = 0.00009615*MEASURED_ARRAY_C5a_O2[:,2]

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

  shape_C3a = (IC3a - min_IC3a)/(max_IC3a - min_IC3a)
  shape_C5a = (IC5a - min_IC5a)/(max_IC5a - min_IC5a)
  shape_MC3a = (MC3a - min_MC3a)/(max_MC3a - min_MC3a)
  shape_MC5a = (MC5a - min_MC5a)/(max_MC5a - min_MC5a)

  # Compute the scale error -
  scale_error_C3a = (max_IC3a - max_MC3a)/(max_MC3a)
  scale_error_C5a = (max_IC5a - max_MC5a)/(max_MC5a)

  error_C3a = sum((shape_C3a - shape_MC3a).^2) + 10*(scale_error_C3a)^2
  error_C5a = sum((shape_C5a - shape_MC5a).^2) + 10*(scale_error_C5a)^2

  # Compute the objective array -
  obj_array[1] = error_C3a
  obj_array[2] = error_C5a

  return (sum(obj_array))
end

function calculate_error_obj_3(T,X)

  obj_array = BIG*ones(2)

  # Need to interpolate the simulation onto the experimental time scale -
  number_of_measurements = 20
  IC3a = np.interp(MEASURED_ARRAY_C3a_O3[:,1],T,X[:,19])
  IC5a = np.interp(MEASURED_ARRAY_C5a_O3[:,1],T,X[:,15])

  # Compute the error values -
  MC3a = 0.11002*MEASURED_ARRAY_C3a_O3[:,2]
  MC5a = 0.00009615*MEASURED_ARRAY_C5a_O3[:,2]
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

  shape_C3a = (IC3a - min_IC3a)/(max_IC3a - min_IC3a)
  shape_C5a = (IC5a - min_IC5a)/(max_IC5a - min_IC5a)
  shape_MC3a = (MC3a - min_MC3a)/(max_MC3a - min_MC3a)
  shape_MC5a = (MC5a - min_MC5a)/(max_MC5a - min_MC5a)

  # Compute the scale error -
  scale_error_C3a = (max_IC3a - max_MC3a)/(max_MC3a)
  scale_error_C5a = (max_IC5a - max_MC5a)/(max_MC5a)

  error_C3a = sum((shape_C3a - shape_MC3a).^2) + (scale_error_C3a)^2
  error_C5a = sum((shape_C5a - shape_MC5a).^2) + (scale_error_C5a)^2

  # Compute the objective array -
  obj_array[1] = error_C3a
  obj_array[2] = error_C5a

  return (sum(obj_array))
end

function objective_function(parameter_array)

  # Setup the timescale -
  tStart = 0.0;
  tStop = 25.0;
  tStep = 0.1;

  # Calculate the objective function array -
  obj_array = BIG*ones(2,1)

  # Setup data dictionary -
  data_dictionary = Dict()
  data_dictionary["PARAMETER_ARRAY"] = parameter_array

  # Experiment 1 = w/o zymosan (run with default ICs)
  initial_condition_array = zeros(19)
  initial_condition_array[1] = 0.0      # 1 zymosan
  initial_condition_array[2] = 1.9      # 2 C4
  initial_condition_array[3] = 0.322    # 3 C2
  initial_condition_array[8] = 7.57     # 4 C3
  initial_condition_array[14] = 0.195   # 5 C5
  initial_condition_array[17] = 2.23    # 6 Factor H
  initial_condition_array[18] = 0.417   # 7 C4BP
  data_dictionary["INITIAL_CONDITION_ARRAY"] = initial_condition_array;
  (T1,X1) = SolveBalances(tStart,tStop,tStep,data_dictionary)

  # Experiment 2 = w 1g zymosan -
  initial_condition_array[1] = 1.0      # 1 zymosan
  data_dictionary["INITIAL_CONDITION_ARRAY"] = initial_condition_array;
  (T2,X2) = SolveBalances(tStart,tStop,tStep,data_dictionary)

  # # Experiment 3 = w 0.01g zymosan -
  # initial_condition_array[1] = 0.01      # 1 zymosan
  # data_dictionary["INITIAL_CONDITION_ARRAY"] = initial_condition_array;
  # (T3,X3) = SolveBalances(tStart,tStop,tStep,data_dictionary)

  # call the experiment functions -
  obj_array[1] = calculate_error_obj_1(T1,X1)
  obj_array[2] = calculate_error_obj_2(T2,X2)
  #obj_array[3] = calculate_error_obj_3(T3,X3)

  # return -
  return obj_array
end


# --------------------------------------- DO NOT EDIT BELOW THIS LINE  --------------------------------------- #
# Evaluates the objective function values -
function local_refienment_step(parameter_array)

  SIGMA = 0.05

  # initialize -
  number_of_parameters = length(parameter_array)

  # Setup the bound constraints -
  LOWER_BOUND = SMALL
  UPPER_BOUND = 10*ones(number_of_parameters)

  # calculate the starting error -
  parameter_array_best = parameter_array
  error_array = BIG*ones(4)
  error_array[1] = sum(objective_function(parameter_array_best))

  # main refinement loop -
  iteration_counter = 1
  iteration_max = 1000
  while (iteration_counter<iteration_max)

    # take a step up -
    parameter_up = parameter_array_best.*(1+SIGMA*rand(number_of_parameters))
    parameter_up = parameter_bounds_function(parameter_up,LOWER_BOUND*ones(number_of_parameters),UPPER_BOUND)

    # take a step down -
    parameter_down = parameter_array_best.*(1-SIGMA*rand(number_of_parameters))
    parameter_down = parameter_bounds_function(parameter_down,LOWER_BOUND*ones(number_of_parameters),UPPER_BOUND)

    # Evaluate the obj function -
    error_array[2] = sum(objective_function(parameter_up))
    error_array[3] = sum(objective_function(parameter_down))

    # Calculate a correction factor -
    a = error_array[2]+error_array[3] - 2.0*error_array[1]
    parameter_corrected = parameter_array_best
    if (a>0.0)
      amda = -0.5*(error_array[3] - error_array[2])/a
      parameter_corrected = parameter_array_best+amda*rand(number_of_parameters)
      parameter_corrected = parameter_bounds_function(parameter_corrected,LOWER_BOUND*ones(number_of_parameters),UPPER_BOUND)
      error_array[4] = sum(objective_function(parameter_corrected))
    end

    # Which step has the min error?
    min_index = indmin(error_array)
    if (min_index == 1)
			parameter_array_best = parameter_array_best
	 	elseif (min_index == 2)
		 	parameter_array_best = parameter_up
	 	elseif (min_index == 3)
		 	parameter_array_best = parameter_down
	 	elseif (min_index == 4)
		 parameter_array_best = parameter_corrected
	 	end

		# Update the local error
		error_array[1] = error_array[min_index]

    @show iteration_counter,error_array[min_index]

    # update local counter -
    iteration_counter = iteration_counter + 1
  end

  return parameter_array_best
end

function neighbor_function(parameter_array)

  SIGMA = 0.05
  number_of_parameters = length(parameter_array)

  # calculate new parameters -
  new_parameter_array = parameter_array.*(1+SIGMA*randn(number_of_parameters))

  # Check the bound constraints -
  LOWER_BOUND = SMALL
  UPPER_BOUND = [
        50.0  ;          # k_Initiator2_C4             = 0.1*parameter_array[1];
        15.0  ;          # Km_Initiator2_C4            = 1.0*parameter_array[2];
        50.0  ;          # k_Initiator2_C2             = 0.1*parameter_array[3];
        15.0  ;          # Km_Initiator2_C2            = 1.0*parameter_array[4];
        2.0   ;          # n_In2_C4                    = parameter_array[5];
        2.0   ;          # n_In2_C2                    = parameter_array[6];
        15.0  ;          # alpha_In2_C4                = parameter_array[7];
        15.0  ;          # alpha_In2_C2                = parameter_array[8];
        2.0   ;          # n_C5a                       = parameter_array[9];
        2.0   ;          # n_C3a                       = parameter_array[10];
        15.0  ;          # k_C4bC2a                    = parameter_array[11];
        15.0  ;          # k_C3b_basal                 = parameter_array[12];
        15.0  ;          # k_C3Convertase2             = parameter_array[13];
        15.0  ;          # k_C5Convertase2             = parameter_array[14];
        15.0  ;          # k_C5_C4bC2a_bind            = parameter_array[15];
        15.0  ;          # k_C5_conversion             = parameter_array[16];
        15.0  ;          # Km_C5_conversion            = parameter_array[17];
        15.0  ;          # k_C5_conversion_alternate   = parameter_array[18];
        15.0  ;          # Km_C5_conversion_alternate  = parameter_array[19];
        15.0  ;          # k_inhibit_C4bC2a            = parameter_array[20];
        15.0  ;          # k_inhibit_C3Convertase2     = parameter_array[21];
        15.0  ;          # k_cat_C4bC2a                = parameter_array[22];
        15.0  ;          # Km_C4bC2a                   = parameter_array[23];
        15.0  ;          # k_cat_C3convertase2         = parameter_array[24];
        15.0  ;          # Km_C3Convertase2            = parameter_array[25];
        1.0   ;          # k_degradationC3a            = parameter_array[26];
        15.0  ;          # kC5inhibit                  = parameter_array[27];
        1.0   ;          # k_degradationC5a            = parameter_array[28];
  ];

  # return the corrected parameter arrays -
  return parameter_bounds_function(new_parameter_array,LOWER_BOUND*ones(number_of_parameters),UPPER_BOUND)
end


function acceptance_probability_function(rank_array,temperature)
  return (exp(-rank_array[end]/temperature))
end

function cooling_function(temperature)

  # define my new temperature -
  alpha = 0.90
  return alpha*temperature
end

function parameter_bounds_function(parameter_array,lower_bound_array,upper_bound_array)

  # reflection_factor -
  epsilon = 0.01

  # iterate through and fix the parameters -
  new_parameter_array = copy(parameter_array)
  for (index,value) in enumerate(parameter_array)

    lower_bound = lower_bound_array[index]
    upper_bound = upper_bound_array[index]

    if (value<lower_bound)
      new_parameter_array[index] = lower_bound
    elseif (value>upper_bound)
      new_parameter_array[index] = upper_bound
    end
  end

  return new_parameter_array
end
# --------------------------------------- DO NOT EDIT ABOVE THIS LINE  --------------------------------------- #
