# script to compute the gray-scale color for the robsutness figure

# Global bounds -
BIG = 2.5
SMALL = 1e-6

# load the robsutness data -
data_array = readdlm("./robustness/alpha_Z1_99.dat")


# transform -
(number_of_samples,number_of_cases) = size(data_array)
gray_scale_array = zeros(number_of_samples,number_of_cases)
for sample_index = 1:number_of_samples

  for case_index = 1:number_of_cases

    value = abs(data_array[sample_index,case_index])
    if (value>BIG)
      value = BIG
    end

    gray_scale_array[sample_index,case_index] = round(255*(1 - (value - SMALL)/(BIG - SMALL)))
  end
end
