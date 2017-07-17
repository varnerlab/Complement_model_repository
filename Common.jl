function write_array_set_to_disk(data_set::Set{Array{Float64,2}},save_directory::String,file_prefix::String)

  counter = 1
  while (length(data_set)>0)

    # Grab a data array from the set -
    data_array = pop!(data_set)

    # write the array -
    full_path_string = save_directory*"/"*file_prefix*"_T$(counter).dat"
    writedlm(full_path_string,data_array)

    # update -
    counter = counter + 1
  end
end

function calculate_mean_array_from_set(data_set::Set{Array{Float64,2}},number_of_rows::Int,number_of_cols::Int)

    # Initialize the mean -
    mean_array = zeros(number_of_rows,number_of_cols)

    counter = 0
    while (length(data_set)>0)

      # Grab a data array from the set -
      data_array = pop!(data_set)

      # iterate -
      for row_index = 1:number_of_rows
        for col_index = 1:number_of_cols

          # calculate the mean value -
          mean_array[row_index,col_index] = (data_array[row_index,col_index]+counter*mean_array[row_index,col_index])/(counter+1)

        end
      end


      # update the counter -
      counter = counter + 1
    end

    # return -
    return mean_array
end

function rand_normal(mean, stdev)

    if stdev <= 0.0
        error("standard deviation must be positive")
    end

    u1 = rand()
    u2 = rand()
    r = sqrt(-2.0*log(u1))

    theta = 2.0*pi*u2
    return mean + stdev*r*sin(theta)
end
