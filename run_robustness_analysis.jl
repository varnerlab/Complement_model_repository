include("Robustness.jl")

index_C3a = 19
index_C5a = 15
index_C5 = 14
index_C3 = 8

alpha_array_final_50 = zeros(65)
alpha_array_final_90 = zeros(65)
alpha_array_final_99 = zeros(65)

# 50% reduction C3a and C5a -
alpha_array_C3a_C3_50 = main(index_C3a,[index_C3],[0.5])
alpha_array_C3a_C5_50 = main(index_C3a,[index_C5],[0.5])

alpha_array_C5a_C3_50 = main(index_C5a,[index_C3],[0.5])
alpha_array_C5a_C5_50 = main(index_C5a,[index_C5],[0.5])

alpha_array_C3a_C3_C5_50 = main(index_C3a,[index_C3,index_C5],[0.5,0.5])
alpha_array_C5a_C3_C5_50 = main(index_C5a,[index_C3,index_C5],[0.5,0.5])

alpha_array_final_50 = [alpha_array_final_50 alpha_array_C3a_C3_50 alpha_array_C3a_C5_50 alpha_array_C5a_C3_50 alpha_array_C5a_C5_50 alpha_array_C3a_C3_C5_50 alpha_array_C5a_C3_C5_50]

# 90% -
alpha_array_C3a_C3_90 = main(index_C3a,[index_C3],[0.1])
alpha_array_C3a_C5_90 = main(index_C3a,[index_C5],[0.1])

alpha_array_C5a_C3_90 = main(index_C5a,[index_C3],[0.1])
alpha_array_C5a_C5_90 = main(index_C5a,[index_C5],[0.1])

alpha_array_C3a_C3_C5_90 = main(index_C3a,[index_C3,index_C5],[0.1,0.1])
alpha_array_C5a_C3_C5_90 = main(index_C5a,[index_C3,index_C5],[0.1,0.1])

alpha_array_final_90 = [alpha_array_final_90 alpha_array_C3a_C3_90 alpha_array_C3a_C5_90 alpha_array_C5a_C3_90 alpha_array_C5a_C5_90 alpha_array_C3a_C3_C5_90 alpha_array_C5a_C3_C5_90]

# 99% -
alpha_array_C3a_C3_99 = main(index_C3a,[index_C3],[0.01])
alpha_array_C3a_C5_99 = main(index_C3a,[index_C5],[0.01])

alpha_array_C5a_C3_99 = main(index_C5a,[index_C3],[0.01])
alpha_array_C5a_C5_99 = main(index_C5a,[index_C5],[0.01])

alpha_array_C3a_C3_C5_99 = main(index_C3a,[index_C3,index_C5],[0.01,0.01])
alpha_array_C5a_C3_C5_99 = main(index_C5a,[index_C3,index_C5],[0.01,0.01])

alpha_array_final_99 = [alpha_array_final_99 alpha_array_C3a_C3_99 alpha_array_C3a_C5_99 alpha_array_C5a_C3_99 alpha_array_C5a_C5_99 alpha_array_C3a_C3_C5_99 alpha_array_C5a_C3_C5_99]

# Write the results files -
writedlm("./robustness/alpha_Z1_50.dat",alpha_array_final_50[:,2:end])
writedlm("./robustness/alpha_Z1_90.dat",alpha_array_final_90[:,2:end])
writedlm("./robustness/alpha_Z1_99.dat",alpha_array_final_99[:,2:end])
