# Library of plot functions -
using PyPlot

function plot_C3a_Control()

  # load the experimental data -
  experimental_data = readdlm("./data/Shaw2015_Fig2a_C3a.txt",',')

  # Plot the data with points -
  plot(experimental_data[:,1],experimental_data[:,2]*0.11002,"ko",LineWidth=4)

end

function plot_C5a_Control()

  # load the experimental data -
  experimental_data = readdlm("./data/Shaw2015_Fig3ai_C5a_original.txt",',')

  # Plot the data with points -
  plot(experimental_data[:,1],experimental_data[:,2]*(0.00009615),"ko",LineWidth=4)

end

function plot_C3a_Z1()

  # load the experimental data -
  experimental_data = readdlm("./data/Shaw2015_Fig2e_C3a.txt",',')

  # Plot the data with points -
  plot(experimental_data[:,1],experimental_data[:,2]*0.11002,"ko",LineWidth=4)

end

function plot_C5a_Z1()

  # load the experimental data -
  experimental_data = readdlm("./data/Shaw2015_Fig3c_C5a_original.txt",',')

  # Plot the data with points -
  plot(experimental_data[:,1],experimental_data[:,2]*(0.00009615),"ko",LineWidth=4)

end

function plot_C3a_Z01()

  # load the experimental data -
  experimental_data = readdlm("./data/Shaw2015_Fig2d_C3a.txt",',')

  # Plot the data with points -
  plot(experimental_data[:,1],experimental_data[:,2]*0.11002,"ko",LineWidth=4)

end

function plot_C5a_Z01()

  # load the experimental data -
  experimental_data = readdlm("./data/Shaw2015_Fig3b_C5a_original.txt",',')

  # Plot the data with points -
  plot(experimental_data[:,1],experimental_data[:,2]*(0.00009615),"ko",LineWidth=4)

end


function plot_C3a_Z001()

  # load the experimental data -
  experimental_data = readdlm("./data/Shaw2015_Fig2c_C3a.txt",',')


  # Plot the data with points -
  plot(experimental_data[:,1],experimental_data[:,2]*0.11002,"ko",LineWidth=4)

end

function plot_C5a_Z001()

  # load the experimental data -
  experimental_data = readdlm("./data/Shaw2015_Fig3aiii_C5a_original.txt",',')

  # Plot the data with points -
  plot(experimental_data[:,1],experimental_data[:,2]*(0.00009615),"ko",LineWidth=4)

end

function plot_C3a_Z0001()

  # load the experimental data -
  experimental_data = readdlm("./data/Shaw2015_Fig2b_C3a.txt",',')


  # Plot the data with points -
  plot(experimental_data[:,1],experimental_data[:,2]*0.11002,"ko",LineWidth=4)

end

function plot_C5a_Z0001()

  # load the experimental data -
  experimental_data = readdlm("./data/Shaw2015_Fig3aii_C5a_original.txt",',')

  # Plot the data with points -
  plot(experimental_data[:,1],experimental_data[:,2]*(0.00009615),"ko",LineWidth=4)

end
