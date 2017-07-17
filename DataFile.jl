# ----------------------------------------------------------------------------------- #
# Copyright (c) 2016 Varnerlab
# School of Chemical Engineering Purdue University
# W. Lafayette IN 46907 USA

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #
function DataFile(TSTART,TSTOP,Ts)
# ----------------------------------------------------------------------------------- #
# DataFile.jl was generated using the Kwatee code generation system.
# DataFile: Stores model parameters as key - value pairs in a Julia Dict()
# Username: jeffreyvarner
# Type: GRN-JULIA
# Version: 1.0
# Generation timestamp: 02-26-2016 14:31:22
#
# Input arguments:
# TSTART  - Time start
# TSTOP  - Time stop
# Ts - Time step
#
# Return arguments:
# data_dictionary  - Data dictionary instance (holds model parameters)
# ----------------------------------------------------------------------------------- #

# Initial conditions (default, no initiator)
initial_condition_array = zeros(19)
initial_condition_array[1] = 0.0      # 1 zymosan
initial_condition_array[2] = 1.9      # 2 C4
initial_condition_array[3] = 0.322    # 3 C2
initial_condition_array[8] = 7.57     # 4 C3
initial_condition_array[14] = 0.195   # 5 C5
initial_condition_array[17] = 2.23    # 6 Factor H
initial_condition_array[18] = 0.417   # 7 C4BP

# Load the parameters from disk -
parameter_array = readdlm("./best_solution.txt")

# ---------------------------- DO NOT EDIT BELOW THIS LINE -------------------------- #
data_dictionary = Dict();
data_dictionary["PARAMETER_ARRAY"] = parameter_array;
data_dictionary["INITIAL_CONDITION_ARRAY"] = initial_condition_array;
# ----------------------------------------------------------------------------------- #
return data_dictionary;
end
