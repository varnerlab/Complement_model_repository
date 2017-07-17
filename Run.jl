include("SolveBalances_Flow.jl")
include("DataFile_Flow.jl")
using PyPlot

TSTART = 0
TSTOP = 24
Ts = 0.01

data_dictionary = DataFile_Flow(TSTART,TSTOP,Ts)

(T,X) = SolveBalances_Flow(TSTART,TSTOP,Ts,data_dictionary)
