module HodgkinHuxleyModel

using Simulation73
using RecursiveArrayTools
using Parameters

include("stimulus.jl")
include("connectivity.jl")
include("neuron.jl")
include("network.jl")

end # module
