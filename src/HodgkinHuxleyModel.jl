module HodgkinHuxleyModel

using Simulation73
using RecursiveArrayTools
using Parameters
using DifferentialEquations: DEDataArray

export HHNetwork, HHPopulation

export RampingBumpStimulus

export DecayingExponentialConnectivity, AMPAandGABASynapses

include("stimulus.jl")
include("connectivity.jl")
include("neuron.jl")
include("network.jl")

end # module
