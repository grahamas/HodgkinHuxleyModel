module HodgkinHuxleyModel

using Simulation73
using RecursiveArrayTools
using Parameters
using DifferentialEquations: DEDataArray
using BioNeuralNetworkModels

export HHNetwork, HHPopulation

export RampingBumpStimulus

export AMPAandGABASynapses

include("stimulus.jl")
include("connectivity.jl")
include("neuron.jl")
include("network.jl")

end # module
