@with_kw struct HHNetwork{T,D} <: AbstractModelwithDelay{T,D,2} 
    space::AbstractSpace{T,D}
    stimulus::AbstractArray{<:AbstractStimulus{T}}
    synapses::AbstractSynapses{T}
    populations::AbstractArray{<:HHPopulation{T}}
end

function Simulation73.initial_value(network::HHNetwork)
    ArrayPartition(
        initial_value(network.neuron_E, network.space),
        initial_value(network.neuron_I, network.space)
    )
end

function Simulation73.history(network::HHNetwork)
    init = initial_value(network); init2 = initial_value(network)
    init .= -∞
    for (i, neuron) ∈ init.x
        neuron.dt_refractory = init2.x[i].dt_refractory
        neuron.threshold = init2.x[i].threshold
    end
end

# FIXME collapse these functions into one with the cardinal difference (V_s vs V) dispatched
# State.last_spike_time is array containing last spike time
# State.x[1] is neuronE
# State[1][1] is V_E
# State[1][2] is n_E
# State[1][3] = V_f (high pass)
# State[1][4] = V_s (low pass)
# State[1][5] is g_AMPA
# State[1][6] is z_AMPA
# State[1][7] is g_GABA
# State[1][8] is z_GABA
# State[2] is neuronI
# ..
function Simulation73.make_system_mutator(network::HHNetwork)
    N_per_point = Tuple(neuron.N_per_point for neuron in network.neurons)
    stimulus_mutator! = make_mutator(network.stimulus, network.space, N_per_point)
    synapses_mutator! = make_mutator(network.synapses, network.space, N_per_point)
    neuron_mutator! = make_mutator(network.populations)
    function system_mutator!(d_state, state, h, p, t)
        # Use nested ArrayPartitions for d_state and state
            # First level, neuron type
            # Second level, conductances, potential etc
        zero!(d_state)
        stimulus_mutator!(d_state, state, t)
        synapses_mutator!(d_state, state, h, p, t)
        neuron_mutator!(d_state, state, t)
    end
end
