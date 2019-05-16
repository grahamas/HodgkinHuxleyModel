@with_kw struct HHNetwork{T} <: AbstractModel{T,D,2}
    space::AbstractSpace{T,D}
    stimulus::AbstractArray{<:AbstractStimulus{T}}
    synapse_AMPA::AMPASynapse{T}
    synapse_GABA::GABASynapse{T}
    neuron_E::HHNeuron{T}
    neuron_I::HHNeuron{T}
end

function initial_value(network::HHNetwork)
    ArrayPartition(
        initial_value(network.neuron_E, network.space),
        initial_value(network.neuron_I, network.space)
    )
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
    stimulus_mutator! = make_mutator(network.stimulus, network.space)
    synapse_AMPA_mutator! = make_mutator(network.synapse_AMPA)
    synapse_GABA_mutator! = make_mutator(network.synapse_GABA)
    neuron_E_mutator! = make_mutator(network.neuron_E)
    neuron_I_mutator! = make_mutator(network.neuron_I)
    function system_mutator!(d_state, state, h, p, t)
        # Use nested ArrayPartitions for d_state and state
            # First level, neuron type
            # Second level, conductances, potential etc
        zero!(d_state)
        stimulus_mutator!(d_state, state, t)
        synapse_AMPA_mutator!(d_state, state, h, p, t)
        synapse_GABA_mutator!(d_state, state, h, p, t)
        neuron_E_mutator!(d_state.x[1], state.x[1], t)
        neuron_I_mutator!(d_state.x[2], state.x[2], t)
    end
end
