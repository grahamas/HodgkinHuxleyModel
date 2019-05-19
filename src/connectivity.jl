
abstract type AbstractSynapse{T,W} end
abstract type AbstractSynapses{T,W} end

struct AMPASynapse{T,W} <: AbstractSynapse{T,W}
    E::T
    τ::T
    α::T
    β::T
    γ::T
    delay::T
    connectivity::W
end
# FIXME recent PR allows abstract constructors
function AMPASynapse(; E::T=nothing, delay::T=nothing, g_max::W=nothing, τ::T=nothing) where {T,W}
    α = 1/τ
    AMPASynapse{T,W}(E, τ, α, -2α, -α^2, delay, g_max)
end

struct GABASynapse{T,W} <: AbstractSynapse{T,W}
    E::T
    τ::T
    α::T
    β::T
    γ::T
    delay::T
    connectivity::T
end
function GABASynapse(; E::T=nothing, delay::T=nothing, g_max::W=nothing, τ::T=nothing) where {T,W}
    α = 1/τ
    GABASynapse{T,W}(E, τ, α, -2α, -α^2, delay, g_max)
end
@with_kw struct AMPAandGABASynapses{T,W} <: AbstractSynapses{T,W}
    synapses::AbstractArray{<:Union{AMPASynapse{T,W},GABASynapse{T,W}}}
end
function AMPAandGABASynapses(;kwargs...)
    pops([AMPASynapse GABASynapse;
          AMPASynapse GABASynapse]; kwargs...)
end
function maximal_conductances(synapses::AMPAandGABASynapses, space::AbstractSpace, N_per_point::Tuple{Int,Int})
    assuming_1_per_point = map((syn) -> g_max(syn, space), synapses.synapses)
    map(CartesianIndices(assuming_1_per_point)) do ix
        repeat(assuming_1_per_point[ix], inner=(N_per_point[ix[1]], N_per_point[ix[2]]))
    end
end


@memoize function g_max(syn::AbstractSynapse, space::AbstractSpace)
    directed_weights(syn.connectivity, space)
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
function detect_spikes(state, history, p, t, dt, presynapse_ix)
    @views begin
        Vf_before = history(p, t-dt).x[presynapse_ix].x.x[3]
        Vf_at = history(p, t).x[presynapse_ix].x.x[3]
        Vf_after = history(p, t+dt).x[presynapse_ix].x.x[3]

        last_spike_time = history(p,t).x[presynapse_ix].last_spike_time
        dt_refractory = history(p,t).x[presynapse_ix].dt_refractory
        threshold = history(p,t).x[presynapse_ix].threshold

        @. spikes_bitarr = ((Vf_at > threshold)
                             & (Vf_before <= Vf_at)
                             & (Vf_after <= Vf_at)
                             & (t - last_spike_time) > dt_refractory)
    end
    return spikes_bitarr
end

"Mutate potential of either E or I neuron due to conductance from AMPA"
function mutate_potential!(d_neuron_state, neuron_state, syn::AMPASynapse)
    @views begin
        dV = d_neuron_state.x.x[1]
        Vs = neuron_state.x.x[4]
        g_AMPA = neuron_stae.x.x[5]
        @. d_neuron_state.x.x[1] += g_AMPA * (syn.E * Vs)
    end
end
function mutate_potential!(d_pop_state, pop_state, syn::GABASynapse)
    @views begin
        dV = d_pop_state.x.x[1]
        Vs = pop_state.x.x[4]
        g_GABA = neuron_state.x.x[7]
        @. d_pop_state.x.x[1] += g_GABA * (syn.E * V)
    end
end
function mutate_synaptic_conductance!(d_pop_state, pop_state, syn::AMPASynapse, g_max, presynaptic_spikes)
    @views begin
        d_g, d_z = d_pop_state.x.x[[5,6]]
        g, z = pop_state.x.x[[5,6]]
        @. x = g_max * presynaptic_spikes
        @. d_z += syn.α * x + syn.β * z + syn.γ * g
        @. d_g += z
    end
end
function mutate_synaptic_conductance!(d_pop_state, pop_state, syn::GABASynapse, g_max, presynaptic_spikes)
    @views begin
        d_g, d_z = d_pop_state.x.x[[7.8]]
        g, z = pop_state.x.x[[7.8]]
        @. x = g_max * presynaptic_spikes
        @. d_z += syn.α * x + syn.β * z + syn.γ * g
        @. d_g += z
    end
end

function make_mutator(synapses::AMPAandGABASynapses, space::AbstractSpace, N_per_point::Tuple{Int,Int})
    g_maxes = maximal_conductances(synapses, space, N_per_point)
    function synapses_mutator!(d_state, state, history, p, t)
        @views for pre_ix in 1:size(connectivity,2)
            presynaptic_spikes = detect_spikes(state, history, t-delay, dt, pre_ix)
            for post_ix in 1:size(connectivity,1)
                d_postsynaptic_population = d_state.x[post_ix]
                postsynaptic_population = state.x[post_ix]
                mutate_potential!(d_postsynaptic_population, postsynaptic_population, synapse[post_ix,pre_ix])
                mutate_conductance!(d_postsynaptic_population, postsynaptic_population, synapse[post_ix,pre_ix], g_maxes[post_ix, pre_ix], presynaptic_spikes)
            end
        end
    end
end
