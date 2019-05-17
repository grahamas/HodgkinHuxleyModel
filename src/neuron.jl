
@with_kw struct HHNeuron{T}
    C::T
    E_K::T
    E_L::T
    E_Na::T
    SA::T
    V_m_half::T
    V_n_half::T
    g_K::T
    g_L::T
    g_Na::T
    k_m::T
    k_n::T
    τ_n::T
    τ_V_f::T
    τ_V_s::T
    dt_refractory::T,
    threshold::T
    N_per_point::T=1
end

mutable struct NeuronData{T,N} <: DEDataArray{T,N}
    x::ArrayPartition{T,N}
    last_spike_time::Array{T,N}
    dt_refractory::T
    threshold::T
end

function zero(neuron::HHNeuron{T}, space::AbstractSpace{T,D}) where {T,D}
    cat([zero(space) for i in 1:neuron.N_per_point]..., dims=D+1)
end
initial_value(neuron::HHNeuron, space::AbstractSpace) = NeuronData(
    ArrayPartition(
        [zero(neuron, space) for i in 1:8]... # [V, n, Vf, Vs, gAMPA, zA, gGABA, gG]
    ), # Differential data
    zero(neuron, space) .- ∞, # last_spike_time
    neuron.dt_refractory,
    neuron.threshold
)

function Θ(V, V_half, K)
    1 / (1 + exp((V_half - V) / K))
end

function make_mutator(neuron::HHNeuron)
    @unpack_HHNeuron neuron # FIXME I am dangerous
    function neuron_mutator!(d_neuron_state, neuron_state, t)
        # Gets neuron_state of single neuron, so neuron_state[1] is V, neuron_state[2] is n,
        # neuron_state[3] = V_f (high pass), neuron_state[4] = V_s (low pass)
        @views begin
            dV, dn, dV_f, dV_s = d_neuron_state.x.x[1:4]
            V, n, V_f, V_s = neuron_state.x.x[1:4]
            @. n_∞ = Θ(V, V_n_half, k_n)
            @. m_∞ = Θ(V, V_m_half, k_m)
            @. dV += g_L * (E_L - V)
            @. dV += g_Na * m_∞ * (E_Na - V)
            @. dV += g_K * n * (E_K - V)
            @. dn += (n_∞ - n) / τ_n
            @. dV_f += dV - (V_f / τ_V_f)
            @. dV_s += (V - V_s) / τ_V_s
        end
    end
end
