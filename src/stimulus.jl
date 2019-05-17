abstract type AbstractStimulus{T,D} <: AbstractParameter{T} end

@memoize Dict function make_mutator(stimulus_arr::AbstractArray{<:AbstractStimulus{T}}, space::AbstractSpace, N_per_point::Tuple{Int,Int}) where T
    stimulus_mutators = [make_stimulus(stim, space, N) for (stim, N) in zip(stimulus_arr, N_per_point)]
    function stimulus_mutator!(d_state::AbstractArray{T,D}, state::AbstractArray{T,D}, t::T) where {T,D}
        @views for (i, stimulus!) in enumerate(stimulus_mutators)
            stimulus!(d_state.x[i].x.x[1], t) # The change in voltage of the ith neuron
        end
    end
end

function distance(x,y) where {T <: Real}
    sqrt(sum((x .- y).^2))
end

struct NoStimulus{T,N} <: AbstractStimulus{T,N} end
function make_stimulus(nostim::NoStimulus{T,N}, space::AbstractSpace{T,N}, P::Int) where {T,N,AT<: AbstractArray{T,N}}
    (val,t) -> return
end

abstract type TransientBumpStimulus{T,N} <: AbstractStimulus{T,N} end
function make_stimulus(bump::TBS, space::AbstractSpace{T,N}, P::Int) where {T, N, TBS<:TransientBumpStimulus{T,N}}
    bump_frame = on_frame(bump, space, P)
    onset = bump.time_window[1]
    offset = bump.time_window[2]
    function stimulus!(val, t)
        if onset <= t < offset
            val .+= bump_frame
        end
    end
end

struct SharpBumpStimulus{T,N} <: TransientBumpStimulus{T,N}
    width::T
    strength::T
    time_window::Tuple{T,T}
    center::NTuple{N,T}
end

function SharpBumpStimulus{T,N}(; strength::T=nothing, width::T=nothing,
        duration=nothing, time_window=nothing, center=NTuple{N,T}(zero(T) for i in 1:N)) where {T,N}
    if time_window == nothing
        return SharpBumpStimulus{T,N}(width, strength, (zero(T), duration), center)
    else
        @assert duration == nothing
        return SharpBumpStimulus{T,N}(width, strength, time_window, center)
    end
end

function on_frame(sbs::SharpBumpStimulus{T,N}, space::AbstractSpace{T,N}, P::Int) where {T,N}
    one_frame = zero(space)
    coords = coordinates(space)
    half_width = sbs.width / 2.0
    one_frame[distance.(coords, Ref(sbs.center)) .<= half_width] .= sbs.strength
    frame = cat((one_frame for i in 1:P)..., dims=N+1)
    return frame
end
