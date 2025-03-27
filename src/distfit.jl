#list of the distance functions used in the fit of the multi-state model
abstract type AbstractDistanceFitRNA end
abstract type AbstractDistanceFitBurst end
abstract type AbstractDistanceFitRNAandBurst end



struct LikelihoodRNA <: AbstractDistanceFitRNA 
    weight::Float32
end
struct LikelihoodConvolutedRNA <: AbstractDistanceFitRNA 
    weight::Float32
end

struct LikelihoodConvolutedConditionalRNA <: AbstractDistanceFitRNA 
    weight::Float32
    condition::Int64 # the condition is #mRNA>=condition
end

struct LsqSurvival <: AbstractDistanceFitBurst 
    weight::Float32
end

"weighted average of the lsq on the linear data and the log data"
struct LsqSurvivalLogLin <: AbstractDistanceFitBurst 
    weight::Float32
    plog::Float32
end


struct LsqProb <: AbstractDistanceFitBurst 
    weight::Float32
end


struct LsqNumber <: AbstractDistanceFitBurst 
    weight::Float32
end

struct LsqIntensity <: AbstractDistanceFitBurst
    weight::Float32
end


function opt_dist(estimate_signal::Vector, ref_signal::Vector, dist::AbstractDistanceFitRNA; kwargs...)
    dist(estimate_signal, ref_signal; kwargs...)
end


function opt_dist(estimate_signal::Vector, ref_signal::Vector, dist::AbstractDistanceFitBurst; kwargs...)
    dist(estimate_signal, ref_signal; kwargs...)
end


function (f::LikelihoodRNA)(estimate_signal::AbstractVector{T}, ref_signal::Vector; kwargs...) where T
    estimate_signal_ = @. max(estimate_signal, 0)
    ind = @. Int(floor(ref_signal) + 1) # because the index is from 0 
    -sum(log.(estimate_signal_[ind] .+ 1e-6))*f.weight # add a small quantity to avoid log(0)
end

function (f::LikelihoodConvolutedRNA)(estimate_signal::AbstractVector{T}, ref_signal::Vector; kwargs...) where T
    estimate_signal_ = @. max(estimate_signal, 0)
     #convolution of two independend distributions
     conv_result = zeros(Float64, 2*length(estimate_signal_) - 1)
     for i in eachindex(estimate_signal_)
         for j in eachindex(estimate_signal_)
             conv_result[i + j - 1] += estimate_signal_[i] * estimate_signal_[j]
         end
     end
    ind = @. Int(floor(ref_signal) + 1) # because the index is from 0 
    -sum(log.(conv_result[ind] .+ 1e-6))*f.weight # add a small quantity to avoid log(0)
end

function (f::LikelihoodConvolutedConditionalRNA)(estimate_signal::AbstractVector{T}, ref_signal::Vector; kwargs...) where T
     estimate_signal_ = @. max(estimate_signal, 0)
     #convolution of two independend distributions
     conv_result = zeros(Float64, 2*length(estimate_signal_) - 1)
     for i in eachindex(estimate_signal_)
         for j in eachindex(estimate_signal_)
             conv_result[i + j - 1] += estimate_signal_[i] * estimate_signal_[j]
         end
     end
    condDist = conv_result[f.condition+1:end]./sum(conv_result[f.condition+1:end]) # +1 is because the first index is 0
    ind = @. Int(floor(ref_signal) + 1) # because the index is from 0 
    -sum(log.(condDist[ind] .+ 1e-6))*f.weight # add a small quantity to avoid log(0)
end

function (f::LsqSurvival)(estimate_signal::AbstractVector{T}, ref_signal::Tuple{Vector,Vector}; kwargs...) where T
    f.weight * sum((log.(estimate_signal .+ 1e-6) - log.(ref_signal[2].+ 1e-6)).^2)/length(estimate_signal)
end

function (f::LsqSurvivalLogLin)(estimate_signal::AbstractVector{T}, ref_signal::Tuple{Vector,Vector}; kwargs...) where T
    f.weight * (f.plog * sum((log.(estimate_signal .+ 1e-6) - log.(ref_signal[2].+ 1e-6)).^2) + (1-f.plog) * sum(((estimate_signal) - (ref_signal[2])).^2) )/length(estimate_signal)
end

function (f::LsqProb)(estimate_signal_tot::T, ref_signal::Float64; kwargs...)where T
    f.weight * (log.(estimate_signal_tot.+ 1e-6) - log.(ref_signal.+ 1e-6)).^2
end

function (f::LsqNumber)(estimate_signal_tot::T, ref_signal::Float64; kwargs...) where T
    f.weight * (estimate_signal_tot - ref_signal)^2
end

function (f::LsqIntensity)(estimate_signal::AbstractVector{T}, ref_signal::Tuple{Vector,Vector}; kwargs...) where T
    f.weight * sum((estimate_signal - ref_signal[2]).^2)/length(estimate_signal)
end



