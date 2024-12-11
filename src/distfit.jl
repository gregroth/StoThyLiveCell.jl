#list of the distance functions used in the fit of the multi-state model
abstract type AbstractDistanceFitRNA end
abstract type AbstractDistanceFitBurst end
abstract type AbstractDistanceFitRNAandBurst end



struct LikelihoodRNA <: AbstractDistanceFitRNA end

struct LsqSurvival <: AbstractDistanceFitBurst end

struct LsqProb <: AbstractDistanceFitBurst end


struct LsqNumber <: AbstractDistanceFitBurst end



function opt_dist(estimate_signal::Vector, ref_signal::Vector, dist::AbstractDistanceFitRNA; kwargs...)
    dist(estimate_signal, ref_signal; kwargs...)
end


function opt_dist(estimate_signal::Vector, ref_signal::Vector, dist::AbstractDistanceFitBurst; kwargs...)
    dist(estimate_signal, ref_signal; kwargs...)
end


function (f::LikelihoodRNA)(estimate_signal::Vector, ref_signal::Vector; kwargs...)
    estimate_signal_ = @. max(estimate_signal, 0)
    ind = @. Int(floor(ref_signal) + 1) # because the index is from 0 
    -sum(log.(estimate_signal_[ind] .+ 1e-6)) # add a small quantity to avoid log(0)
end

function (f::AbstractDistanceFitBurst)(estimate_signal_tot::Vector, ref_signal::Array, upperbound::Int; kwargs...)
    data_cutoff = findlast(ref_signal[:,1] .<=upperbound)
    estimate_signal = estimate_signal_tot[Int.(ref_signal[1:data_cutoff,1])]
    sum((log.(estimate_signal) - log.(ref_signal[1:data_cutoff,2])).^2)/length(estimate_signal)
end

function (f::LsqProb)(estimate_signal_tot::Float64, ref_signal::Float64, upperbound::Int; kwargs...)
    (log.(estimate_signal_tot) - log.(ref_signal)).^2
end

function (f::LsqNumber)(estimate_signal_tot::Float64, ref_signal::Float64, upperbound::Int; kwargs...)
    (estimate_signal_tot - ref_signal)^2
end





