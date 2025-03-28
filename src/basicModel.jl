

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct StandardStoModel
    "number of 'promoter' states"
    nbstate::Int64
    "number of parameter, excluding initiation and degradation"
    nbparameters::Int64
    "number of initiation rates"
    nbkini::Int64
    "list of the cartesian indices for which the rate matrix is non zero"
    ParamToRate_idx::Vector{CartesianIndex{2}}
    "list of the parameter indices corresponding to the ParamToRate_idx entries of the rate matrix"
    ParamToRate_val::Vector{Int64}
    "list of the 'promoter' states from which the initiation rate is non zero"
    TrState::Vector{Int64}
    "list of indices of initiation rate in the active states; same size than TrState"
    KiniToRate_idx::Vector{Int64}
    "index of the degradation rate in the parameter vector"
    deltaToRate_idx::Int64
end




"""
    StoModel(model::StandardStoModel, parameters::Vector{Float64},maxrna::Int)

Define a model instance based on the topology of the model.
    
#Arguments
- `parameters::Vector{Float64}`: list of the all the rate parameters
- `maxrna::Int`: maximum number of mRNA considered in the model
"""
function StoModel(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int)
    Qstate = zeros(model.nbstate,model.nbstate) 
    Q = zeros(model.nbstate*(maxrna+1),model.nbstate*(maxrna+1))
    for i in eachindex(model.ParamToRate_val)
        Qstate[model.ParamToRate_idx[i]] = parameters[model.ParamToRate_val[i]]
    end
    for i=0:maxrna
        Q[i*model.nbstate+1:(i+1)*model.nbstate,i*model.nbstate+1:(i+1)*model.nbstate] .= Qstate
    end

    for i=0:maxrna-1
        for j in eachindex(model.TrState)
            Q[i*model.nbstate+model.TrState[j],i*model.nbstate+model.TrState[j]+model.nbstate] = parameters[model.KiniToRate_idx[j]]
        end
        for j = 1:model.nbstate
            Q[(i+1)*model.nbstate+j,i*model.nbstate+j] = parameters[model.deltaToRate_idx]*(i+1)
        end
    end
    for i in axes(Q,1)
        Q[i,i] = Q[i,i] -sum(Q[i,:])
    end
    
    return exp(Q)
end

"""
    StoModelAD(model::StandardStoModel, parameters::Vector{Float64},maxrna::Int)

Define a model instance based on the topology of the model, specifically for AutoDiff optimization fcts
    
#Arguments
- `parameters::Vector{Float64}`: list of the all the rate parameters
- `maxrna::Int`: maximum number of mRNA considered in the model
"""
function StoModelAD(model::StandardStoModel, parameters::AbstractVector{T}, maxrna::Int) where T
    Qstate = zeros(T, model.nbstate,model.nbstate) 
    Q = zeros(T, model.nbstate*(maxrna+1),model.nbstate*(maxrna+1))
    for i in eachindex(model.ParamToRate_val)
        Qstate[model.ParamToRate_idx[i]] = parameters[model.ParamToRate_val[i]]
    end
    for i=0:maxrna
        Q[i*model.nbstate+1:(i+1)*model.nbstate,i*model.nbstate+1:(i+1)*model.nbstate] .= Qstate
    end

    for i=0:maxrna-1
        for j in eachindex(model.TrState)
            Q[i*model.nbstate+model.TrState[j],i*model.nbstate+model.TrState[j]+model.nbstate] = parameters[model.KiniToRate_idx[j]]
        end
        for j = 1:model.nbstate
            Q[(i+1)*model.nbstate+j,i*model.nbstate+j] = parameters[model.deltaToRate_idx]*(i+1)
        end
    end
    for i in axes(Q,1)
        Q[i,i] = Q[i,i] -sum(Q[i,:])
    end
    
    return exponential!(copyto!(similar(Q), Q), ExpMethodGeneric())
end

"""
    StoModel!(P::Array{Float64,2}, model::StandardStoModel, parameters::Vector{Float64},kini::Float64,delta::Float64, maxrna::Int)

 Modify P in the transition matrix of the model model with maximum mRNA number = maxrna, 
    
#Arguments
- `parameters::Vector{Float64}`: list of the all the rate parameters
- `maxrna::Int`: maximum number of mRNA considered in the model
"""
function StoModel!(P::Array{Float64,2}, model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int)
    Qstate = zeros(model.nbstate,model.nbstate) 
    for i in eachindex(model.ParamToRate_val)
        Qstate[model.ParamToRate_idx[i]] = parameters[model.ParamToRate_val[i]]
    end
    for i=0:maxrna
        P[i*model.nbstate+1:(i+1)*model.nbstate,i*model.nbstate+1:(i+1)*model.nbstate] .= Qstate
    end

    for i=0:maxrna-1
        for j in eachindex(model.TrState)
            P[i*model.nbstate+model.TrState[j],i*model.nbstate+model.TrState[j]+model.nbstate] = parameters[model.KiniToRate_idx[j]]
        end
        for j = 1:model.nbstate
            P[(i+1)*model.nbstate+j,i*model.nbstate+j] = parameters[model.deltaToRate_idx]*(i+1)
        end
    end
    for i in axes(P,1)
        P[i,i] = P[i,i] -sum(P[i,:])
    end
    P .= exp(P)
end



"""
    StoModel_RateMat(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int)

Return the rate matrix of a model instance based on the topology of the model.

#Arguments
- `parameters::Vector{Float64}`: list of the all the rate parameters
- `maxrna::Int`: maximum number of mRNA considered in the model
"""
function StoModel_RateMat(model::StandardStoModel, parameters::AbstractVector{T}, maxrna::Int) where T
    Qstate = zeros(T, model.nbstate,model.nbstate) 
    Q = zeros(T, model.nbstate*(maxrna+1),model.nbstate*(maxrna+1))
    for i in eachindex(model.ParamToRate_val)
        Qstate[model.ParamToRate_idx[i]] = parameters[model.ParamToRate_val[i]]
    end
    for i=0:maxrna
        Q[i*model.nbstate+1:(i+1)*model.nbstate,i*model.nbstate+1:(i+1)*model.nbstate] .= Qstate
    end

    for i=0:maxrna-1
        for j in eachindex(model.TrState)
            Q[i*model.nbstate+model.TrState[j],i*model.nbstate+model.TrState[j]+model.nbstate] = parameters[model.KiniToRate_idx[j]]
        end
        for j = 1:model.nbstate
            Q[(i+1)*model.nbstate+j,i*model.nbstate+j] = parameters[model.deltaToRate_idx]*(i+1)
        end
    end
    for i in axes(Q,1)
        Q[i,i] = Q[i,i] -sum(Q[i,:])
    end
    
    return Q
end