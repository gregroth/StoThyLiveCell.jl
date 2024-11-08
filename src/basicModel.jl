

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct StandardStoModel
    "number of 'promoter' states"
    nbstate::Int64
    "number of parameter, excluding initiation and degradation"
    nbparameters::Int64
    "list of the cartesian indices for which the rate matrix is non zero"
    ParamToRate_idx::Vector{CartesianIndex{2}}
    "list of the parameter indices corresponding to the ParamToRate_idx entries of the rate matrix"
    ParamToRate_val::Vector{Int64}
    "list of the 'promoter' states from which the initiation rate is non zero"
    TrState::Vector{Int64}
end




"""
    StoModel(model::StandardStoModel, parameters::Vector{Float64},kini::Vector{Float64},delta::Float64,maxrna::Int)

Define a model instance based on the topology of the model.
    
#Arguments
- `parameters::Vector{Float64}`: list of the rates of the transition between the 'promoter' states
- `kini::Vector{Float64}`: list of the initiation rates corresponding the list of active 'promoter' states described in model.TrState 
- `delta::Float64`: degradation/release rate
- `maxrna::Int`: maximum number of mRNA considered in the model
"""
function StoModel(model::StandardStoModel, parameters::Vector{Float64},kini::Vector{Float64},delta::Float64, maxrna::Int)
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
            Q[i*model.nbstate+model.TrState[j],i*model.nbstate+model.TrState[j]+model.nbstate] = kini[j]
        end
        for j = 1:model.nbstate
            Q[(i+1)*model.nbstate+j,i*model.nbstate+j] = delta*(i+1)
        end
    end
    for i in axes(Q,1)
        Q[i,i] = Q[i,i] -sum(Q[i,:])
    end
    
    return exp(Q)
end

"""
    StoModel!(P::Array{Float64,2}, model::StandardStoModel, parameters::Vector{Float64},kini::Float64,delta::Float64, maxrna::Int)

 Modify P in the transition matrix of the model model with maximum mRNA number = maxrna, 
    
#Arguments
- `parameters::Vector{Float64}`: list of the rates of the transition between the 'promoter' states
- `kini::Vector{Float64}`: list of the initiation rates corresponding the list of active 'promoter' states described in model.TrState 
- `delta::Float64`: degradation/release rate
- `maxrna::Int`: maximum number of mRNA considered in the model
"""
function StoModel!(P::Array{Float64,2}, model::StandardStoModel, parameters::Vector{Float64},kini::Float64,delta::Float64, maxrna::Int)
    Qstate = zeros(model.nbstate,model.nbstate) 
    for i in eachindex(model.ParamToRate_val)
        Qstate[model.ParamToRate_idx[i]] = parameters[model.ParamToRate_val[i]]
    end
    for i=0:maxrna
        P[i*model.nbstate+1:(i+1)*model.nbstate,i*model.nbstate+1:(i+1)*model.nbstate] .= Qstate
    end

    for i=0:maxrna-1
        for j in eachindex(model.TrState)
            P[i*model.nbstate+model.TrState[j],i*model.nbstate+model.TrState[j]+model.nbstate] = kini[j]
        end
        for j = 1:model.nbstate
            P[(i+1)*model.nbstate+j,i*model.nbstate+j] = delta*(i+1)
        end
    end
    for i in axes(P,1)
        P[i,i] = P[i,i] -sum(P[i,:])
    end
    P .= exp(P)
end



