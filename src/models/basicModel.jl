using LinearAlgebra
using ExponentialUtilities


"""
$(TYPEDEF)
"""
struct StandardStoModel
    nbstate::Int
    nbparameters::Int
    ParamToRate::Array{Float64,2}
    TrState::Vector{Int}
end


"""
$(TYPEDEF)
"""
struct StoModel
    P::Array{Float64,2} #transition probability matrix (nbstate x (maxrna+1))^2
end

"""
    StoModel(model::StandardStoModel, parameters::Vector{Float64},kini::Float64,delta::Float64)

Define a model instance based on the structured of the model, 
the parameters (transition rate between states), and the initiation and degrradation rate
"""
function StoModel(model::StandardStoModel, parameters::Vector{Float64},kini::Float64,delta::Float64, maxrna::Int)
    Qstate = zeros(model.nbstate) 
    Q = zeros(model.nbstate*(maxrna+1))
    for (idx,val) in enumerate(model.ParamToRate)
        Qstate[idx] = parameters[val]
    end
    for i=0:maxrna
        Q[i*model.nbstate+1:(i+1)*model.nbstate,i*model.nbstate+1:(i+1)*model.nbstate] .= Qstate
    end
    for i=1:model.nbstate*(maxrna)
        Q[i,model.nbstate+i] = kini
        Q[model.nbstate+i,i] = delta
    end 
    
    Q = exp(Q)
 
    StoModel(kini,delta,P)
end


