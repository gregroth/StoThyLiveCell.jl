#list of the distance functions used in the fit of the multi-state model
abstract type AbstractData end
abstract type AbstractDataOptim end

struct DataFit{DT, DG, D} <: AbstractDataOptim
    datatypes::DT
    datagroups::DG
    data::D
    detectionLimitLC::Int
    detectionLimitNS::Int
end

struct Survival_Burst <: AbstractData
end

struct Survival_InterBurst <: AbstractData

end

struct Survival_NextBurst <: AbstractData

end

struct Prob_Burst <: AbstractData

end

struct Mean_Nascent <: AbstractData

end

struct Intensity_Burst <: AbstractData

end

struct Correlation_InterBurst <: AbstractData

end



