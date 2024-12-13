#list of the distance functions used in the fit of the multi-state model
abstract type AbstractDataLiveCell end
abstract type AbstractDataFixedCell end

abstract type AbstractDataOptim end

struct DataFit{DT, DG, D} <: AbstractDataOptim
    datatypes::DT
    data::D
    detectionLimitLC::Int
    detectionLimitNS::Int
end

struct Survival_Burst <: AbstractDataLiveCell
end

struct Survival_InterBurst <: AbstractDataLiveCell

end

struct Survival_NextBurst <: AbstractDataLiveCell

end

struct Prob_Burst <: AbstractDataLiveCell

end

struct Mean_Nascent <: AbstractDataLiveCell

end

struct Intensity_Burst <: AbstractDataLiveCell

end

struct Correlation_InterBurst <: AbstractDataLiveCell

end

struct Distribution_RNA <: AbstractDataFixedCell

end


