#list of the distance functions used in the fit of the multi-state model
abstract type AbstractData end
abstract type AbstractDataGroup end

abstract type AbstractDataOptim end

struct DataFit{DT, D} <: AbstractDataOptim
    datatypes::DT
    datagroup::AbstractDataGroup
    data::D
    detectionLimitLC::Int
    detectionLimitNS::Int
    burstsinglet::Symbol
end

struct LiveCellData <: AbstractDataGroup end

struct FixedCellData <: AbstractDataGroup end

struct FixedAndLiveCellData <: AbstractDataGroup end

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

struct Distribution_RNA <: AbstractData

end


