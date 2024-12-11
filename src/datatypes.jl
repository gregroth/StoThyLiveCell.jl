#list of the distance functions used in the fit of the multi-state model
abstract type AbstractData end
abstract type AbstractDataOptim end

struct DataFit{DT, DG, D, DI} <: AbstractDataOptim
    datatypes::DT
    datagroups::DG
    data::D
    dist::DI
    maxrnaLC::Int
    maxrnaFC::Int
    detectionLimitLC::Int
    detectionLimitFC::Int
end

struct Survival_Burst{D} <: AbstractData
    data::D
end

struct Survival_InterBurst{D} <: AbstractData
    data::D
end

struct Survival_NextBurst{D} <: AbstractData
    data::D
end

struct Prob_Burst{D} <: AbstractData
    data::D
end

struct Intensity_Burs{D} <: AbstractData
    data::D
end

struct Correlation_InterBurst{D} <: AbstractData
    data::D
end





function (f::Survival_Burst)(tmax::Int, optimstuct::OptimStructWrapper)
    @unpack utilemats = optimstruct
    
    @unpack stateTr, stateTr_on, stateAbs_on, weightsTr_off, P, ssp, PabsOff = utilemats


    weightsAbs_on = ssp[stateAbs_on]./sum(ssp[stateAbs_on])     
    weightsTr_on = weightsAbs_on' * P[stateAbs_on,stateTr_on]./sum(weightsAbs_on' * P[stateAbs_on,stateTr_on])

    PabsOn = P[stateTr_on,stateTr_on]
    tempdist = weightsTr_on
    survivalspot_model_full = Vector{Float64}(undef,tmax)
    for i in 1:tmax
        tempdist = tempdist* PabsOn
        survivalspot_model_full[i] = sum(tempdist)
    end 
    return survivalspot_model_full
end

function (f::Survival_InterBurst)(tmax::Int, optimstuct::OptimStructWrapper)
    @unpack utilemats = optimstruct
    
    @unpack stateTr, stateTr_on, stateAbs_on, weightsTr_off, P, ssp, PabsOff = utilemats

    tempdist = weightsTr_off
    survivaldark_model_full = Vector{Float64}(undef,tmax)
    for i in 1:tmax
        tempdist = PabsOff'*tempdist
        survivaldark_model_full[i] = sum(tempdist)
    end 
    return survivaldark_model_full
end

function (f::Survival_NextBurst)(tmax::Int, optimstuct::OptimStructWrapper)
    @unpack utilemats = optimstruct
    
    @unpack stateTr, stateTr_on, stateAbs_on, weightsTr_off, P, ssp, PabsOff = utilemats

    tempdist = sspTr_off'
    survivalnextburst_model = Vector{Float64}(undef,tmax)
    for i in 1:tmax
        tempdist = tempdist* PabsOff
        survivalnextburst_model[i] = sum(tempdist)
    end 
    return survivalnextburst_model
end


function (f::Prob_Burst)(optimstuct::OptimStructWrapper)
    @unpack utilemats = optimstruct
    
    @unpack stateTr, stateTr_on, stateAbs_on, weightsTr_off, P, ssp, PabsOff = utilemats

    return sum(ssp[stateTr_on])

end


function (f::Intensity_Burst)(tmax::Int, optimstuct::OptimStructWrapper)
    @unpack utilemats = optimstruct
    
    @unpack stateTr, stateTr_on, stateAbs_on, weightsTr_off, P, ssp, PabsOff = utilemats

    weightsON = normalizemat!(P[stateTr_off,stateAbs_off]'*sspTr_off)
 
    rnanbvec_on = vcat(kron([x for x in detectionlimitLC : maxrna],ones(nbstate)))

    intensitytemp = weightsON
    intensity_model = Vector{Float64}(undef,tmax)
    for i in 1:tmax
        intensity_model[i] = (rnanbvec_on'*intensitytemp)[1]
        intensitytemp =  Pabs'*intensitytemp
    end 
    return intensity_model./maximum(intensity_model)
end


function (f::Correlation_InterBurst)(tmax::Int, optimstuct::OptimStructWrapper)
    @unpack utilemats = optimstruct
    
    @unpack stateTr, stateTr_on, stateAbs_on, weightsTr_off, P, ssp, PabsOff = utilemats

    #correlation of the interburst durations
    Qn = P[stateAbs_on,stateAbs_on]
    Rn = P[stateAbs_on,stateTr_on]

    Qb = P[stateTr_on,stateTr_on]
    Rb = P[stateTr_on,stateAbs_on]
    c = ones(length(stateAbs_on))

    Nn = (I - Qn)^(-1)
    Nb = (I - Qb)^(-1)

    NR = Nb*Rb
    Nc = Nn*c

    cortemp=0
    wpre = weightsTr_off
    wpre2 = Rn'*wpre./sum(wpre)
    wpre3 = NR'*wpre2/sum(wpre2)
    ET2t = Nc'*wpre3 

    for t=1:tmax
        cortemp = cortemp + t*ET2t[1]*sum(Rn'*wpre)
        wpre = Qn'*wpre
        wpre2 = Rn'*wpre./sum(wpre)
        wpre3 = NR'*wpre2./sum(wpre2)
        ET2t = Nc'*wpre3
        if sum(wpre)<1e-6
            break
        end
    end
    Et1 = Nc'weightsTr_off 
    M2T = Nc'*(2*Nn'-I)*weightsTr_off
    VarT = M2T[1] - Et1[1]^2

    return (cortemp-Et1[1]^2)/VarT
end

