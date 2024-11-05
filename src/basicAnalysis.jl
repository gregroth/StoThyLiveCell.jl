"""
    ModelOutput(model::StandardStoModel, parameters::Vector{Float64},kini::Float64,delta::Float64, maxrna::Int,tmaxon,tmaxoff,tmaxnextburst,tmaxintensity)

Return the main outputs of the model: mean nb of nascent mrna, probability to observe a burst, ON time survival,
OFF time survival,next burst time survival, correlation of consecutive inter burst evts, average intensity track
"""
function ModelOutput(model::StandardStoModel, parameters::Vector{Float64},kini::Float64,delta::Float64, maxrna::Int,tmaxon,tmaxoff,tmaxnextburst,tmaxintensity) 
    timevec_on = 1:1:tmaxon
    timevec_off = 1:1:tmaxoff
    timevec_nextburst = 1:1:tmaxnextburst
    timevec_intensity = 1:1:tmaxintensity
    #model instance
    P = StoModel(model::StandardStoModel, parameters::Vector{Float64},kini::Float64,delta::Float64, maxrna::Int)
    
    #nascent mrna
    evs = eigvecs(P')
    ssp = real.(evs[:,end]./sum(evs[:,end]))
    stateBns = [x for x in 2*model.nbstate+1 :(maxrna+1)*model.nbstate]

    pB = sum(ssp[stateBns])
    prna = ssp'kron(diagm(ones(maxrna+1)), ones(model.nbstate))
    mnascentmrna_model = [x for x in 2 : maxrna]'prna[3:end]./pB
    
    #on times
    stateTr_on = [x for x in model.nbstate+1 :(maxrna+1)*model.nbstate]
    stateAbs_on = [x for x in 1 : model.nbstate]
    
    weightsAbs_on = ssp[stateAbs_on]./sum(ssp[stateAbs_on])     
    weightsTr_on = weightsAbs_on' * P[stateAbs_on,stateTr_on]./sum(weightsAbs_on' * P[stateAbs_on,stateTr_on])
    
    PabsOn = P[stateTr_on,stateTr_on]
    tempdist = weightsTr_on
    survivalspot_model_full = Vector{Float64}(undef,tmaxon)
    for i in 1:maximum(timevec_on)
        tempdist = tempdist* PabsOn
        survivalspot_model_full[i] = sum(tempdist)
    end 
  
    #off times    
    stateAbs_off = [x for x in model.nbstate+1 :(maxrna+1)*model.nbstate]
    stateTr_off = [x for x in 1 : model.nbstate]
    
    weightsAbs_off = ssp[stateAbs_off]./sum(ssp[stateAbs_off])
    weightsTr_off = weightsAbs_off' * P[stateAbs_off,stateTr_off]./sum(weightsAbs_off' * P[stateAbs_off,stateTr_off])
    
    PabsOff = P[stateTr_off,stateTr_off]
    tempdist = weightsTr_off
    survivaldark_model_full = Vector{Float64}(undef,tmaxoff)
    for i in 1:maximum(timevec_off)
        tempdist = tempdist* PabsOff
        survivaldark_model_full[i] = sum(tempdist)
    end 

    #next burst times
    sspTr_off =ssp[stateTr_off]./sum(ssp[stateTr_off]) 
    tempdist = sspTr_off'
    survivalnextburst_model = Vector{Float64}(undef,tmaxnextburst)
    for i in 1:maximum(timevec_nextburst)
        tempdist = tempdist* PabsOff
        survivalnextburst_model[i] = sum(tempdist)
    end 

    #probability to observer a burst
    pburst_model = sum(ssp[stateTr_on]) 

    #correlation of the interburst durations
    Qn = P[stateAbs_on,stateAbs_on]
    Rn = P[stateAbs_on,stateTr_on]

    Qb = P[stateTr_on,stateTr_on]
    Rb = P[stateTr_on,stateAbs_on]
    c = ones(length(stateAbs_on))

    Nn = (I - Qn)^(-1)
    Nb = (I - Qb)^(-1)

    cortemp=0
    wpre = weightsTr_off
    for t=1:15000
        wpre2 = wpre*Rn./sum(wpre)
        wpre3 = wpre2*Nb*Rb./sum(wpre2)
        ET2t = wpre3*Nn*c 
        cortemp = cortemp + t*ET2t[1]*sum(wpre*Rn)
        wpre = wpre*Qn
        if sum(wpre)<1e-6
            break
        end
    end
    Et1 = weightsTr_off*Nn*c 
    M2T = weightsTr_off*(2*Nn-I)*Nn*c
    VarT = M2T[1] - Et1[1]^2
    
    corr_interburst = (cortemp-Et1[1]^2)/VarT

    #average intensity track 
    weightsON = sspTr_off' * P[stateTr_off,stateAbs_off]./sum(sspTr_off' * P[stateTr_off,stateAbs_off])
 
    rnanbvec_on = vcat(kron([x for x in 1:maxrna],ones(model.nbstate)))

    Pabs = Qb

    intensitytemp = weightsON
    intensity_model = Vector{Float64}(undef,length(timevec_intensity))
    for i in eachindex(timevec_intensity)
        intensity_model[i] = (intensitytemp*rnanbvec_on)[1]
        intensitytemp = intensitytemp* Pabs
    end 

    return  mnascentmrna_model, pburst_model, survivalspot_model_full,survivaldark_model_full, survivalnextburst_model, corr_interburst, intensity_model
end
    
