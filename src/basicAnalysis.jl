"""
    normalizemat!(A::Array{Float64,2})

return the matrix normilized by its sum
"""
function normalizemat!(A::Array{Float64,2})
    A .= A./sum(A)
end

"""
    normalizemat!(A::Vector{Float64})

return the vector normilized by its sum
"""
function normalizemat!(A::Vector{Float64})
    A .= A./sum(A)
end
"""
    ModelOutput(model::StandardStoModel, parameters::Vector{Float64},kini::Float64,delta::Float64, maxrna::Int64,tmaxon,tmaxoff,tmaxnextburst,tmaxInt64ensity)

Return the main outputs of the model: mean nb of nascent mrna, probability to observe a burst, ON time survival,
OFF time survival,next burst time survival, correlation of consecutive inter burst evts, average intensity track
"""
function ModelOutput(model::StandardStoModel, parameters::Vector{Float64},kini::Float64,delta::Float64, maxrna::Int64,tmaxon::Int64,tmaxoff::Int64,tmaxnextburst::Int64,tmaxintensity::Int64) 
    timevec_on = 1:1:tmaxon
    timevec_off = 1:1:tmaxoff
    timevec_nextburst = 1:1:tmaxnextburst
    timevec_intensity = 1:1:tmaxintensity
    #model instance
    P = StoModel(model::StandardStoModel, parameters::Vector{Float64},kini::Float64,delta::Float64, maxrna::Int64)
    
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
    intensity_model .= intensity_model./intensity_model[1] 
    return  mnascentmrna_model, pburst_model, survivalspot_model_full,survivaldark_model_full, survivalnextburst_model, corr_interburst, intensity_model
end


"""
    mo_basics(model::StandardStoModel, parameters::Vector{Float64},kini::Float64,delta::Float64, maxrna::Int64)

return important vectors and matrices used in the analysis of the model
"""
function mo_basics(model::StandardStoModel, parameters::Vector{Float64},kini::Float64,delta::Float64, maxrna::Int64) 
    P = StoModel(model::StandardStoModel, parameters::Vector{Float64},kini::Float64,delta::Float64, maxrna::Int64)
    evs = eigvecs(P')
    ssp = normalizemat(real.(evs[:,end]))
    stateTr = [x for x in 2*model.nbstate+1 :(maxrna+1)*model.nbstate]
    stateTr_on = [x for x in model.nbstate+1 :(maxrna+1)*model.nbstate]
    stateAbs_on = [x for x in 1 : model.nbstate]
    weightsAbs_off = normalizemat(ssp[stateTr_on])
    weightsTr_off = normalizemat(weightsAbs_off' * P[stateTr_on,stateAbs_on])
    PabsOff = P[stateAbs_on,stateAbs_on]
    sspTr_off =normalizemat(ssp[stateAbs_on])
    Pabs = P[stateTr_on,stateTr_on]
    return P,ssp, stateTr, stateTr_on, stateAbs_on, weightsTr_off,PabsOff, sspTr_off, Pabs
end



"""
    mo_basics!(model::StandardStoModel, parameters::Vector{Float64},kini::Float64,delta::Float64, maxrna::Int64, P::Array{Float64,2},ssp::Vector{Float64}, stateTr::Vector{Int64}, stateTr_on::Vector{Int64}, stateAbs_on::Vector{Int64}, weightsTr_off::Array{Float64,2},PabsOff::Array{Float64,2})

change vector and matrices used in on and off time; mean nascent rna, p_on, correlation
"""
function mo_basics!(model::StandardStoModel, parameters::Vector{Float64},kini::Float64,delta::Float64, maxrna::Int64, P::Array{Float64,2},ssp::Vector{Float64}, stateTr_on::Vector{Int64}, stateAbs_on::Vector{Int64}, weightsTr_off::Array{Float64,2},PabsOff::Array{Float64,2}) 
    P .= StoModel(model, parameters,kini,delta, maxrna)
    evs = eigvecs(P')
    ssp .= normalizemat(real.(evs[:,end]))
    weightsAbs_off = normalizemat(ssp[stateTr_on])
    weightsTr_off .= normalizemat(weightsAbs_off' * P[stateTr_on,stateAbs_on])
    PabsOff .= P[stateAbs_on,stateAbs_on]
end

"""
    mo_basics!(model::StandardStoModel, parameters::Vector{Float64},kini::Float64,delta::Float64, maxrna::Int64, P::Array{Float64,2},ssp::Vector{Float64}, stateTr::Vector{Int64}, stateTr_on::Vector{Int64}, stateAbs_on::Vector{Int64}, weightsTr_off::Array{Float64,2},PabsOff::Array{Float64,2}, sspTr_off::Vecotr{Float64})

change vector and matrices used in on and off time; mean nascent rna, p_on, correlation, next burst survival
"""
function mo_basics!(model::StandardStoModel, parameters::Vector{Float64},kini::Float64,delta::Float64, maxrna::Int64, P::Array{Float64,2},ssp::Vector{Float64}, stateTr_on::Vector{Int64}, stateAbs_on::Vector{Int64}, weightsTr_off::Array{Float64,2},PabsOff::Array{Float64,2}, sspTr_off::Vector{Float64}) 
    P .= StoModel(model, parameters,kini,delta, maxrna)
    evs = eigvecs(P')
    ssp .= normalizemat(real.(evs[:,end]))
    weightsAbs_off = normalizemat(ssp[stateTr_on])
    weightsTr_off .= normalizemat(weightsAbs_off' * P[stateTr_on,stateAbs_on])
    PabsOff .= P[stateAbs_on,stateAbs_on]
    sspTr_off .= normalizemat(ssp[stateAbs_on])
end

"""
    mo_basics!(model::StandardStoModel, parameters::Vector{Float64},kini::Float64,delta::Float64, maxrna::Int64, P::Array{Float64,2},ssp::Vector{Float64}, stateTr::Vector{Int64}, stateTr_on::Vector{Int64}, stateAbs_on::Vector{Int64}, weightsTr_off::Array{Float64,2},PabsOff::Array{Float64,2}, sspTr_off::Vecotr{Float64}, Pabs::Array{Float64,2})

change vector and matrices used in ALL the statistics
"""
function mo_basics!(model::StandardStoModel, parameters::Vector{Float64},kini::Float64,delta::Float64, maxrna::Int64, P::Array{Float64,2},ssp::Vector{Float64},  stateTr_on::Vector{Int64}, stateAbs_on::Vector{Int64}, weightsTr_off::Array{Float64,2},PabsOff::Array{Float64,2}, sspTr_off::Vector{Float64}, Pabs::Array{Float64,2}) 
    P .= StoModel(model, parameters,kini,delta, maxrna)
    evs = eigvecs(P')
    ssp .= normalizemat(real.(evs[:,end]))
    weightsAbs_off = normalizemat(ssp[stateTr_on])
    weightsTr_off .= normalizemat(weightsAbs_off' * P[stateTr_on,stateAbs_on])
    PabsOff .= P[stateAbs_on,stateAbs_on]
    sspTr_off .= normalizemat(ssp[stateAbs_on])
    Pabs .= P[stateTr_on,stateTr_on]
end

"""
    mo_mnascent(ssp::Vector{Float64}, maxrna::Int64, stateTr::Vector{Int64}, nbstate::Int64)

return the mean number of nascent mrna
"""
function mo_mnascent(ssp::Vector{Float64}, maxrna::Int64, stateTr::Vector{Int64}, nbstate::Int64) 
    pB = sum(ssp[stateTr])
    prna = ssp'kron(diagm(ones(maxrna+1)),ones(nbstate))
    return [x for x in 2 : maxrna]'prna[3:end]./pB
end


"""
    mo_ontime(P::Array{Float64,2}, ssp::Vector{Float64},stateTr_on::Vector{Int64}, stateAbs_on::Vector{Int64},timevec_on::Vector{Int64})

return the on time survival probabilities
"""
function mo_ontime(P::Array{Float64,2}, ssp::Vector{Float64},stateTr_on::Vector{Int64}, stateAbs_on::Vector{Int64},timevec_on::StepRange{Int64,Int64}) 
    weightsAbs_on = ssp[stateAbs_on]./sum(ssp[stateAbs_on])     
    weightsTr_on = weightsAbs_on' * P[stateAbs_on,stateTr_on]./sum(weightsAbs_on' * P[stateAbs_on,stateTr_on])

    PabsOn = P[stateTr_on,stateTr_on]
    tempdist = weightsTr_on
    survivalspot_model_full = Vector{Float64}(undef,maximum(timevec_on))
    for i in 1:maximum(timevec_on)
        tempdist = tempdist* PabsOn
        survivalspot_model_full[i] = sum(tempdist)
    end 
    return survivalspot_model_full[timevec_on]

end

"""
    mo_offtime(P::Array{Float64,2}, ssp::Vector{Float64},stateTr_off::Vector{Int64}, stateAbs_off::Vector{Int64},timevec_off::Vector{Int64})

return the off time survival probabilities
"""
function mo_offtime(PabsOff::Array{Float64,2}, weightsTr_off::Array{Float64,2},timevec_off::StepRange{Int64,Int64}) 
    tempdist = weightsTr_off
    survivaldark_model_full = Vector{Float64}(undef,maximum(timevec_off))
    for i in 1:maximum(timevec_off)
        tempdist = tempdist* PabsOff
        survivaldark_model_full[i] = sum(tempdist)
    end 
    return survivaldark_model_full[timevec_off]
end


"""
    mo_nextbursttime(ssp::Vector{Float64},stateTr_off::Vector{Int64},timevec_nextburst::Vector{Int64})

return the time to next burst survival probabilities
"""
function mo_nextbursttime(sspTr_off::Vector{Float64},PabsOff::Array{Float64,2}, timevec_nextburst::StepRange{Int64,Int64})  
    tempdist = sspTr_off'
    survivalnextburst_model = Vector{Float64}(undef,maximum(timevec_nextburst))
    for i in 1:maximum(timevec_nextburst)
        tempdist = tempdist* PabsOff
        survivalnextburst_model[i] = sum(tempdist)
    end 
    return survivalnextburst_model
end


"""
    mo_pon(ssp::Vector{Float64},stateTr_on::Vector{Int64})

return the probability to observe a burst in steady-state
"""
function mo_pon(ssp::Vector{Float64},stateTr_on::Vector{Int64}) 
    return sum(ssp[stateTr_on])
end


"""
    mo_interburstcorr(P::Array{Float64,2}, weightsTr_off::Vector{Float64},stateAbs_on::Vector{Int64}, stateTr_on::Vector{Int64}, timehorizon::Int64)

return the correlation between two consecutive inter-burst events
"""
function mo_interburstcorr(P::Array{Float64,2}, weightsTr_off::Array{Float64,2},stateAbs_on::Vector{Int64}, stateTr_on::Vector{Int64}, timehorizon::Int64) 
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
    wpre2 = wpre*Rn./sum(wpre)
    wpre3 = wpre2*Nb*Rb./sum(wpre2)
    ET2t = wpre3*Nn*c 
    for t=1:timehorizon
        cortemp = cortemp + t*ET2t[1]*sum(wpre*Rn)
        wpre = wpre*Qn
        wpre2 = wpre*Rn./sum(wpre)
        wpre3 = wpre2*Nb*Rb./sum(wpre2)
        ET2t = wpre3*Nn*c 
        if sum(wpre)<1e-6
            break
        end
    end
    Et1 = weightsTr_off*Nn*c 
    M2T = weightsTr_off*(2*Nn-I)*Nn*c
    VarT = M2T[1] - Et1[1]^2

    return (cortemp-Et1[1]^2)/VarT
end


"""
    mo_avgintensity(P::Array{Float64,2},Pabs::Array{Float64,2}, sspTr_off::Vector{Float64},stateTr_off::Vector{Int64}, stateAbs_off::Vector{Int64},timevec_Int64ensity::Vector{Int64}, nbstate::Int64, maxrna::Int64)

return the mean track intensity, normalized to 1
"""
function mo_avgintensity(P::Array{Float64,2}, Pabs::Array{Float64,2}, sspTr_off::Vector{Float64},stateTr_off::Vector{Int64}, stateAbs_off::Vector{Int64},timevec_intensity::StepRange{Int64,Int64}, nbstate::Int64, maxrna::Int64) 
    weightsON = normalizemat(sspTr_off' * P[stateTr_off,stateAbs_off])
 
    rnanbvec_on = vcat(kron([x for x in 1:maxrna],ones(nbstate)))

    intensitytemp = weightsON
    intensity_model = Vector{Float64}(undef,length(timevec_intensity))
    for i in eachindex(timevec_intensity)
        intensity_model[i] = (intensitytemp*rnanbvec_on)[1]
        intensitytemp = intensitytemp* Pabs
    end 
    return intensity_model./maximum(intensity_model)
end


"""
    mo_rna(model::StandardStoModel, parameters::Vector{Float64},kini::Float64,delta::Float64, maxrna::Int64)

return mRNA distribution for model with parameters
"""
function mo_rna(model::StandardStoModel, parameters::Vector{Float64},kini::Float64,delta::Float64, maxrna::Int64) 
    P = StoModel(model, parameters,kini,delta, maxrna)
    evs = eigvecs(P')
    ssp = normalizemat(real.(evs[:,end]))
    
    ssd_rna = ssp'kron(diagm(ones(maxrna+1)), ones(model.nbstate))
    ssd_rna[ssd_rna .<=0] .= 1e-9 
    return ssd_rna
end

