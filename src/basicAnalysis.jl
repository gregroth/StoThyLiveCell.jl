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
    ModelOutput(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64, detectionlimitLC::Int64, detectionlimitNS::Int64,tmaxon::Int64,tmaxoff::Int64,tmaxnextburst::Int64,tmaxintensity::Int64)

Return the main outputs of the model: mean nb of nascent mrna, probability to observe a burst, ON time survival,
OFF time survival,next burst time survival, correlation of consecutive inter burst evts, average intensity track
#Arguments
- `model::StandardStoModel:` a model topology
- `parameters::Vector{Float64}`: list of the rates of the transition between the 'promoter' states
- `maxrna::Int`: maximum number of mRNA considered in the model
- `detectionlimitLC::Int64`: minimum number of mRNA detectable in live cell experiment
- `detectionlimitNS::Int64`: minimum number of nascent mRNA detectable in FISH experiment 
- `tmaxon::Int64`: maximum time for the ON time suvival probabilities 
- `tmaxoff::Int64`: maximum time for the OFF time suvival probabilities 
- `tmaxnextburst::Int64`: maximum time for the Next burst time suvival probabilities 
- `tmaxintensity::Int64`: maximum time for the averaged track intensity 
"""
function ModelOutput(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64, detectionlimitLC::Int64, detectionlimitNS::Int64, tmaxon::Int64,tmaxoff::Int64,tmaxnextburst::Int64,tmaxintensity::Int64) 
    timevec_on = 1:1:tmaxon
    timevec_off = 1:1:tmaxoff
    timevec_nextburst = 1:1:tmaxnextburst
    timevec_intensity = 1:1:tmaxintensity
    #model instance
    P = StoModel(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64)
    
    #nascent mrna
    evs = eigvecs(P')
    ssp = real.(evs[:,end]./sum(evs[:,end]))

    stateBns = [x for x in detectionlimitNS*model.nbstate+1 :(maxrna+1)*model.nbstate]

    pB = sum(ssp[stateBns])
    prna = ssp'kron(diagm(ones(maxrna+1)), ones(model.nbstate))
    mnascentmrna_model = [x for x in detectionlimitNS : maxrna]'prna[detectionlimitNS+1:end]./pB
    
    #on times
    stateTr_on = [x for x in detectionlimitLC*model.nbstate+1 :(maxrna+1)*model.nbstate]
    stateAbs_on = [x for x in 1 : detectionlimitLC*model.nbstate]
    
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
    stateAbs_off = [x for x in detectionlimitLC*model.nbstate+1 :(maxrna+1)*model.nbstate]
    stateTr_off = [x for x in 1 : detectionlimitLC*model.nbstate]
    
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
 
    rnanbvec_on = vcat(kron([x for x in detectionlimitLC:maxrna],ones(model.nbstate)))

    Pabs = Qb

    intensitytemp = weightsON
    intensity_model = Vector{Float64}(undef,length(timevec_intensity))
    for i in eachindex(timevec_intensity)
        intensity_model[i] = (intensitytemp*rnanbvec_on)[1]
        intensitytemp = intensitytemp* Pabs
    end 
    intensity_model .= intensity_model./maximum(intensity_model) 
    return  survivalspot_model_full,survivaldark_model_full, survivalnextburst_model,  pburst_model, mnascentmrna_model, corr_interburst, intensity_model
end


"""
    mo_basics(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64, detectionlimitLC::Int64, detectionlimitNS::Int64)

return important vectors and matrices used in the analysis of the model
"""
function mo_basics(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64, detectionlimitLC::Int64, detectionlimitNS::Int64) 
    P = StoModel(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64)
    evs = eigvecs(P')
    ssp = normalizemat!(real.(evs[:,end]))
    stateTr = [x for x in detectionlimitNS*model.nbstate+1 :(maxrna+1)*model.nbstate]
    stateTr_on = [x for x in detectionlimitLC*model.nbstate+1 :(maxrna+1)*model.nbstate]
    stateAbs_on = [x for x in 1 : detectionlimitLC*model.nbstate]
    weightsAbs_off = normalizemat!(ssp[stateTr_on])
    weightsTr_off = normalizemat!(P[stateTr_on,stateAbs_on]'*weightsAbs_off)
    PabsOff = P[stateAbs_on,stateAbs_on]
    sspTr_off =normalizemat!(ssp[stateAbs_on])
    Pabs = P[stateTr_on,stateTr_on]
    return P,ssp, stateTr, stateTr_on, stateAbs_on, weightsTr_off,PabsOff, sspTr_off, Pabs
end


"""
    ModelOutput_wosinglet(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64, detectionlimitLC::Int64, detectionlimitNS::Int64,tmaxon::Int64,tmaxoff::Int64,tmaxnextburst::Int64,tmaxintensity::Int64)

Return the main outputs of the model: mean nb of nascent mrna, probability to observe a burst, ON time survival,
OFF time survival,next burst time survival, correlation of consecutive inter burst evts, average intensity track
#Arguments
- `model::StandardStoModel:` a model topology
- `parameters::Vector{Float64}`: list of the rates of the transition between the 'promoter' states
- `maxrna::Int`: maximum number of mRNA considered in the model
- `detectionlimitLC::Int64`: minimum number of mRNA detectable in live cell experiment
- `detectionlimitNS::Int64`: minimum number of nascent mRNA detectable in FISH experiment 
- `tmaxon::Int64`: maximum time for the ON time suvival probabilities 
- `tmaxoff::Int64`: maximum time for the OFF time suvival probabilities 
- `tmaxnextburst::Int64`: maximum time for the Next burst time suvival probabilities 
- `tmaxintensity::Int64`: maximum time for the averaged track intensity 
"""
function ModelOutput_wosinglet(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64, detectionLimitLC::Int64, detectionLimitNS::Int64, tmaxon::Int64,tmaxoff::Int64,tmaxnextburst::Int64,tmaxintensity::Int64) 
    timevec_on = 1:1:tmaxon
    timevec_off = 1:1:tmaxoff
    timevec_nextburst = 1:1:tmaxnextburst
    timevec_intensity = 1:1:tmaxintensity

    (nascentbin, P, ssp, stateTr, stateTr_on, stateAbs_on, totnbs, Pwos, stateAbs_on_wos, statePre_on_wos, weightsAbs_off_wos, sspTr_off_wos, weightsAbs_on, sspPreB, weightsTr_on, PabsOff, weightsTr_on_wos, weightsAbsorbed_off_wos, sspwos, weightsPre_on_and_on, Rn, NR, Nc, Qn, Nn, weightsTr_off_wos, Pabs_wos, weightsON_wos, rnanbvec_on) = mo_basics_wosinglet(model, parameters, maxrna, detectionLimitLC, detectionLimitNS) 
   
    mean_nascentrna = mo_mnascent_wosinglet(ssp, nascentbin, stateTr, maxrna, detectionLimitNS) 
    survival_burst = mo_ontime_wosinglet( P, stateTr_on, weightsTr_on,timevec_on) 
    survival_interburst = mo_offtime_wosinglet(PabsOff, weightsTr_on_wos,timevec_off)  
    survival_nextburst = mo_nextbursttime_wosinglet(weightsAbsorbed_off_wos,PabsOff, timevec_nextburst)  
    prob_burst = mo_pon_wosinglet(sspwos,weightsPre_on_and_on, stateTr_on) 
    correlation_interburst = mo_interburstcorr_wosinglet(Rn, NR, Nc, Qn, Nn, weightsTr_off_wos, 15000) 
    intensity_burst = mo_avgintensity_wosinglet(rnanbvec_on, Pabs_wos, weightsON_wos,timevec_intensity)  

    return  survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst
end



"""
    mo_basics_wosinglet(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64, detectionlimitLC::Int64, detectionlimitNS::Int64)

return important vectors and matrices used in the analysis of the model without burst singlets
"""
function mo_basics_wosinglet(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64, detectionlimitLC::Int64, detectionlimitNS::Int64) 
    P = StoModel(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64)
    evs = eigvecs(P')
    ssp = normalizemat!(real.(evs[:,end]))
    stateTr = [x for x in detectionlimitNS*model.nbstate+1 :(maxrna+1)*model.nbstate]
    stateTr_on = [x for x in detectionlimitLC*model.nbstate+1 :(maxrna+1)*model.nbstate]
    stateAbs_on = [x for x in 1 : detectionlimitLC*model.nbstate]
    nascentbin = kron(diagm(ones(maxrna+1)),ones(model.nbstate))
#calculation of transition matrix adding the singlet states
    totnbs = (maxrna+1)*model.nbstate
    Pwos = zeros(totnbs+length(stateTr_on),totnbs+length(stateTr_on))
    for i in axes(P,1)
        if i in stateTr_on
            for j in axes(P)
                Pwos[i,j] = P[i,j]
            end
        else
            for (j,s) in enumerate(stateTr_on)
                Pwos[i,j+totnbs] = P[i,s]
            end
            for (j,s) in enumerate(stateAbs_on)
                Pwos[i,s] = P[i,s]
            end
        end
    end
    for (i,s) in enumerate(stateTr_on)
        for j in axes(Pwos)
            Pwos[i+totnbs,j] = Pwos[s,j]
        end
    end

#survival burst
    weightsAbs_on = ssp[stateAbs_on]./sum(ssp[stateAbs_on])     
    weightsTr_on = normalizemat!(P[stateAbs_on,stateTr_on]'*weightsAbs_on)

#=     weightsAbs_off = normalizemat!(ssp[stateTr_on])
    weightsTr_off = normalizemat!(P[stateTr_on,stateAbs_on]'*weightsAbs_off)
 =#

 #survival inter-burst
    evswos = eigvecs(Pwos')
    sspwos = normalizemat!(real.(evswos[:,end]))
    
    stateAbs_on_wos = vcat(stateAbs_on,[x for x in totnbs+1 : totnbs+length(stateTr_on)])
    statePre_on_wos = [x for x in totnbs+1 : totnbs+length(stateTr_on)]

    weightsAbs_off_wos = normalizemat!(sspwos[stateTr_on])
    weightsTr_off_wos =  normalizemat!(Pwos[stateTr_on,stateAbs_on_wos]'*weightsAbs_off_wos)


    PabsOff =  Pwos[stateAbs_on_wos,stateAbs_on_wos]

    weightsTr_on_wos = PabsOff'*weightsTr_off_wos

#survival next burst

    sspTr_off_wos =normalizemat!(sspwos[stateAbs_on_wos])
    weightsAbsorbed_off_wos = normalizemat!(PabsOff'*sspTr_off_wos)

#probability burst

    sspPreB = sspwos[statePre_on_wos]
    weightsPre_on_and_on = Pwos[statePre_on_wos,stateTr_on]'*sspPreB


    
#correlation
    Qn = Pwos[stateAbs_on_wos,stateAbs_on_wos]
    Rn = Pwos[stateAbs_on_wos,stateTr_on]

    Qb = Pwos[stateTr_on,stateTr_on]
    Rb = Pwos[stateTr_on,stateAbs_on_wos]
    c = ones(length(stateAbs_on_wos))

    Nn = (I - Qn)^(-1)
    Nb = (I - Qb)^(-1)

    NR = Nb*Rb
    Nc = Nn*c

#intensity 

    rnanbvec_on = vcat(kron([x for x in detectionlimitLC : maxrna],ones(model.nbstate)))
    weightsPre_on_wos = normalizemat!(sspwos[statePre_on_wos])
    weightsON_wos = Pwos[statePre_on_wos,stateTr_on]'*weightsPre_on_wos
    Pabs_wos = Pwos[stateTr_on,stateTr_on]

    return nascentbin, P, ssp, stateTr, stateTr_on, stateAbs_on, totnbs, Pwos, stateAbs_on_wos, statePre_on_wos, weightsAbs_off_wos, sspTr_off_wos, weightsAbs_on, sspPreB, weightsTr_on, PabsOff, weightsTr_on_wos, weightsAbsorbed_off_wos, sspwos, weightsPre_on_and_on, Rn, NR, Nc, Qn, Nn, weightsTr_off_wos, Pabs_wos, weightsON_wos, rnanbvec_on
end


"""
    mo_basics!(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64,  P::Array{Float64,2},ssp::Vector{Float64}, stateTr::Vector{Int64}, stateTr_on::Vector{Int64}, stateAbs_on::Vector{Int64}, weightsTr_off::Vector{Float64},PabsOff::Array{Float64,2})

change vector and matrices used in on and off time; mean nascent rna, p_on, correlation
"""
function mo_basics!(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64, P::Array{Float64,2},ssp::Vector{Float64}, stateTr_on::Vector{Int64}, stateAbs_on::Vector{Int64}, weightsTr_off::Vector{Float64},PabsOff::Array{Float64,2}) 
    P .= StoModel(model, parameters, maxrna)
    evs = eigvecs(P')
    ssp .= normalizemat!(real.(evs[:,end]))
    weightsAbs_off = normalizemat!(ssp[stateTr_on])
    weightsTr_off .= normalizemat!(P[stateTr_on,stateAbs_on]'*weightsAbs_off)
    PabsOff .= P[stateAbs_on,stateAbs_on]
end

"""
    mo_basics!(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64,  P::Array{Float64,2},ssp::Vector{Float64}, stateTr::Vector{Int64}, stateTr_on::Vector{Int64}, stateAbs_on::Vector{Int64}, weightsTr_off::Vector{Float64},PabsOff::Array{Float64,2}, sspTr_off::Vecotr{Float64})

change vector and matrices used in on and off time; mean nascent rna, p_on, correlation, next burst survival
"""
function mo_basics!(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64,  P::Array{Float64,2},ssp::Vector{Float64}, stateTr_on::Vector{Int64}, stateAbs_on::Vector{Int64}, weightsTr_off::Vector{Float64},PabsOff::Array{Float64,2}, sspTr_off::Vector{Float64}) 
    P .= StoModel(model, parameters, maxrna)
    evs = eigvecs(P')
    ssp .= normalizemat!(real.(evs[:,end]))
    weightsAbs_off = normalizemat!(ssp[stateTr_on])
    weightsTr_off .= normalizemat!(P[stateTr_on,stateAbs_on]'*weightsAbs_off)
    PabsOff .= P[stateAbs_on,stateAbs_on]
    sspTr_off .= normalizemat!(ssp[stateAbs_on])
end

"""
    mo_basics!(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64, P::Array{Float64,2},ssp::Vector{Float64}, stateTr::Vector{Int64}, stateTr_on::Vector{Int64}, stateAbs_on::Vector{Int64}, weightsTr_off::Vector{Float64},PabsOff::Array{Float64,2}, sspTr_off::Vecotr{Float64}, Pabs::Array{Float64,2})

change vector and matrices used in ALL the statistics
"""
function mo_basics!(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64, P::Array{Float64,2},ssp::Vector{Float64},  stateTr_on::Vector{Int64}, stateAbs_on::Vector{Int64}, weightsTr_off::Vector{Float64},PabsOff::Array{Float64,2}, sspTr_off::Vector{Float64}, Pabs::Array{Float64,2}) 
    P .= StoModel(model, parameters, maxrna)
    evs = eigvecs(P')
    ssp .= normalizemat!(real.(evs[:,end]))
    weightsAbs_off = normalizemat!(ssp[stateTr_on])
    weightsTr_off .= normalizemat!(P[stateTr_on,stateAbs_on]'*weightsAbs_off)
    PabsOff .= P[stateAbs_on,stateAbs_on]
    sspTr_off .= normalizemat!(ssp[stateAbs_on])
    Pabs .= P[stateTr_on,stateTr_on]
end

"""
    mo_basics_wosinglet!(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64, P::Array{Float64,2},ssp::Vector{Float64},  stateTr_on::Vector{Int64}, stateAbs_on::Vector{Int64}, weightsTr_off::Vector{Float64},PabsOff::Array{Float64,2}, sspTr_off::Vector{Float64}, Pabs::Array{Float64,2})

change vector and matrices used in ALL the statistics, model without burst singlets
"""
function mo_basics_wosinglet!(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64, P::Array{Float64,2}, ssp::Vector{Float64}, stateTr_on::Vector{Int64}, stateAbs_on::Vector{Int64}, totnbs::Int, Pwos::Array{Float64,2}, stateAbs_on_wos::Vector{Int64}, statePre_on_wos::Vector{Int64}, weightsAbs_off_wos::Vector{Float64}, sspTr_off_wos::Vector{Float64}, weightsAbs_on::Vector{Float64}, sspPreB::Vector{Float64}, weightsTr_on::Vector{Float64}, PabsOff::Array{Float64,2}, weightsTr_on_wos::Vector{Float64}, weightsAbsorbed_off_wos::Vector{Float64}, sspwos::Vector{Float64}, weightsPre_on_and_on::Vector{Float64}, Rn::Array{Float64,2}, NR::Array{Float64,2}, Nc::Vector{Float64}, Qn::Array{Float64,2}, Nn::Array{Float64,2}, weightsTr_off_wos::Vector{Float64}, Pabs_wos::Array{Float64,2}, weightsON_wos::Vector{Float64}) 
    P .= StoModel(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64)
    evs = eigvecs(P')
    ssp .= normalizemat!(real.(evs[:,end]))
#calculation of transition matrix adding the singlet states
    for i in axes(P,1)
        if i in stateTr_on
            for j in axes(P)
                Pwos[i,j] = P[i,j]
            end
        else
            for (j,s) in enumerate(stateTr_on)
                Pwos[i,j+totnbs] = P[i,s]
            end
            for (j,s) in enumerate(stateAbs_on)
                Pwos[i,s] = P[i,s]
            end
        end
    end
    for (i,s) in enumerate(stateTr_on)
        for j in axes(Pwos)
            Pwos[i+totnbs,j] = Pwos[s,j]
        end
    end

#survival burst
    weightsAbs_on .= ssp[stateAbs_on]./sum(ssp[stateAbs_on])     
    weightsTr_on .= normalizemat!(P[stateAbs_on,stateTr_on]'*weightsAbs_on)

#=     weightsAbs_off = normalizemat!(ssp[stateTr_on])
    weightsTr_off = normalizemat!(P[stateTr_on,stateAbs_on]'*weightsAbs_off)
 =#

 #survival inter-burst
    evswos = eigvecs(Pwos')
    sspwos .= normalizemat!(real.(evswos[:,end]))
    
    weightsAbs_off_wos .= normalizemat!(sspwos[stateTr_on])
    weightsTr_off_wos .=  normalizemat!(Pwos[stateTr_on,stateAbs_on_wos]'*weightsAbs_off_wos)


    PabsOff .=  Pwos[stateAbs_on_wos,stateAbs_on_wos]

    weightsTr_on_wos .= PabsOff'*weightsTr_off_wos

#survival next burst

    sspTr_off_wos .=normalizemat!(sspwos[stateAbs_on_wos])
    weightsAbsorbed_off_wos .= normalizemat!(PabsOff'*sspTr_off_wos)

#probability burst

    sspPreB .= sspwos[statePre_on_wos]
    weightsPre_on_and_on .= Pwos[statePre_on_wos,stateTr_on]'*sspPreB


    
#correlation
    Qn .= Pwos[stateAbs_on_wos,stateAbs_on_wos]
    Rn .= Pwos[stateAbs_on_wos,stateTr_on]

    Qb = Pwos[stateTr_on,stateTr_on]
    Rb = Pwos[stateTr_on,stateAbs_on_wos]
    c = ones(length(stateAbs_on_wos))

    Nn .= (I - Qn)^(-1)
    Nb = (I - Qb)^(-1)

    NR .= Nb*Rb
    Nc .= Nn*c

#intensity 
    weightsPre_on_wos = normalizemat!(sspwos[statePre_on_wos])
    weightsON_wos .= Pwos[statePre_on_wos,stateTr_on]'*weightsPre_on_wos
    Pabs_wos .= Pwos[stateTr_on,stateTr_on]
end

"""
    mo_mnascent(ssp::Vector{Float64}, maxrna::Int64, stateTr::Vector{Int64}, nbstate::Int64)

return the mean number of nascent mrna
"""
function mo_mnascent(ssp::Vector{Float64}, maxrna::Int64, stateTr::Vector{Int64}, nbstate::Int64, detectionlimitNS::Int) 
    pB = sum(ssp[stateTr])
    prna = ssp'kron(diagm(ones(maxrna+1)),ones(nbstate))
    return [x for x in detectionlimitNS : maxrna]'prna[detectionlimitNS+1:end]./pB
end


"""
    mo_ontime( P::Array{Float64,2}, ssp::Vector{Float64},stateTr_on::Vector{Int64}, stateAbs_on::Vector{Int64},timevec_on::Vector{Int64})

return the on time survival probabilities
"""
function mo_ontime( P::Array{Float64,2}, ssp::Vector{Float64},stateTr_on::Vector{Int64}, stateAbs_on::Vector{Int64},timevec_on::StepRange{Int64,Int64}) 
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
    mo_offtime( P::Array{Float64,2}, ssp::Vector{Float64},stateTr_off::Vector{Int64}, stateAbs_off::Vector{Int64},timevec_off::Vector{Int64})

return the off time survival probabilities
"""
function mo_offtime(PabsOff::Array{Float64,2}, weightsTr_off::Vector{Float64},timevec_off::StepRange{Int64,Int64}) 
    tempdist = weightsTr_off
    survivaldark_model_full = Vector{Float64}(undef,maximum(timevec_off))
    for i in 1:maximum(timevec_off)
        tempdist = PabsOff'*tempdist
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
    mo_interburstcorr( P::Array{Float64,2}, weightsTr_off::Vector{Float64},stateAbs_on::Vector{Int64}, stateTr_on::Vector{Int64}, timehorizon::Int64)

return the correlation between two consecutive inter-burst events
"""
function mo_interburstcorr( P::Array{Float64,2}, weightsTr_off::Vector{Float64},stateAbs_on::Vector{Int64}, stateTr_on::Vector{Int64}, timehorizon::Int64) 
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

    for t=1:timehorizon
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


"""
    mo_avgintensity(detectionlimitLC::Int64, P::Array{Float64,2},Pabs::Array{Float64,2}, sspTr_off::Vector{Float64},stateTr_off::Vector{Int64}, stateAbs_off::Vector{Int64},timevec_Int64ensity::Vector{Int64}, nbstate::Int64, maxrna::Int64)

return the mean track intensity, normalized to 1
"""
function mo_avgintensity(detectionlimitLC::Int64, P::Array{Float64,2}, Pabs::Array{Float64,2}, sspTr_off::Vector{Float64},stateTr_off::Vector{Int64}, stateAbs_off::Vector{Int64},timevec_intensity::StepRange{Int64,Int64}, nbstate::Int64, maxrna::Int64) 
    weightsON = normalizemat!(P[stateTr_off,stateAbs_off]'*sspTr_off)
 
    rnanbvec_on = vcat(kron([x for x in detectionlimitLC : maxrna],ones(nbstate)))

    intensitytemp = weightsON
    intensity_model = Vector{Float64}(undef,length(timevec_intensity))
    for i in eachindex(timevec_intensity)
        intensity_model[i] = (rnanbvec_on'*intensitytemp)[1]
        intensitytemp =  Pabs'*intensitytemp
    end 
    return intensity_model./maximum(intensity_model)
end


"""
    mo_rna(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64)

return mRNA distribution for model with parameters
"""
function mo_rna(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64) 
    Q = StoModel_RateMat(model, parameters, maxrna)
    Q[:,end] = ones(model.nbstate* (maxrna+1))
    b = zeros(model.nbstate* (maxrna+1))
    b[end] = 1
    ssp = Q' \ b
    ssd_rna= ssp'kron(diagm(ones(maxrna+1)), ones(model.nbstate))
    ssd_rna[ssd_rna .<=0] .= 1e-9 
    return ssd_rna
#=     evs = eigvecs(P')
    ssp = normalizemat!(real.(evs[:,end]))
    
    ssd_rna = ssp'kron(diagm(ones(maxrna+1)), ones(model.nbstate))
    ssd_rna[ssd_rna .<=0] .= 1e-9 
    return ssd_rna =#
end


function distrna_basic!(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64, Q::Array{Float64,2}) 
    Q .= StoModel_RateMat(model, parameters, maxrna)
end




#-----------------------
#wo singlets
#--------



function mo_mnascent_wosinglet(ssp::Vector{Float64}, nascentbin::Array{Float64,2}, stateTr::Vector{Int64}, maxrna::Int, detectionLimitNS::Int) 
    pB = sum(ssp[stateTr])
    prna = ssp'nascentbin
    return [x for x in detectionLimitNS : maxrna]'prna[detectionLimitNS+1:end]./pB
end


"""
    mo_ontime_wosinglet( P::Array{Float64,2}, ssp::Vector{Float64},stateTr_on::Vector{Int64}, stateAbs_on::Vector{Int64},timevec_on::Vector{Int64})

return the on time survival probabilities, without burst singlets
"""
function mo_ontime_wosinglet( P::Array{Float64,2}, stateTr_on::Vector{Int64}, weightsTr_on::Vector{Float64},timevec_on::StepRange{Int64,Int64}) 
    PabsOn = P[stateTr_on,stateTr_on]'
    survivalspot_model_full = Vector{Float64}(undef,maximum(timevec_on))
    for i in 1:maximum(timevec_on)
        weightsTr_on = PabsOn*weightsTr_on
        survivalspot_model_full[i] = sum(weightsTr_on)
    end 
    return survivalspot_model_full[timevec_on]./survivalspot_model_full[1]
end

"""
    mo_offtime_wosinglet(PabsOff::Array{Float64,2}, weightsTr_off_wos::Vector{Float64},timevec_off::StepRange{Int64,Int64})

return the off time survival probabilities, without burst singlets
"""
function mo_offtime_wosinglet(PabsOff::Array{Float64,2}, weightsTr_on_wos::Vector{Float64},timevec_off::StepRange{Int64,Int64}) 
    survivaldark_model_full = Vector{Float64}(undef,maximum(timevec_off))
    for i in 1:maximum(timevec_off)
        weightsTr_on_wos = PabsOff'*weightsTr_on_wos
        survivaldark_model_full[i] = sum(weightsTr_on_wos)
    end 
    return survivaldark_model_full[timevec_off]
end


"""
    mo_nextbursttime_wosingelts(ssp::Vector{Float64},stateTr_off::Vector{Int64},timevec_nextburst::Vector{Int64})

return the time to next burst survival probabilities, without burst singlets
"""
function mo_nextbursttime_wosinglet(weightsAbsorbed_off_wos::Vector{Float64},PabsOff::Array{Float64,2}, timevec_nextburst::StepRange{Int64,Int64})  
    survivalnextburst_model = Vector{Float64}(undef,maximum(timevec_nextburst))
    for i in 1:maximum(timevec_nextburst)
        weightsAbsorbed_off_wos = PabsOff'*weightsAbsorbed_off_wos
        survivalnextburst_model[i] = sum(weightsAbsorbed_off_wos)
    end 
    return survivalnextburst_model
end


"""
    mo_pon_wosinglet(ssp::Vector{Float64},stateTr_on::Vector{Int64})

return the probability to observe a burst in steady-state
"""
function mo_pon_wosinglet(sspwos::Vector{Float64},weightsPre_on_and_on::Vector{Float64}, stateTr_on::Vector{Int64}) 
    return sum(sspwos[stateTr_on]) + sum(weightsPre_on_and_on)
end


"""
    mo_interburstcorr_wosinglet( P::Array{Float64,2}, weightsTr_off::Vector{Float64},stateAbs_on::Vector{Int64}, stateTr_on::Vector{Int64}, timehorizon::Int64)

return the correlation between two consecutive inter-burst events, without burst singlets
"""
function mo_interburstcorr_wosinglet(Rn::Array{Float64,2}, NR::Array{Float64,2}, Nc::Vector{Float64}, Qn::Array{Float64,2}, Nn::Array{Float64,2}, weightsTr_off_wos::Vector{Float64}, timehorizon::Int64) 
    #correlation of the interburst durations
    cortemp=0
    wpre = weightsTr_off_wos
    
    
    wpre2 = Rn'*wpre./sum(wpre)
    wpre3 = NR'*wpre2/sum(wpre2)
    ET2t = Nc'*wpre3 .-1

    for t=1:timehorizon
        cortemp = cortemp + t*ET2t[1]*sum(Rn'*wpre)
        wpre = Qn'*wpre
        wpre2 = Rn'*wpre./sum(wpre)
        wpre3 = NR'*wpre2./sum(wpre2)
        ET2t = Nc'*wpre3 .-1
        if sum(wpre)<1e-6
            break
        end
    end
    Et1 = Nc'weightsTr_off_wos .-1
    M2T = Nc'*(2*Nn'-3I)*weightsTr_off_wos .+1
    VarT = M2T[1] - Et1[1]^2

    return (cortemp-Et1[1]^2)/VarT
end


"""
    mo_avgintensity_wosinglet(detectionlimitLC::Int64, P::Array{Float64,2},Pabs::Array{Float64,2}, sspTr_off::Vector{Float64},stateTr_off::Vector{Int64}, stateAbs_off::Vector{Int64},timevec_Int64ensity::Vector{Int64}, nbstate::Int64, maxrna::Int64)

return the mean track intensity, normalized to 1, without burst singlets
"""
function mo_avgintensity_wosinglet(rnanbvec_on::Vector{Float64}, Pabs_wos::Array{Float64,2}, weightsON_wos::Vector{Float64},timevec_intensity::StepRange{Int64,Int64})  
    intensity_model = Vector{Float64}(undef,length(timevec_intensity))
    for i in eachindex(timevec_intensity)
        intensity_model[i] = (rnanbvec_on'*weightsON_wos)[1]
        weightsON_wos =  Pabs_wos'*weightsON_wos
    end 
    return intensity_model./maximum(intensity_model)
end