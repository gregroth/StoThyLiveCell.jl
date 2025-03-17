"""
    normalizemat!(A::Array{Float64,2})

return the matrix normilized by its sum
"""
function normalizemat!(A::AbstractArray{T,2}) where T
    A .= A./sum(A)
end

"""
    normalizemat!(A::Vector{Float64})

return the vector normilized by its sum
"""
function normalizemat!(A::AbstractVector{T}) where T
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
    timevec_on = collect(1:1:tmaxon)
    timevec_off = collect(1:1:tmaxoff)
    timevec_nextburst = collect(1:1:tmaxnextburst)
    timevec_intensity = collect(1:1:tmaxintensity)

    (P,ssp, stateTr, stateTr_on, stateAbs_on, weightsTr_off,PabsOff, sspTr_off, Pabs) = mo_basics(model, parameters, maxrna, detectionlimitLC, detectionlimitNS) 
   
    mean_nascentrna = StoThyLiveCell.mean_nascentrna(ssp, maxrna, stateTr, model.nbstate, detectionlimitNS) 
    survival_burst = StoThyLiveCell.survival_burst(P, ssp,stateTr_on, stateAbs_on,timevec_on) 
    survival_interburst = StoThyLiveCell.survival_interburst(PabsOff, weightsTr_off, timevec_off)  
    survival_nextburst = StoThyLiveCell.survival_nextburst(sspTr_off,PabsOff, timevec_nextburst)  
    prob_burst = StoThyLiveCell.prob_burst(ssp,stateTr_on) 
    correlation_interburst = StoThyLiveCell.correlation_interburst(P, weightsTr_off,stateAbs_on, stateTr_on, 15000) 
    intensity_burst = StoThyLiveCell.intensity_burst(detectionlimitLC, P, Pabs, sspTr_off, stateTr_on, stateAbs_on,timevec_intensity, model.nbstate, maxrna)  

    return  survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst
end


function ModelOutputTest(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64, detectionlimitLC::Int64, detectionlimitNS::Int64, tmaxon::Int64,tmaxoff::Int64,tmaxnextburst::Int64,tmaxintensity::Int64) 
    timevec_on = collect(1:1:tmaxon)
    timevec_off = collect(1:1:tmaxoff)
    timevec_nextburst = collect(1:1:tmaxnextburst)
    timevec_intensity = collect(1:1:tmaxintensity)

    (P,ssp, stateTr, stateTr_on, stateAbs_on, weightsTr_off,PabsOff, sspTr_off, Pabs) = mo_basics(model, parameters, maxrna, detectionlimitLC, detectionlimitNS) 
   
    mean_nascentrna = StoThyLiveCell.mean_nascentrna(ssp, maxrna, stateTr, model.nbstate, detectionlimitNS) 
    survival_burst = StoThyLiveCell.survival_burst(P, ssp,stateTr_on, stateAbs_on,timevec_on) 
    survival_interburst = StoThyLiveCell.survival_interburst(PabsOff, weightsTr_off, timevec_off)  
    survival_nextburst = StoThyLiveCell.survival_nextburst(sspTr_off,PabsOff, timevec_nextburst)  
    prob_burst = StoThyLiveCell.prob_burst(ssp,stateTr_on) 
    correlation_interburst = StoThyLiveCell.correlation_interburstTest(P, weightsTr_off,stateAbs_on, stateTr_on, 15000) 
    intensity_burst = StoThyLiveCell.intensity_burst(detectionlimitLC, P, Pabs, sspTr_off, stateTr_on, stateAbs_on,timevec_intensity, model.nbstate, maxrna)  

    return  survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst
end

"""
    mo_basics(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64, detectionlimitLC::Int64, detectionlimitNS::Int64)

return important vectors and matrices used in the analysis of the model
"""
function mo_basics(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64, detectionlimitLC::Int64, detectionlimitNS::Int64) 
    P = StoModel(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64)
    Q = P-I
    Q[:,end] = ones(size(Q,1))
    b = zeros(size(Q,1))
    b[end] = 1
    ssp = Q' \ b

    #evs = eigvecs(P')
    #ssp = normalizemat!(real.(evs[:,end]))
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
    mo_basics(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64, detectionlimitLC::Int64, detectionlimitNS::Int64)

return important vectors and matrices used in the analysis of the model
"""
function mo_basics_nonparam(model::StandardStoModel, maxrna::Int64, detectionlimitLC::Int64, detectionlimitNS::Int64) 
    stateTr = [x for x in detectionlimitNS*model.nbstate+1 :(maxrna+1)*model.nbstate]
    stateTr_on = [x for x in detectionlimitLC*model.nbstate+1 :(maxrna+1)*model.nbstate]
    stateAbs_on = [x for x in 1 : detectionlimitLC*model.nbstate]
    return stateTr, stateTr_on, stateAbs_on
end

"""
    mo_basics!(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64, stateTr::Vector{Int64}, stateTr_on::Vector{Int64}, stateAbs_on::Vector{Int64})

change vector and matrices used in on and off time; mean nascent rna, p_on, correlation
"""
function mo_basics!(model::StandardStoModel, parameters::AbstractVector{T}, maxrna::Int64, stateTr_on::Vector{Int64}, stateAbs_on::Vector{Int64}) where T
    P = StoModelAD(model, parameters, maxrna)
#=     evs = eigvecs(P')
    ssp = normalizemat!(real.(evs[:,end])) =#
    Q = P-I
    Q[:,end] = ones(T, size(Q,1))
    b = zeros(T, size(Q,1))
    b[end] = 1
    ssp = Q' \ b
    weightsAbs_off = normalizemat!(ssp[stateTr_on])
    weightsTr_off = normalizemat!(P[stateTr_on,stateAbs_on]'*weightsAbs_off)
    PabsOff = P[stateAbs_on,stateAbs_on]
    sspTr_Off =normalizemat!(ssp[stateAbs_on])
    Pabs = P[stateTr_on,stateTr_on]
    return P, ssp, weightsTr_off, PabsOff, sspTr_Off, Pabs
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
    timevec_on = collect(1:1:tmaxon)
    timevec_off = collect(1:1:tmaxoff)
    timevec_nextburst = collect(1:1:tmaxnextburst)
    timevec_intensity = collect(1:1:tmaxintensity)

    (nascentbin, P, ssp, stateTr, stateTr_on, stateAbs_on, totnbs, Pwos, stateAbs_on_wos, statePre_on_wos, weightsAbs_off_wos, sspTr_off_wos, weightsAbs_on, sspPreB, weightsTr_on, PabsOff, weightsTr_on_wos, weightsAbsorbed_off_wos, sspwos, weightsPre_on_and_on, Rn, NR, Nc, Qn, Nn, weightsTr_off_wos, Pabs_wos, weightsON_wos, rnanbvec_on, weightsPre_on_wos) = mo_basics_wosinglet(model, parameters, maxrna, detectionLimitLC, detectionLimitNS) 
   
    mean_nascentrna = mean_nascentrna_wosinglet(ssp, nascentbin, stateTr, maxrna, detectionLimitNS) 
    survival_burst = survival_burst_wosinglet( P, stateTr_on, weightsTr_on,timevec_on) 
    survival_interburst = survival_interburst_wosinglet(PabsOff, weightsTr_on_wos,timevec_off)  
    survival_nextburst = survival_nextburst_wosinglet(weightsAbsorbed_off_wos,PabsOff, timevec_nextburst)  
    prob_burst = prob_burst_wosinglet(sspwos,weightsPre_on_and_on, stateTr_on) 
    correlation_interburst = correlation_interburst_wosinglet(Rn, NR, Nc, Qn, Nn, weightsTr_off_wos, 15000) 
    intensity_burst = intensity_burst_wosinglet(rnanbvec_on, Pwos, Pabs_wos, weightsPre_on_wos, statePre_on_wos, stateTr_on, weightsON_wos,timevec_intensity)  

    return  survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst
end


function ModelOutput_wosingletTest(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64, detectionLimitLC::Int64, detectionLimitNS::Int64, tmaxon::Int64,tmaxoff::Int64,tmaxnextburst::Int64,tmaxintensity::Int64) 
    timevec_on = collect(1:1:tmaxon)
    timevec_off = collect(1:1:tmaxoff)
    timevec_nextburst = collect(1:1:tmaxnextburst)
    timevec_intensity = collect(1:1:tmaxintensity)

    (nascentbin, P, ssp, stateTr, stateTr_on, stateAbs_on, totnbs, Pwos, stateAbs_on_wos, statePre_on_wos, weightsAbs_off_wos, sspTr_off_wos, weightsAbs_on, sspPreB, weightsTr_on, PabsOff, weightsTr_on_wos, weightsAbsorbed_off_wos, sspwos, weightsPre_on_and_on, Rn, NR, Nc, Qn, Nn, weightsTr_off_wos, Pabs_wos, weightsON_wos, rnanbvec_on, weightsPre_on_wos) = mo_basics_wosinglet(model, parameters, maxrna, detectionLimitLC, detectionLimitNS) 
   
    mean_nascentrna = mean_nascentrna_wosinglet(ssp, nascentbin, stateTr, maxrna, detectionLimitNS) 
    survival_burst = survival_burst_wosinglet( P, stateTr_on, weightsTr_on,timevec_on) 
    survival_interburst = survival_interburst_wosinglet(PabsOff, weightsTr_on_wos,timevec_off)  
    survival_nextburst = survival_nextburst_wosinglet(weightsAbsorbed_off_wos,PabsOff, timevec_nextburst)  
    prob_burst = prob_burst_wosinglet(sspwos,weightsPre_on_and_on, stateTr_on) 
    correlation_interburst = correlation_interburstTest_wosinglet(Pwos, stateAbs_on_wos, stateTr_on, weightsTr_off_wos, 15000) 
    intensity_burst = intensity_burst_wosinglet(rnanbvec_on, Pwos, Pabs_wos, weightsPre_on_wos, statePre_on_wos, stateTr_on, weightsON_wos,timevec_intensity)  

    return  survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst
end

"""
    mo_basics_wosinglet(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64, detectionlimitLC::Int64, detectionlimitNS::Int64)

return important vectors and matrices used in the analysis of the model without burst singlets
"""
function mo_basics_wosinglet(model::StandardStoModel, parameters::AbstractVector{T}, maxrna::Int64, detectionlimitLC::Int64, detectionlimitNS::Int64) where T
    P = StoModel(model::StandardStoModel, parameters, maxrna)
    Q = P-I
    Q[:,end] = ones(T, size(Q,1))
    b = zeros(T, size(Q,1))
    b[end] = 1
    ssp = Q' \ b
    #evs = eigvecs(P')
    #ssp = normalizemat!(real.(evs[:,end]))
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
    Qwos = Pwos-I
    Qwos[:,end] = ones(size(Qwos,1))
    bwos = zeros(size(Qwos,1))
    bwos[end] = 1
    sspwos = Qwos' \ bwos
    #evswos = eigvecs(Pwos')
    #sspwos = normalizemat!(real.(evswos[:,end]))
    
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
    weightsON_wos = normalizemat!(Pwos[statePre_on_wos,stateTr_on]'*weightsPre_on_wos)
    Pabs_wos = Pwos[stateTr_on,stateTr_on]

    return nascentbin, P, ssp, stateTr, stateTr_on, stateAbs_on, totnbs, Pwos, stateAbs_on_wos, statePre_on_wos, weightsAbs_off_wos, sspTr_off_wos, weightsAbs_on, sspPreB, weightsTr_on, PabsOff, weightsTr_on_wos, weightsAbsorbed_off_wos, sspwos, weightsPre_on_and_on, Rn, NR, Nc, Qn, Nn, weightsTr_off_wos, Pabs_wos, weightsON_wos, rnanbvec_on, weightsPre_on_wos
end

"""
    mo_basics_wosinglet(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64, detectionlimitLC::Int64, detectionlimitNS::Int64)

return important vectors and matrices used in the analysis of the model without burst singlets
"""
function mo_basics_wosinglet_nonparam(model::StandardStoModel, maxrna::Int64, detectionlimitLC::Int64, detectionlimitNS::Int64) 
    stateTr = [x for x in detectionlimitNS*model.nbstate+1 :(maxrna+1)*model.nbstate]
    stateTr_on = [x for x in detectionlimitLC*model.nbstate+1 :(maxrna+1)*model.nbstate]
    stateAbs_on = [x for x in 1 : detectionlimitLC*model.nbstate]
    nascentbin = kron(diagm(ones(maxrna+1)),ones(model.nbstate))
    totnbs = (maxrna+1)*model.nbstate
    stateAbs_on_wos = vcat(stateAbs_on,[x for x in totnbs+1 : totnbs+length(stateTr_on)])
    statePre_on_wos = [x for x in totnbs+1 : totnbs+length(stateTr_on)]
    rnanbvec_on = vcat(kron([x for x in detectionlimitLC : maxrna],ones(model.nbstate)))
   
    return nascentbin, stateTr, stateTr_on, stateAbs_on, totnbs, stateAbs_on_wos, statePre_on_wos,  rnanbvec_on
end



"""
    mo_basics_wosinglet(model::StandardStoModel, parameters::Vector{Float64}, maxrna::Int64, detectionlimitLC::Int64, detectionlimitNS::Int64)

return important vectors and matrices used in the analysis of the model without burst singlets
"""
function mo_basics_wosingletAD(model::StandardStoModel, parameters::AbstractVector{T}, maxrna::Int64, stateTr_on::Vector{Int64}, stateAbs_on::Vector{Int64}, totnbs::Int,stateAbs_on_wos::Vector{Int64}, statePre_on_wos::Vector{Int64}  ) where T
    P = StoModelAD(model, parameters, maxrna)
    Q = P-I
    Q[:,end] = ones(T, size(Q,1))
    b = zeros(T, size(Q,1))
    b[end] = 1
    ssp = Q' \ b
#calculation of transition matrix adding the singlet states
    Pwos = zeros(T, totnbs+length(stateTr_on),totnbs+length(stateTr_on))
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
 #survival inter-burst
    Qwos = Pwos-I
    Qwos[:,end] = ones(T, size(Qwos,1))
    bwos = zeros(T, size(Qwos,1))
    bwos[end] = 1
    sspwos = Qwos' \ bwos

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
    weightsPre_on_wos = normalizemat!(sspwos[statePre_on_wos])
    weightsON_wos = normalizemat!(Pwos[statePre_on_wos,stateTr_on]'*weightsPre_on_wos)
    Pabs_wos = Pwos[stateTr_on,stateTr_on]

    return  P, ssp, Pwos,  weightsAbs_off_wos, sspTr_off_wos, weightsAbs_on, sspPreB, weightsTr_on, PabsOff, weightsTr_on_wos, weightsAbsorbed_off_wos, sspwos, weightsPre_on_and_on, Rn, NR, Nc, Qn, Nn, weightsTr_off_wos, Pabs_wos, weightsON_wos, weightsPre_on_wos
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
    Q = P-I
    Q[:,end] = ones(size(Q,1))
    b = zeros(size(Q,1))
    b[end] = 1
    ssp .= Q' \ b
    #evs = eigvecs(P')
    #ssp .= normalizemat!(real.(evs[:,end]))
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
    Q = P-I
    Q[:,end] = ones(size(Q,1))
    b = zeros(size(Q,1))
    b[end] = 1
    ssp .= Q' \ b
    #evs = eigvecs(P')
    #ssp .= normalizemat!(real.(evs[:,end]))

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
    Qwos = Pwos-I
    Qwos[:,end] = ones(size(Qwos,1))
    bwos = zeros(size(Qwos,1))
    bwos[end] = 1
    sspwos .= Qwos' \ bwos
    
    #evswos = eigvecs(Pwos')
    #sspwos .= normalizemat!(real.(evswos[:,end]))
    
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
    weightsON_wos .= normalizemat!(Pwos[statePre_on_wos,stateTr_on]'*weightsPre_on_wos)
    Pabs_wos .= Pwos[stateTr_on,stateTr_on]
end


"""
mean_nascentrna(ssp::AbstractVector{T}, maxrna::Int64, stateTr::Vector{Int64}, nbstate::Int64, detectionlimitNS::Int)

return the mean number of nascent mrna
"""
function mean_nascentrna(ssp::AbstractVector{T}, maxrna::Int64, stateTr::Vector{Int64}, nbstate::Int64, detectionlimitNS::Int) where T
    pB = sum(ssp[stateTr])
    prna = ssp'kron(diagm(ones(T,maxrna+1)),ones(T,nbstate))
    return [x for x in detectionlimitNS : maxrna]'prna[detectionlimitNS+1:end]./pB
end


"""
    survival_burst( P::AbstractArray{T,2}, ssp::AbstractVector{T},stateTr_on::Vector{Int64}, stateAbs_on::Vector{Int64},timevec_on::Vector{Int64}) 
    
return the on time survival probabilities
"""
function survival_burst( P::AbstractArray{T,2}, ssp::AbstractVector{T},stateTr_on::Vector{Int64}, stateAbs_on::Vector{Int64},timevec_on::Vector{Int64}) where T
    weightsAbs_on = ssp[stateAbs_on]./sum(ssp[stateAbs_on])     
    weightsTr_on = weightsAbs_on' * P[stateAbs_on,stateTr_on]./sum(weightsAbs_on' * P[stateAbs_on,stateTr_on])

    PabsOn = P[stateTr_on,stateTr_on]
    tempdist = weightsTr_on
    survivalspot_model_full = Vector{T}(undef,timevec_on[end])
    for i in 1:timevec_on[end]
        tempdist = tempdist* PabsOn
        survivalspot_model_full[i] = sum(tempdist)
    end 
    return survivalspot_model_full[timevec_on]

end

"""
survival_interburst(PabsOff::AbstractArray{T,2}, weightsTr_off::AbstractVector{T},timevec_off::Vector{Int64}) 

return the off time survival probabilities
"""
function survival_interburst(PabsOff::AbstractArray{T,2}, weightsTr_off::AbstractVector{T},timevec_off::Vector{Int64}) where T
    tempdist = weightsTr_off
    survivaldark_model_full = Vector{T}(undef,timevec_off[end])
    for i in 1:timevec_off[end]
        tempdist = PabsOff'*tempdist
        survivaldark_model_full[i] = sum(tempdist)
    end 
    return survivaldark_model_full[timevec_off]
end


"""
    survival_nextburst(sspTr_off::AbstractVector{T},PabsOff::AbstractArray{T,2}, timevec_nextburst::Vector{Int64}) 

return the time to next burst survival probabilities
"""
function survival_nextburst(sspTr_off::AbstractVector{T},PabsOff::AbstractArray{T,2}, timevec_nextburst::Vector{Int64})  where T
    tempdist = sspTr_off'
    survivalnextburst_model = Vector{T}(undef,timevec_nextburst[end])
    for i in 1:timevec_nextburst[end]
        tempdist = tempdist* PabsOff
        survivalnextburst_model[i] = sum(tempdist)
    end 
    return survivalnextburst_model[timevec_nextburst]
end


"""
prob_burst(ssp::AbstractVector{T},stateTr_on::Vector{Int64})

return the probability to observe a burst in steady-state
"""
function prob_burst(ssp::AbstractVector{T},stateTr_on::Vector{Int64}) where T
    return sum(ssp[stateTr_on])
end


"""
correlation_interburst( P::AbstractArray{T,2}, weightsTr_off::AbstractVector{T},stateAbs_on::Vector{Int64}, stateTr_on::Vector{Int64}, timehorizon::Int64

return the correlation between two consecutive inter-burst events
"""
function correlation_interburst( P::AbstractArray{T,2}, weightsTr_off::AbstractVector{T},stateAbs_on::Vector{Int64}, stateTr_on::Vector{Int64}, timehorizon::Int64) where T
    #correlation of the interburst durations
    Qn = P[stateAbs_on,stateAbs_on]
    Rn = P[stateAbs_on,stateTr_on]

    Qb = P[stateTr_on,stateTr_on]
    Rb = P[stateTr_on,stateAbs_on]
    c = ones(T,length(stateAbs_on))

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

function correlation_interburstTest( P::AbstractArray{T,2}, weightsTr_off::AbstractVector{T},stateAbs_on::Vector{Int64}, stateTr_on::Vector{Int64}, timehorizon::Int64) where T
    #correlation of the interburst durations
    PNN = P[stateAbs_on,stateAbs_on]
    PNB = P[stateAbs_on,stateTr_on]

    PBB = P[stateTr_on,stateTr_on]
    PBN = P[stateTr_on,stateAbs_on]
    c = ones(T,length(stateAbs_on))

    NN = (I - PNN)^(-1)
    NB = (I - PBB)^(-1)

    Nc = NN*c
    
    cortemp=0
    pivec = weightsTr_off'
    M = PNB*NB*PBN*NN*c
    corvec = M
    for t=2:timehorizon
        M = PNN*M
        corvec= corvec + t*M
        if sum(M)<1e-9
            break
        end
    end
    #cortemp = sum(pivec .* M)
    cortemp = pivec*corvec
    Et1 = Nc'weightsTr_off 
    M2T = Nc'*(2*NN'-I)*weightsTr_off
    VarT = M2T[1] - Et1[1]^2

    return (cortemp-Et1[1]^2)/VarT
end

"""
intensity_burst(detectionlimitLC::Int64, P::AbstractArray{T,2}, Pabs::AbstractArray{T,2}, sspTr_off::AbstractVector{T},stateTr_off::Vector{Int64}, stateAbs_off::Vector{Int64},timevec_intensity::Vector{Int64}St, nbstate::Int64, maxrna::Int64) 

return the mean track intensity, normalized to 1
"""
function intensity_burst(detectionlimitLC::Int64, P::AbstractArray{T,2}, Pabs::AbstractArray{T,2}, sspTr_off::AbstractVector{T},stateTr_on::Vector{Int64}, stateAbs_on::Vector{Int64},timevec_intensity::Vector{Int64}, nbstate::Int64, maxrna::Int64) where T
    weightsON = normalizemat!(P[stateAbs_on,stateTr_on]'*sspTr_off)
 
    rnanbvec_on = vcat(kron([x for x in detectionlimitLC : maxrna],ones(nbstate)))

    intensitytemp = weightsON
    intensity_model = Vector{T}(undef,length(timevec_intensity))
    for i in eachindex(timevec_intensity)
        intensity_model[i] = (rnanbvec_on'*intensitytemp)[1]
        intensitytemp =  Pabs'*intensitytemp
    end 
    return intensity_model./maximum(intensity_model)
end


"""
    mo_rna(model::StandardStoModel, parameters::AbstractVector{T}, maxrna::Int64)

return mRNA distribution for model with parameters
"""
function distribution_mrna(model::StandardStoModel, parameters::AbstractVector{T}, maxrna::Int64) where T
    Q = StoModel_RateMat(model, parameters, maxrna)
    Q[:,end] = ones(model.nbstate* (maxrna+1))
    b = zeros(model.nbstate* (maxrna+1))
    b[end] = 1
    ssp = Q' \ b
    ssd_rna= ssp'kron(diagm(ones(maxrna+1)), ones(model.nbstate))
    ssd_rna[ssd_rna .<=0] .= 1e-9 
    return ssd_rna'
#=     evs = eigvecs(P')
    ssp = normalizemat!(real.(evs[:,end]))
    
    ssd_rna = ssp'kron(diagm(ones(maxrna+1)), ones(model.nbstate))
    ssd_rna[ssd_rna .<=0] .= 1e-9 
    return ssd_rna =#
end

"""
    mo_rna(model::StandardStoModel, parameters::AbstractVector{T}, maxrna::Int64)

return mRNA distribution for model with parameters
"""
function distribution_mrna(Q::AbstractArray{T,2}, maxrna::Int64, nbstate::Int64) where T
    Q[:,end] = ones(size(Q,1))
    b = zeros(size(Q,1))
    b[end] = 1
    ssp = Q' \ b
    ssd_rna= ssp'kron(diagm(ones(maxrna+1)), ones(nbstate))
    ssd_rna[ssd_rna .<=0] .= 1e-9 
    return ssd_rna'
#=     evs = eigvecs(P')
    ssp = normalizemat!(real.(evs[:,end]))
    
    ssd_rna = ssp'kron(diagm(ones(maxrna+1)), ones(model.nbstate))
    ssd_rna[ssd_rna .<=0] .= 1e-9 
    return ssd_rna =#
end

function distrna_basic!(model::StandardStoModel, parameters::AbstractVector{T}, maxrna::Int64, Q::AbstractArray{T,2}) where T
    Q .= StoModel_RateMat(model, parameters, maxrna)
end

function distrna_basic(model::StandardStoModel, parameters::AbstractVector{T}, maxrna::Int64) where T
    Q = StoModel_RateMat(model, parameters, maxrna)
end


#-----------------------
#wo singlets
#--------



function mean_nascentrna_wosinglet(ssp::AbstractVector{T}, nascentbin::Array{Float64,2}, stateTr::Vector{Int64}, maxrna::Int, detectionLimitNS::Int) where T
    pB = sum(ssp[stateTr])
    prna = ssp'nascentbin
    return [x for x in detectionLimitNS : maxrna]'prna[detectionLimitNS+1:end]./pB
end


"""
survival_burst_wosinglet( P::AbstractArray{T,2}, ssp::AbstractVector{T},stateTr_on::Vector{Int64}, stateAbs_on::Vector{Int64},timevec_on::Vector{Int64})

return the on time survival probabilities, without burst singlets
"""
function survival_burst_wosinglet( P::AbstractArray{T,2}, stateTr_on::Vector{Int64}, weightsTr_on::AbstractVector{T},timevec_on::Vector{Int64}) where T
    PabsOn = P[stateTr_on,stateTr_on]'
    survivalspot_model_full = Vector{T}(undef,timevec_on[end])
    for i in 1:timevec_on[end]
        weightsTr_on = PabsOn*weightsTr_on
        survivalspot_model_full[i] = sum(weightsTr_on)
    end 
    return survivalspot_model_full[timevec_on]./survivalspot_model_full[1]
end

"""
    survival_interburst_wosinglet(PabsOff::AbstractArray{T,2}, weightsTr_off_wos::AbstractVector{T},timevec_off::StepRange{Int64,Int64})

return the off time survival probabilities, without burst singlets
"""
function survival_interburst_wosinglet(PabsOff::AbstractArray{T,2}, weightsTr_on_wos::AbstractVector{T},timevec_off::Vector{Int64}) where T
    survivaldark_model_full = Vector{T}(undef,timevec_off[end])
    for i in 1:timevec_off[end]
        weightsTr_on_wos = PabsOff'*weightsTr_on_wos
        survivaldark_model_full[i] = sum(weightsTr_on_wos)
    end 
    return survivaldark_model_full[timevec_off]
end


"""
    survival_nextburst_wosingelts(ssp::AbstractVector{T},stateTr_off::Vector{Int64},timevec_nextburst::Vector{Int64})

return the time to next burst survival probabilities, without burst singlets
"""
function survival_nextburst_wosinglet(weightsAbsorbed_off_wos::AbstractVector{T},PabsOff::AbstractArray{T,2}, timevec_nextburst::Vector{Int64})  where T
    survivalnextburst_model = Vector{T}(undef,timevec_nextburst[end])
    for i in 1:timevec_nextburst[end]
        weightsAbsorbed_off_wos = PabsOff'*weightsAbsorbed_off_wos
        survivalnextburst_model[i] = sum(weightsAbsorbed_off_wos)
    end 
    return survivalnextburst_model[timevec_nextburst]
end


"""
prob_burst_wosinglet(sspwos::AbstractVector{T},weightsPre_on_and_on::AbstractVector{T}, stateTr_on::Vector{Int64})_wosinglet(ssp::AbstractVector{T},stateTr_on::Vector{Int64})

return the probability to observe a burst in steady-state
"""
function prob_burst_wosinglet(sspwos::AbstractVector{T},weightsPre_on_and_on::AbstractVector{T}, stateTr_on::Vector{Int64}) where T
    return sum(sspwos[stateTr_on]) + sum(weightsPre_on_and_on)
end


"""
correlation_interburst_wosinglet( P::AbstractArray{T,2}, weightsTr_off::AbstractVector{T},stateAbs_on::Vector{Int64}, stateTr_on::Vector{Int64}, timehorizon::Int64)

return the correlation between two consecutive inter-burst events, without burst singlets
"""
function correlation_interburst_wosinglet(Rn::AbstractArray{T,2}, NR::AbstractArray{T,2}, Nc::AbstractVector{T}, Qn::AbstractArray{T,2}, Nn::AbstractArray{T,2}, weightsTr_off_wos::AbstractVector{T}, timehorizon::Int64) where T
    #correlation of the interburst durations
    cortemp=0
    wpre = Qn'*weightsTr_off_wos
    
    
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

function correlation_interburstTest_wosinglet( Pwos::AbstractArray{T,2},stateAbs_on_wos::Vector{Int64},stateTr_on_wos::Vector{Int64}, weightsTr_off_wos::AbstractVector{T}, timehorizon::Int64) where T
    #correlation of the interburst durations
    PNN = Pwos[stateAbs_on_wos,stateAbs_on_wos]
    PNB = Pwos[stateAbs_on_wos,stateTr_on_wos]

    PBB = Pwos[stateTr_on_wos,stateTr_on_wos]
    PBN = Pwos[stateTr_on_wos,stateAbs_on_wos]
    c = ones(T,length(stateAbs_on_wos))

    NN = (I - PNN)^(-1)
    NB = (I - PBB)^(-1)

    Nc = NN*c
    
    cortemp=0
    pivec = weightsTr_off_wos'
    M = PNB*NB*PBN*NN*c
    corvec = M
    for t=2:timehorizon
        M = PNN*M
        corvec= corvec + t*M
        if sum(M)<1e-9
            break
        end
    end
    #cortemp = sum(pivec .* M)
    cortemp = pivec*corvec
    Et1 = Nc'weightsTr_off_wos
    M2T = Nc'*(2*NN'-I)*weightsTr_off_wos
    VarT = M2T[1] - (Et1[1])^2

    return (cortemp-Et1[1]^2)/VarT
end

"""
intensity_burst_wosinglet(detectionlimitLC::Int64, P::AbstractArray{T,2},Pabs::AbstractArray{T,2}, sspTr_off::AbstractVector{T},stateTr_off::Vector{Int64}, stateAbs_off::Vector{Int64},timevec_Int64ensity::Vector{Int64}, nbstate::Int64, maxrna::Int64)

return the mean track intensity, normalized to 1, without burst singlets
"""
function intensity_burst_wosinglet(rnanbvec_on::Vector{Float64}, Pwos::AbstractArray{T,2}, Pabs_wos::AbstractArray{T,2}, weightsPre_on_wos::AbstractVector{T}, statePre_on_wos::Vector{Int64}, stateTr_on::Vector{Int64}, weightsON_wos::AbstractVector{T},timevec_intensity::Vector{Int64})   where T
    intensity_model = Vector{T}(undef,timevec_intensity[end])
    intensity_model[1] = (rnanbvec_on'*(weightsPre_on_wos.*sum(Pwos[statePre_on_wos,stateTr_on], dims=2))./sum((weightsPre_on_wos.*sum(Pwos[statePre_on_wos,stateTr_on], dims=2))))[1]
    for i in 2:timevec_intensity[end]
        intensity_model[i] = (rnanbvec_on'*weightsON_wos)[1]
        weightsON_wos =  Pabs_wos'*weightsON_wos
    end 
    return intensity_model[timevec_intensity]./maximum(intensity_model)
end
