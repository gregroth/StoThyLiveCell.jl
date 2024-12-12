#functions used in the optimization procedure

struct OptimStruct{DF,DI,M}
    data::DF
    dist::DI
    model::M
end

function OptimStruct(data::Vector, dist::D, model::M) where {D,M}
    return OptimStruct{typeof(dist),typeof(model)}(data, dist, model)
end

struct OptimStructWrapper{DF,DI,M,EF}
    data::DF
    FRange::Vector{Tuple{Int64, Int64}}
    dist::DI
    model::M
    SRange::Vector{Tuple{Int64, Int64}}
    maxrnaLC::Int
    maxrnaFC::Int
    freeparametersidx::Vector{Int}
    fixedparam::Vector{Float32}
    utileMat::NamedTuple{(:stateTr, :stateTr_on, :stateAbs_on, :weightsTr_off, :P, :ssp, :PabsOff, :sspTr_Off, :Pabs),Tuple{Vector{Int64}, Vector{Int64}, Vector{Int64}, Vector{Float64}, Array{Float64,2}, Vector{Float64}, Array{Float64,2}, Vector{Float64}, Array{Float64,2}}}
    err_func::EF
end


function optim_function(SRange, FRange, optim_struct::OptimStruct, args...; maxrnaLC = 10, maxrnaFC = 60, freeparametersidx =[x for x in eachindex(SRange)], fixedparameters =[-1],  kwargs...)
    @unpack data, dist, model = optim_struct

    err_func = ini_optim(optim_struct)

    (P,ssp, stateTr, stateTr_on, stateAbs_on, weightsTr_off,PabsOff, sspTr_Off, Pabs ) = StoThyLiveCell.mo_basics(model, zeros(model.nbparameters+model.nbkini+1), maxrnaLC, data.detectionLimitLC, data.detectionLimitFC) 
     
    utileMat = (stateTr=stateTr, stateTr_on=stateTr_on, stateAbs_on=stateAbs_on, weightsTr_off=weightsTr_off, P=P, ssp=ssp, PabsOff=PabsOff, sspTr_Off=sspTr_Off, Pabs=Pabs)
    
    optim_struct_wrapper = OptimStructWrapper{typeof(optim_struct.data),typeof(optim_struct.dist), typeof(optim_struct.model),typeof(err_func)}(optim_struct.data, FRange, optim_struct.dist, optim_struct.model, SRange, maxrnaLC, maxrnaFC, freeparametersidx,fixedparameters, utileMat, err_func)

    sol = start_optim(optim_struct_wrapper, args...; kwargs...)

    #model, stage = optim_struct.model, optim_struct.stage
    ##params = thetax[1:end-1]
    #estimate_data = compute_distribution(params, NT, model, stage, infer_counts, filter=filter_uniform)
    #return thetax, hcat(estimate_data, reference_data)
end


function start_optim(optim_struct_wrapper::OptimStructWrapper, args...; maxtime::Int=1, maxiters::Int=1 , Method=BBO_adaptive_de_rand_1_bin_radiuslimited(), kwargs...)
    @unpack SRange, err_func = optim_struct_wrapper
    lbfull = [SRange[i][1] for i in eachindex(SRange)]
    ubfull = [SRange[i][2] for i in eachindex(SRange)]
    lb = lbfull[optim_struct_wrapper.freeparametersidx]
    ub = ubfull[optim_struct_wrapper.freeparametersidx]
    db = ub - lb
    u0 = lb .+ rand(length(lb)).*db
    optprob = OptimizationFunction(err_func);
    prob = OptimizationProblem(optprob, u0, optim_struct_wrapper, lb = lb, ub = ub)
    # Import a solver package and solve the optimization problem
    sol = solve(prob, Method; maxtime = maxtime, maxiters = maxiters);
    return sol
end

function ini_optim(optim_struct::OptimStruct; kwargs...)
    if :burst == optim_struct.data.datagroups 
        function err_func(params,optim_struct_wrapper::OptimStructWrapper)
            parameters = utiles.mergeparameter_base(optim_struct_wrapper.fixedparam, params, optim_struct_wrapper.freeparametersidx)
            #@unpack utileMat = optim_struct_wrapper
            #model outputs
            StoThyLiveCell.mo_basics!(optim_struct_wrapper.model, parameters, optim_struct_wrapper.maxrnaLC, optim_struct_wrapper.utileMat.P, optim_struct_wrapper.utileMat.ssp, optim_struct_wrapper.utileMat.stateTr_on, optim_struct_wrapper.utileMat.stateAbs_on, optim_struct_wrapper.utileMat.weightsTr_off, optim_struct_wrapper.utileMat.PabsOff) 
            error = 0.  
            for i in eachindex(optim_struct_wrapper.data.datatypes)
                estimate_signal_tot = optim_struct_wrapper.data.datatypes[i](optim_struct_wrapper.FRange[i][2],optim_struct_wrapper)
                error = error + optim_struct_wrapper.dist[i](estimate_signal_tot,optim_struct_wrapper.data.data[i], optim_struct_wrapper.FRange[i][2])
            end
            return error
        end
    end
    return err_func
end



function err_func_basic(params,optim_struct_wrapper::OptimStructWrapper)
    return 0
end






function (f::Survival_Burst)(tmax::Int, optimstruct::OptimStructWrapper)
    @unpack utileMat = optimstruct
    
    @unpack stateTr_on, stateAbs_on, P, ssp = utileMat


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

function (f::Survival_InterBurst)(tmax::Int, optimstruct::OptimStructWrapper)
    @unpack utileMat = optimstruct
    
    @unpack weightsTr_off, PabsOff = utileMat

    tempdist = weightsTr_off
    survivaldark_model_full = Vector{Float64}(undef,tmax)
    for i in 1:tmax
        tempdist = PabsOff'*tempdist
        survivaldark_model_full[i] = sum(tempdist)
    end 
    return survivaldark_model_full
end


function (f::Survival_NextBurst)(tmax::Int, optimstruct::OptimStructWrapper)
    @unpack utileMat = optimstruct
    
    @unpack sspTr_off, PabsOff = utileMat

    tempdist = sspTr_off'
    survivalnextburst_model = Vector{Float64}(undef,tmax)
    for i in 1:tmax
        tempdist = tempdist* PabsOff
        survivalnextburst_model[i] = sum(tempdist)
    end 
    return survivalnextburst_model
end


function (f::Mean_Nascent)(tmax::Int, optimstruct::OptimStructWrapper)
    @unpack utileMat = optimstruct
    
    @unpack stateTr, ssp = utileMat

    pB = sum(ssp[stateTr])
    prna = ssp'kron(diagm(ones(optimstruct.maxrnaLC+1)),ones(optimstruct.model.nbstate))
    return [x for x in optimstruct.data.detectionLimitFC : optimstruct.maxrnaLC]'prna[optimstruct.data.detectionLimitFC+1:end]./pB
end



function (f::Prob_Burst)(tmax::Int, optimstruct::OptimStructWrapper)
    @unpack utileMat = optimstruct
    
    @unpack stateTr_on, ssp = utileMat

    return sum(ssp[stateTr_on])

end


function (f::Intensity_Burst)(tmax::Int, optimstruct::OptimStructWrapper)
    @unpack utileMat = optimstruct
    
    @unpack stateTr_on, stateAbs_on, Pabs, P, sspTr_off = utileMat

    weightsON = normalizemat!(P[stateAbs_on,stateTr_on]'*sspTr_off)
 
    rnanbvec_on = vcat(kron([x for x in detectionlimitLC : maxrna],ones(nbstate)))

    intensitytemp = weightsON
    intensity_model = Vector{Float64}(undef,tmax)
    for i in 1:tmax
        intensity_model[i] = (rnanbvec_on'*intensitytemp)[1]
        intensitytemp =  Pabs'*intensitytemp
    end 
    return intensity_model./maximum(intensity_model)
end


function (f::Correlation_InterBurst)(tmax::Int, optimstruct::OptimStructWrapper)
    @unpack utileMat = optimstruct
    
    @unpack stateTr_on, stateAbs_on, weightsTr_off, P = utileMat

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

