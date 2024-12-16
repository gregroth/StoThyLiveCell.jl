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
    data_ref::DF
    data_fit::DF
    dist::DI
    model::M
    SRange::Vector{Tuple{Float32, Float32}}
    maxrnaLC::Int
    maxrnaFC::Int
    freeparametersidx::Vector{Int}
    fixedparam::Vector{Float32}
    utileMat::NamedTuple{(:stateTr, :stateTr_on, :stateAbs_on, :weightsTr_off, :P, :ssp, :PabsOff, :sspTr_Off, :Pabs, :Qrna),Tuple{Vector{Int64}, Vector{Int64}, Vector{Int64}, Vector{Float64}, Array{Float64,2}, Vector{Float64}, Array{Float64,2}, Vector{Float64}, Array{Float64,2},  Array{Float64,2}}}
    err_func::EF
end


function optim_function(SRange, FRange, optim_struct::OptimStruct, args...; maxrnaLC = 10, maxrnaFC = 60, freeparametersidx =[x for x in eachindex(SRange)], fixedparameters =[-1],  kwargs...)
    @unpack data, dist, model = optim_struct

    err_func = ini_optim(optim_struct, optim_struct.data.datagroup)


    data_fit = ini_data(optim_struct, FRange)

    #allocate memory for the utiles matrices
    (P,ssp, stateTr, stateTr_on, stateAbs_on, weightsTr_off,PabsOff, sspTr_Off, Pabs ) = StoThyLiveCell.mo_basics(model, zeros(model.nbparameters+model.nbkini+1), maxrnaLC, data.detectionLimitLC, data.detectionLimitNS) 
    Qrna = zeros(model.nbstate*(maxrnaFC+1),model.nbstate*(maxrnaFC+1))
    utileMat = (stateTr=stateTr, stateTr_on=stateTr_on, stateAbs_on=stateAbs_on, weightsTr_off=weightsTr_off, P=P, ssp=ssp, PabsOff=PabsOff, sspTr_Off=sspTr_Off, Pabs=Pabs, Qrna=Qrna)
    
    optim_struct_wrapper = OptimStructWrapper{typeof(optim_struct.data),typeof(optim_struct.dist), typeof(optim_struct.model),typeof(err_func)}(optim_struct.data, data_fit, optim_struct.dist, optim_struct.model, SRange, maxrnaLC, maxrnaFC, freeparametersidx,fixedparameters, utileMat, err_func)

    sol = start_optim(optim_struct_wrapper, args...; kwargs...)

    #bestfit parameters
    fvals = [sol[i].objective for i in eachindex(sol)] #collect all the optimization objective fct values
    minval, minidx = findmin(fvals)
    bfparameters= utiles.mergeparameter_base(fixedparameters, sol[minidx].u, freeparametersidx)
    #bestfit signal
    (mnascentmrna_model, pburst_model, survivalspot_model,survivaldark_model, survivalnextburst_model, corr_interburst, intensity_model) =    StoThyLiveCell.ModelOutput(optim_struct_wrapper.model, bfparameters, maxrnaLC, optim_struct_wrapper.data_fit.detectionLimitLC, optim_struct_wrapper.data_fit.detectionLimitNS, 10,200,200,10) 
    estimate_signal = (survival_burst = survivalspot_model, survival_interburst = survivaldark_model, survival_nextburst = survivalnextburst_model, prop_burst = pburst_model, mean_nascentrna = mnascentmrna_model, correlation_interburst = corr_interburst, intensity_burst = intensity_model)
    return sol, bfparameters, minval, minidx, estimate_signal
end


function start_optim(optim_struct_wrapper::OptimStructWrapper, args...; NbOptim::Int=1, maxtime::Int=1, maxiters::Int=1 , Method=BBO_adaptive_de_rand_1_bin_radiuslimited(), kwargs...)
    @unpack SRange, err_func = optim_struct_wrapper
    lbfull = [SRange[i][1] for i in eachindex(SRange)]
    ubfull = [SRange[i][2] for i in eachindex(SRange)]
    lb = lbfull[optim_struct_wrapper.freeparametersidx]
    ub = ubfull[optim_struct_wrapper.freeparametersidx]
    db = ub - lb
    sol = []
    for i = 1: NbOptim
        u0 = lb .+ rand(length(lb)).*db
        optprob = OptimizationFunction(err_func);
        prob = OptimizationProblem(optprob, u0, optim_struct_wrapper, lb = lb, ub = ub)
        # Import a solver package and solve the optimization problem
        push!(sol, solve(prob, Method; maxtime = maxtime, maxiters = maxiters));
    end
    return sol
end

function ini_data(optim_struct::OptimStruct, FRange; kwargs...)
    @unpack data = optim_struct
    datafit = []
    for i in eachindex(FRange)
        if (FRange[i][2]>0) & !(typeof(data.datatypes[i]) == StoThyLiveCell.Distribution_RNA)
            llimit = findfirst(data.data[i][1] .>=FRange[i][1])
            ulimit = findlast(data.data[i][1] .<=FRange[i][2])
            datatemp = (data.data[i][1][llimit:ulimit], data.data[i][2][llimit:ulimit],) 
        elseif typeof(optim_struct.dist[i]) == Distribution_RNA
            println(i)
            println(FRange[i][1])
            println(FRange[i][2])
            datatemp = data.data[i][(data.data[i] .>=FRange[i][1]) .& (data.data[i] .<=FRange[i][2])]
        else
            datatemp = data.data[i]
        end
        push!(datafit, datatemp)
    end
    return StoThyLiveCell.DataFit{typeof(data.datatypes),typeof(data.data)}(data.datatypes, data.datagroup, Tuple(datafit), data.detectionLimitLC, data.detectionLimitNS)

end

function ini_optim(optim_struct::OptimStruct, datagroup::LiveCellData; kwargs...)
    function err_func(params,optim_struct_wrapper::OptimStructWrapper)
        parameters = utiles.mergeparameter_base(optim_struct_wrapper.fixedparam, params, optim_struct_wrapper.freeparametersidx)
        #@unpack utileMat = optim_struct_wrapper
        #model outputs
        StoThyLiveCell.mo_basics!(optim_struct_wrapper.model, parameters, optim_struct_wrapper.maxrnaLC, optim_struct_wrapper.utileMat.P, optim_struct_wrapper.utileMat.ssp, optim_struct_wrapper.utileMat.stateTr_on, optim_struct_wrapper.utileMat.stateAbs_on, optim_struct_wrapper.utileMat.weightsTr_off, optim_struct_wrapper.utileMat.PabsOff) 
        error = 0.  
        for i in eachindex(optim_struct_wrapper.data_fit.datatypes)
            estimate_signal = optim_struct_wrapper.data_fit.datatypes[i](i, optim_struct_wrapper)
            error = error + optim_struct_wrapper.dist[i](estimate_signal, optim_struct_wrapper.data_fit.data[i])
        end
        return error
    end
    return err_func
end

function ini_optim(optim_struct::OptimStruct, datagroup::FixedCellData; kwargs...)
    function err_func(params,optim_struct_wrapper::OptimStructWrapper)
        parameters = utiles.mergeparameter_base(optim_struct_wrapper.fixedparam, params, optim_struct_wrapper.freeparametersidx)
        #@unpack utileMat = optim_struct_wrapper
        #model outputs
        StoThyLiveCell.distrna_basic!(optim_struct_wrapper.model, parameters, optim_struct_wrapper.maxrnaFC, optim_struct_wrapper.utileMat.Qrna) 
        estimate_signal = optim_struct_wrapper.data_fit.datatypes[1](1,optim_struct_wrapper)
        error = optim_struct_wrapper.dist[1](estimate_signal,optim_struct_wrapper.data_fit.data[1])
        return error
    end
    return err_func
end

function ini_optim(optim_struct::OptimStruct, datagroup::FixedAndLiveCellData; kwargs...)
    @warn "A mixture of live cell and fixed cell data is used in the error function"
    @warn " The  last parameter is interpreted as the degradation rate in the calculations of the mRNA number distribution"
    function err_func(params,optim_struct_wrapper::OptimStructWrapper)
        parameters = utiles.mergeparameter_base(optim_struct_wrapper.fixedparam, params, optim_struct_wrapper.freeparametersidx)
        #@unpack utileMat = optim_struct_wrapper
        #model outputs 
        StoThyLiveCell.mo_basics!(optim_struct_wrapper.model, parameters[1:end-1], optim_struct_wrapper.maxrnaLC, optim_struct_wrapper.utileMat.P, optim_struct_wrapper.utileMat.ssp, optim_struct_wrapper.utileMat.stateTr_on, optim_struct_wrapper.utileMat.stateAbs_on, optim_struct_wrapper.utileMat.weightsTr_off, optim_struct_wrapper.utileMat.PabsOff) 
        StoThyLiveCell.distrna_basic!(optim_struct_wrapper.model, vcat(parameters[1:end-2],parameters[end]), optim_struct_wrapper.maxrnaFC, optim_struct_wrapper.utileMat.Qrna) 
        error = 0.  
        for i in eachindex(optim_struct_wrapper.data_fit.datatypes)
            estimate_signal = optim_struct_wrapper.data_fit.datatypes[i](i, optim_struct_wrapper)
            error = error + optim_struct_wrapper.dist[i](estimate_signal,optim_struct_wrapper.data_fit.data[i])
        end
        return error
    end 
    return err_func
end
  


function err_func_basic(params,optim_struct_wrapper::OptimStructWrapper)
    return 0
end






function (f::Survival_Burst)(dataidx::Int, optimstruct::OptimStructWrapper)
    @unpack utileMat, data_fit = optimstruct
    
    @unpack stateTr_on, stateAbs_on, P, ssp = utileMat


    weightsAbs_on = ssp[stateAbs_on]./sum(ssp[stateAbs_on])     
    weightsTr_on = weightsAbs_on' * P[stateAbs_on,stateTr_on]./sum(weightsAbs_on' * P[stateAbs_on,stateTr_on])

    PabsOn = P[stateTr_on,stateTr_on]
    tempdist = weightsTr_on
    survivalspot_model_full = Vector{Float64}(undef,data_fit.data[dataidx][1][end])
    for i in 1:data_fit.data[dataidx][1][end]
        tempdist = tempdist* PabsOn
        survivalspot_model_full[i] = sum(tempdist)
    end 
    return survivalspot_model_full[data_fit.data[dataidx][1]]
end

function (f::Survival_InterBurst)(dataidx::Int, optimstruct::OptimStructWrapper)
    @unpack utileMat, data_fit = optimstruct
    
    @unpack weightsTr_off, PabsOff = utileMat

    tempdist = weightsTr_off
    survivaldark_model_full = Vector{Float64}(undef,data_fit.data[dataidx][1][end])
    for i in 1:data_fit.data[dataidx][1][end]
        tempdist = PabsOff'*tempdist
        survivaldark_model_full[i] = sum(tempdist)
    end 
    return survivaldark_model_full[data_fit.data[dataidx][1]]
end


function (f::Survival_NextBurst)(dataidx::Int, optimstruct::OptimStructWrapper)
    @unpack utileMat, data_fit = optimstruct
    
    @unpack sspTr_off, PabsOff = utileMat

    tempdist = sspTr_off'
    survivalnextburst_model = Vector{Float64}(undef,data_fit.data[dataidx][1][end])
    for i in 1:data_fit.data[dataidx][1][end]
        tempdist = tempdist* PabsOff
        survivalnextburst_model[i] = sum(tempdist)
    end 
    return survivalnextburst_model
end


function (f::Mean_Nascent)(dataidx::Int, optimstruct::OptimStructWrapper)
    @unpack utileMat = optimstruct
    
    @unpack stateTr, ssp = utileMat

    pB = sum(ssp[stateTr])
    prna = ssp'kron(diagm(ones(optimstruct.maxrnaLC+1)),ones(optimstruct.model.nbstate))
    return [x for x in optimstruct.data_fit.detectionLimitNS : optimstruct.maxrnaLC]'prna[optimstruct.data_fit.detectionLimitNS+1:end]./pB
end



function (f::Prob_Burst)(dataidx::Int, optimstruct::OptimStructWrapper)
    @unpack utileMat = optimstruct
    
    @unpack stateTr_on, ssp = utileMat

    return sum(ssp[stateTr_on])

end


function (f::Intensity_Burst)(dataidx::Int, optimstruct::OptimStructWrapper)
    @unpack utileMat, data_fit = optimstruct
    
    @unpack stateTr_on, stateAbs_on, Pabs, P, sspTr_off = utileMat

    weightsON = normalizemat!(P[stateAbs_on,stateTr_on]'*sspTr_off)
 
    rnanbvec_on = vcat(kron([x for x in detectionlimitLC : maxrna],ones(nbstate)))

    intensitytemp = weightsON
    intensity_model = Vector{Float64}(undef,data_fit.data[dataidx][1][end])
    for i in 1:data_fit.data[dataidx][1][end]
        intensity_model[i] = (rnanbvec_on'*intensitytemp)[1]
        intensitytemp =  Pabs'*intensitytemp
    end 
    return intensity_model./maximum(intensity_model)
end


function (f::Correlation_InterBurst)(dataidx::Int, optimstruct::OptimStructWrapper)
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

    for t=1:15000
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



function (f::Distribution_RNA)(dataidx::Int, optimstruct::OptimStructWrapper)
    @unpack utileMat, model, maxrnaFC = optimstruct
    @unpack Qrna = utileMat
    Qrna[:,end] = ones(model.nbstate* (maxrnaFC+1))
    b = zeros(model.nbstate* (maxrnaFC+1))
    b[end] = 1
    ssp = Qrna' \ b
    ssd_rna= ssp'kron(diagm(ones(maxrnaFC+1)), ones(model.nbstate))
    ssd_rna[ssd_rna .<=0] .= 1e-9 
    return ssd_rna'
end

