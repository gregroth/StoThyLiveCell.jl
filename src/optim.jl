#functions used in the optimization procedure

struct OptimStruct{DF,DI,M}
    data::DF
    dist::DI
    model::M
end

function OptimStruct(data::Vector, dist::D, model::M) where {D,M}
    return OptimStruct{typeof(dist),typeof(model)}(data, dist, model)
end

struct OptimStructWrapper{DF,DI,M,EF,UM}
    data_ref::DF
    data_fit::DF
    dist::DI
    model::M
    SRange::Vector{Tuple{Float32, Float32}}
    maxrnaLC::Int
    maxrnaFC::Int
    freeparametersidx::Vector{Int}
    fixedparam::Vector{Float32}
    utileMat::UM #NamedTuple{(:stateTr, :stateTr_on, :stateAbs_on, :weightsTr_off, :P, :ssp, :PabsOff, :sspTr_Off, :Pabs, :Qrna),Tuple{Vector{Int64}, Vector{Int64}, Vector{Int64}, Vector{Float64}, Array{Float64,2}, Vector{Float64}, Array{Float64,2}, Vector{Float64}, Array{Float64,2},  Array{Float64,2}}}
    err_func::EF
end


function optim_function(SRange, FRange, optim_struct::OptimStruct, args...; maxrnaLC = 10, maxrnaFC = 60, freeparametersidx =[x for x in eachindex(SRange)], fixedparameters =[-1],  kwargs...)
    @unpack data, dist, model = optim_struct

    if data.burstsinglet == :with
        err_func = ini_optim(optim_struct, optim_struct.data.datagroup)
    elseif data.burstsinglet == :without
        err_func = ini_optim_withoutsinglet(optim_struct, optim_struct.data.datagroup)
    end

    data_fit = ini_data(optim_struct, FRange)

    #allocate memory for the utiles matrices
    if data.burstsinglet == :with
        (P,ssp, stateTr, stateTr_on, stateAbs_on, weightsTr_off,PabsOff, sspTr_Off, Pabs ) = StoThyLiveCell.mo_basics(model, ones(model.nbparameters+model.nbkini+1), maxrnaLC, data.detectionLimitLC, data.detectionLimitNS) 
        Qrna = zeros(model.nbstate*(maxrnaFC+1),model.nbstate*(maxrnaFC+1))
        utileMat = (stateTr=stateTr, stateTr_on=stateTr_on, stateAbs_on=stateAbs_on, weightsTr_off=weightsTr_off, P=P, ssp=ssp, PabsOff=PabsOff, sspTr_Off=sspTr_Off, Pabs=Pabs, Qrna=Qrna)
        optim_struct_wrapper = OptimStructWrapper{typeof(optim_struct.data),typeof(optim_struct.dist), typeof(optim_struct.model),typeof(err_func), typeof(utileMat)}(optim_struct.data, data_fit, optim_struct.dist, optim_struct.model, SRange, maxrnaLC, maxrnaFC, freeparametersidx,fixedparameters, utileMat, err_func)
    elseif data.burstsinglet == :without
        (nascentbin, P, ssp, stateTr, stateTr_on, stateAbs_on, totnbs, Pwos, stateAbs_on_wos, statePre_on_wos, weightsAbs_off_wos, sspTr_off_wos, weightsAbs_on, sspPreB, weightsTr_on, PabsOff, weightsTr_on_wos, weightsAbsorbed_off_wos, sspwos, weightsPre_on_and_on, Rn, NR, Nc, Qn, Nn, weightsTr_off_wos, Pabs_wos, weightsON_wos, rnanbvec_on, weightsPre_on_wos) = StoThyLiveCell.mo_basics_wosinglet(model, ones(model.nbparameters+model.nbkini+1), maxrnaLC, data.detectionLimitLC, data.detectionLimitNS) 
        Qrna = zeros(model.nbstate*(maxrnaFC+1),model.nbstate*(maxrnaFC+1))
        utileMat = (nascentbin=nascentbin, P=P, ssp=ssp, stateTr=stateTr, stateTr_on=stateTr_on, stateAbs_on=stateAbs_on, totnbs=totnbs, Pwos=Pwos, stateAbs_on_wos=stateAbs_on_wos, statePre_on_wos=statePre_on_wos, weightsAbs_off_wos=weightsAbs_off_wos, sspTr_off_wos=sspTr_off_wos, weightsAbs_on=weightsAbs_on, sspPreB=sspPreB, weightsTr_on=weightsTr_on, PabsOff=PabsOff, weightsTr_on_wos=weightsTr_on_wos, weightsAbsorbed_off_wos=weightsAbsorbed_off_wos, sspwos=sspwos, weightsPre_on_and_on=weightsPre_on_and_on, Rn=Rn, NR=NR, Nc=Nc, Qn=Qn, Nn=Nn, weightsTr_off_wos=weightsTr_off_wos, Pabs_wos=Pabs_wos, weightsON_wos=weightsON_wos, rnanbvec_on=rnanbvec_on, weightsPre_on_wos=weightsPre_on_wos, Qrna=Qrna)
        optim_struct_wrapper = OptimStructWrapper{typeof(optim_struct.data),typeof(optim_struct.dist), typeof(optim_struct.model),typeof(err_func), typeof(utileMat)}(optim_struct.data, data_fit, optim_struct.dist, optim_struct.model, SRange, maxrnaLC, maxrnaFC, freeparametersidx,fixedparameters, utileMat, err_func)
    end

    #run the optimization
    sol = start_optim(optim_struct_wrapper, args...; kwargs...)

    #bestfit parameters
    fvals = [sol[i].objective for i in eachindex(sol)] #collect all the optimization objective fct values
    minval, minidx = findmin(fvals)
    bfparameters= utiles.mergeparameter_base(fixedparameters, sol[minidx].u, freeparametersidx)
    #bestfit signal
    if data.burstsinglet == :with
        if data.datagroup == LiveCellData()
            (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =    StoThyLiveCell.ModelOutput(optim_struct_wrapper.model, bfparameters, maxrnaLC, optim_struct_wrapper.data_fit.detectionLimitLC, optim_struct_wrapper.data_fit.detectionLimitNS, 10,200,200,10) 
            bfparameters_mrna = vcat(bfparameters[1:end-1],.01)
            mrna_distribution_model = StoThyLiveCell.distribution_mrna(optim_struct_wrapper.model, bfparameters_mrna, optim_struct_wrapper.maxrnaFC) 
        elseif data.datagroup == FixedCellData()
            bfparameters_lc = vcat(bfparameters[1:end-1],1.)
            (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =    StoThyLiveCell.ModelOutput(optim_struct_wrapper.model,bfparameters_lc, maxrnaLC, optim_struct_wrapper.data_fit.detectionLimitLC, optim_struct_wrapper.data_fit.detectionLimitNS, 10,200,200,10) 
            mrna_distribution_model = StoThyLiveCell.distribution_mrna(optim_struct_wrapper.model, bfparameters, optim_struct_wrapper.maxrnaFC) 
        elseif data.datagroup == FixedAndLiveCellData()
            (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =    StoThyLiveCell.ModelOutput(optim_struct_wrapper.model, bfparameters[1:end-1], maxrnaLC, optim_struct_wrapper.data_fit.detectionLimitLC, optim_struct_wrapper.data_fit.detectionLimitNS, 10,200,200,10) 
            bfparameters_mrna = vcat(bfparameters[1:end-1],bfparameters[end])
            mrna_distribution_model = StoThyLiveCell.distribution_mrna(optim_struct_wrapper.model, bfparameters_mrna, optim_struct_wrapper.maxrnaFC) 
        end
    elseif data.burstsinglet == :without
        if data.datagroup == LiveCellData()
            (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =    StoThyLiveCell.ModelOutput_wosinglet(optim_struct_wrapper.model, bfparameters, maxrnaLC, optim_struct_wrapper.data_fit.detectionLimitLC, optim_struct_wrapper.data_fit.detectionLimitNS, 10,200,200,10) 
            bfparameters_mrna = vcat(bfparameters[1:end-1],.01)
            mrna_distribution_model = StoThyLiveCell.distribution_mrna(optim_struct_wrapper.model, bfparameters_mrna, optim_struct_wrapper.maxrnaFC) 
        elseif data.datagroup == FixedCellData()
            bfparameters_lc = vcat(bfparameters[1:end-1],1.)
            (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =    StoThyLiveCell.ModelOutput(optim_struct_wrapper.model,bfparameters_lc, maxrnaLC, optim_struct_wrapper.data_fit.detectionLimitLC, optim_struct_wrapper.data_fit.detectionLimitNS, 10,200,200,10) 
            mrna_distribution_model = StoThyLiveCell.distribution_mrna(optim_struct_wrapper.model, bfparameters, optim_struct_wrapper.maxrnaFC) 
        elseif data.datagroup == FixedAndLiveCellData()
            (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =    StoThyLiveCell.ModelOutput_wosinglet(optim_struct_wrapper.model, bfparameters[1:end-1], maxrnaLC, optim_struct_wrapper.data_fit.detectionLimitLC, optim_struct_wrapper.data_fit.detectionLimitNS, 10,200,200,10) 
            bfparameters_mrna = vcat(bfparameters[1:end-1],bfparameters[end])
            mrna_distribution_model = StoThyLiveCell.distribution_mrna(optim_struct_wrapper.model, bfparameters_mrna, optim_struct_wrapper.maxrnaFC) 
        end
    end
    estimate_signal = (survival_burst = survival_burst, survival_interburst = survival_interburst, survival_nextburst = survival_nextburst, prob_burst = prob_burst, mean_nascentrna = mean_nascentrna, correlation_interburst = correlation_interburst, intensity_burst = intensity_burst, mrna_distribution = mrna_distribution_model)    
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
        if (FRange[i][2]>0) & !(data.datatypes[i] == StoThyLiveCell.Distribution_RNA())
            llimit = findfirst(data.data[i][1] .>=FRange[i][1])
            ulimit = findlast(data.data[i][1] .<=FRange[i][2])
            datatemp = (data.data[i][1][llimit:ulimit], data.data[i][2][llimit:ulimit],) 
        elseif data.datatypes[i] == StoThyLiveCell.Distribution_RNA()
            datatemp = data.data[i][(data.data[i] .>=FRange[i][1]) .& (data.data[i] .<=FRange[i][2])]
        else
            datatemp = data.data[i]
        end
        push!(datafit, datatemp)
    end
    return StoThyLiveCell.DataFit{typeof(data.datatypes),typeof(data.data)}(data.datatypes, data.datagroup, Tuple(datafit), data.detectionLimitLC, data.detectionLimitNS, data.burstsinglet)

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
  
function ini_optim_withoutsinglet(optim_struct::OptimStruct, datagroup::LiveCellData; kwargs...)
    function err_func(params,optim_struct_wrapper::OptimStructWrapper)
        parameters = utiles.mergeparameter_base(optim_struct_wrapper.fixedparam, params, optim_struct_wrapper.freeparametersidx)
        #@unpack utileMat = optim_struct_wrapper
        #model outputs
        StoThyLiveCell.mo_basics_wosinglet!(optim_struct_wrapper.model, parameters, optim_struct_wrapper.maxrnaLC, optim_struct_wrapper.utileMat.P, optim_struct_wrapper.utileMat.ssp, optim_struct_wrapper.utileMat.stateTr_on, optim_struct_wrapper.utileMat.stateAbs_on, optim_struct_wrapper.utileMat.totnbs, optim_struct_wrapper.utileMat.Pwos, optim_struct_wrapper.utileMat.stateAbs_on_wos, optim_struct_wrapper.utileMat.statePre_on_wos, optim_struct_wrapper.utileMat.weightsAbs_off_wos, optim_struct_wrapper.utileMat.sspTr_off_wos, optim_struct_wrapper.utileMat.weightsAbs_on, optim_struct_wrapper.utileMat.sspPreB, optim_struct_wrapper.utileMat.weightsTr_on, optim_struct_wrapper.utileMat.PabsOff, optim_struct_wrapper.utileMat.weightsTr_on_wos, optim_struct_wrapper.utileMat.weightsAbsorbed_off_wos, optim_struct_wrapper.utileMat.sspwos, optim_struct_wrapper.utileMat.weightsPre_on_and_on, optim_struct_wrapper.utileMat.Rn, optim_struct_wrapper.utileMat.NR, optim_struct_wrapper.utileMat.Nc, optim_struct_wrapper.utileMat.Qn, optim_struct_wrapper.utileMat.Nn, optim_struct_wrapper.utileMat.weightsTr_off_wos, optim_struct_wrapper.utileMat.Pabs_wos, optim_struct_wrapper.utileMat.weightsON_wos) 
        error = 0.  
        for i in eachindex(optim_struct_wrapper.data_fit.datatypes)
            estimate_signal = optim_struct_wrapper.data_fit.datatypes[i](i, optim_struct_wrapper, :without)
            error = error + optim_struct_wrapper.dist[i](estimate_signal, optim_struct_wrapper.data_fit.data[i])
        end
        return error
    end
    return err_func
end

function ini_optim_withoutsinglet(optim_struct::OptimStruct, datagroup::FixedAndLiveCellData; kwargs...)
    @warn "A mixture of live cell and fixed cell data is used in the error function"
    @warn " The  last parameter is interpreted as the degradation rate in the calculations of the mRNA number distribution"
    function err_func(params,optim_struct_wrapper::OptimStructWrapper)
        parameters = utiles.mergeparameter_base(optim_struct_wrapper.fixedparam, params, optim_struct_wrapper.freeparametersidx)
        #@unpack utileMat = optim_struct_wrapper
        #model outputs 
        StoThyLiveCell.mo_basics_wosinglet!(optim_struct_wrapper.model, parameters, optim_struct_wrapper.maxrnaLC, optim_struct_wrapper.utileMat.P, optim_struct_wrapper.utileMat.ssp, optim_struct_wrapper.utileMat.stateTr_on, optim_struct_wrapper.utileMat.stateAbs_on, optim_struct_wrapper.utileMat.totnbs, optim_struct_wrapper.utileMat.Pwos, optim_struct_wrapper.utileMat.stateAbs_on_wos, optim_struct_wrapper.utileMat.statePre_on_wos, optim_struct_wrapper.utileMat.weightsAbs_off_wos, optim_struct_wrapper.utileMat.sspTr_off_wos, optim_struct_wrapper.utileMat.weightsAbs_on, optim_struct_wrapper.utileMat.sspPreB, optim_struct_wrapper.utileMat.weightsTr_on, optim_struct_wrapper.utileMat.PabsOff, optim_struct_wrapper.utileMat.weightsTr_on_wos, optim_struct_wrapper.utileMat.weightsAbsorbed_off_wos, optim_struct_wrapper.utileMat.sspwos, optim_struct_wrapper.utileMat.weightsPre_on_and_on, optim_struct_wrapper.utileMat.Rn, optim_struct_wrapper.utileMat.NR, optim_struct_wrapper.utileMat.Nc, optim_struct_wrapper.utileMat.Qn, optim_struct_wrapper.utileMat.Nn, optim_struct_wrapper.utileMat.weightsTr_off_wos, optim_struct_wrapper.utileMat.Pabs_wos, optim_struct_wrapper.utileMat.weightsON_wos) 
        StoThyLiveCell.distrna_basic!(optim_struct_wrapper.model, vcat(parameters[1:end-2],parameters[end]), optim_struct_wrapper.maxrnaFC, optim_struct_wrapper.utileMat.Qrna) 
        error = 0.  
        for i in eachindex(optim_struct_wrapper.data_fit.datatypes)
            estimate_signal = optim_struct_wrapper.data_fit.datatypes[i](i, optim_struct_wrapper, :without)
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

    return survival_burst(P, ssp, stateTr_on, stateAbs_on, data_fit.data[dataidx][1]) 
   
end

function (f::Survival_InterBurst)(dataidx::Int, optimstruct::OptimStructWrapper)
    @unpack utileMat, data_fit = optimstruct
    @unpack weightsTr_off, PabsOff = utileMat
    return survival_interburst(PabsOff, weightsTr_off,data_fit.data[dataidx][1])
end


function (f::Survival_NextBurst)(dataidx::Int, optimstruct::OptimStructWrapper)
    @unpack utileMat, data_fit = optimstruct
    @unpack sspTr_off, PabsOff = utileMat

    return survival_nextburst(sspTr_off::Vector{Float64},PabsOff,data_fit.data[dataidx][1])
end


function (f::Mean_Nascent)(dataidx::Int, optimstruct::OptimStructWrapper)
    @unpack utileMat, data_fit = optimstruct
    @unpack stateTr, ssp = utileMat

    return mean_nascentrna(ssp, optimstruct.maxrnaLC, stateTr, optimstruct.model.nbstate, data_fit.detectionLimitNS)
end



function (f::Prob_Burst)(dataidx::Int, optimstruct::OptimStructWrapper)
    @unpack utileMat = optimstruct
    @unpack stateTr_on, ssp = utileMat

    return prob_burst(ssp,stateTr_on) 

end


function (f::Intensity_Burst)(dataidx::Int, optimstruct::OptimStructWrapper)
    @unpack utileMat, data_fit, model, maxrnaLC = optimstruct
    @unpack stateTr_on, stateAbs_on, Pabs, P, sspTr_off = utileMat

    return intensity_burst(data_fit.detectionLimitLC, P, Pabs, sspTr_off,stateTr_on, stateAbs_on,data_fit.data[dataidx][1], model.nbstate, maxrnaLC)
end


function (f::Correlation_InterBurst)(dataidx::Int, optimstruct::OptimStructWrapper)
    @unpack utileMat = optimstruct
    
    @unpack stateTr_on, stateAbs_on, weightsTr_off, P = utileMat

  return correlation_interburst(P, weightsTr_off,stateAbs_on, stateTr_on, 15000)
end



function (f::Distribution_RNA)(dataidx::Int, optimstruct::OptimStructWrapper)
    @unpack utileMat, maxrnaFC, model = optimstruct
    @unpack Qrna = utileMat
    return distribution_mrna(Qrna, maxrnaFC, model.nbstate)
end



#WITHOUT burst singlets


function (f::Survival_Burst)(dataidx::Int, optimstruct::OptimStructWrapper, singlets::Symbol)
    @unpack utileMat, data_fit = optimstruct
    @unpack weightsTr_on, stateTr_on, P = utileMat
    survival_burst_wosinglet(P, stateTr_on, weightsTr_on,data_fit.data[dataidx][1])  
end


function (f::Survival_InterBurst)(dataidx::Int, optimstruct::OptimStructWrapper, singlets::Symbol)
    @unpack utileMat, data_fit = optimstruct
    @unpack PabsOff, weightsTr_on_wos = utileMat
    survival_interburst_wosinglet(PabsOff, weightsTr_on_wos,data_fit.data[dataidx][1])   
end


function (f::Survival_NextBurst)(dataidx::Int, optimstruct::OptimStructWrapper, singlets::Symbol)
    @unpack utileMat, data_fit = optimstruct
    @unpack weightsAbsorbed_off_wos, PabsOff = utileMat
    survival_nextburst_wosinglet(weightsAbsorbed_off_wos,PabsOff, data_fit.data[dataidx][1])  
end



function (f::Mean_Nascent)(dataidx::Int, optimstruct::OptimStructWrapper, singlets::Symbol)
    @unpack utileMat = optimstruct
    @unpack ssp, nascentbin, stateTr = utileMat
    return mean_nascentrna_wosinglet(ssp, nascentbin, stateTr, optimstruct.maxrnaLC, optimstruct.data_fit.detectionLimitNS) 
end

function (f::Prob_Burst)(dataidx::Int, optimstruct::OptimStructWrapper, singlets::Symbol)
    @unpack utileMat = optimstruct
    @unpack sspwos,weightsPre_on_and_on, stateTr_on = utileMat
    return prob_burst_wosinglet(sspwos,weightsPre_on_and_on, stateTr_on)
end



function  (f::Correlation_InterBurst)(dataidx::Int, optimstruct::OptimStructWrapper, singlets::Symbol)
    @unpack utileMat = optimstruct
    @unpack Rn, NR, Nc, Qn, Nn, weightsTr_off_wos = utileMat
    correlation_interburst_wosinglet(Rn, NR, Nc, Qn, Nn, weightsTr_off_wos, 15000) 
end


function (f::Intensity_Burst)(dataidx::Int, optimstruct::OptimStructWrapper,singlets::Symbol)
    @unpack utileMat, data_fit = optimstruct
    @unpack rnanbvec_on, Pabs_wos, weightsON_wos, Pwos, weightsPre_on_wos, statePre_on_wos, stateTr_on  = utileMat
    return intensity_burst_wosinglet(rnanbvec_on, Pwos, Pabs_wos, weightsPre_on_wos, statePre_on_wos, stateTr_on, weightsON_wos,data_fit.data[dataidx][1])
end

function (f::Distribution_RNA)(dataidx::Int, optimstruct::OptimStructWrapper, singlets::Symbol)
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