#functions used in the optimization procedure

struct OptimStruct{DF,DI,M}
    data::DF
    dist::DI
    model::M
    AutoDiff::Bool
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


"""
    optim_function(SRange, FRange, optim_struct::OptimStruct, args...; maxrnaLC = 10, maxrnaFC = 60, freeparametersidx =[x for x in eachindex(SRange)], fixedparameters =[-1],  kwargs...)

run a parameter optimization based on several runs and return
            -   sol = solution from the optimization package
            -   bfparameters = vector of the best fit parameters
            -   minval = evalution of the optimized funct at best fit parameter values
            -   minidx = index of the optimization run that return the minval
            -   estimate_signal = namedTupled with all the observable predicted by the model at best fit parameter values

important
            -   the default optimizer method is BBO_adaptive_de_rand_1_bin_radiuslimited()
            -   the default AD method is AutoForwardDiff()
"""
function optim_function(SRange, FRange, optim_struct::OptimStruct, args...; maxrnaLC = 10, maxrnaFC = 60, freeparametersidx =[x for x in eachindex(SRange)], fixedparameters =[-1], multipledistribution=false, kwargs...)
    @unpack data, dist, model, AutoDiff = optim_struct
    if multipledistribution
        err_func = ini_optim_multipleDistribution(optim_struct, optim_struct.data.datagroup)

        data_fit = ini_data(optim_struct, FRange)

        (P,ssp, stateTr, stateTr_on, stateAbs_on, weightsTr_off,PabsOff, sspTr_Off, Pabs ) = StoThyLiveCell.mo_basics(model, ones(model.nbparameters+model.nbkini+1), maxrnaLC, data.detectionLimitLC, data.detectionLimitNS) 
        Qrna = zeros(model.nbstate*(maxrnaFC+1),model.nbstate*(maxrnaFC+1))
        utileMat = (stateTr=stateTr, stateTr_on=stateTr_on, stateAbs_on=stateAbs_on, weightsTr_off=weightsTr_off, P=P, ssp=ssp, PabsOff=PabsOff, sspTr_Off=sspTr_Off, Pabs=Pabs, Qrna=Qrna)
        optim_struct_wrapper = OptimStructWrapper{typeof(optim_struct.data),typeof(optim_struct.dist), typeof(optim_struct.model),typeof(err_func), typeof(utileMat)}(optim_struct.data, data_fit, optim_struct.dist, optim_struct.model, SRange, maxrnaLC, maxrnaFC, freeparametersidx,fixedparameters, utileMat, err_func)
    
        sol = start_optim(optim_struct_wrapper, args...; kwargs...)

        fvals = [sol[i].objective for i in eachindex(sol)] #collect all the optimization objective fct values
        minval, minidx = findmin(fvals)
        bfparameters= utiles.mergeparameter_base(fixedparameters, sol[minidx].u, freeparametersidx)
        if data.datagroup == FixedCellData()
            bfparameters_lc = vcat(bfparameters[1:3],1.)
            (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =    StoThyLiveCell.ModelOutput(optim_struct_wrapper.model,bfparameters_lc, maxrnaLC, optim_struct_wrapper.data_fit.detectionLimitLC, optim_struct_wrapper.data_fit.detectionLimitNS,10,400,400,10)  
            mrna_distribution_model = bfparameters[end-1]*bfparameters[end].*utiles.convolution(distribution_mrna(optim_struct_wrapper.model, bfparameters[1:4], optim_struct_wrapper.maxrnaFC)) .+ bfparameters[end-1]*(1-bfparameters[end]).*utiles.convolution(distribution_mrna(optim_struct_wrapper.model, vcat(bfparameters[5:7],bfparameters[4]), optim_struct_wrapper.maxrnaFC))
            mrna_distribution_model[1] = mrna_distribution_model[1] + (1-bfparameters[end-1])
        else
            @warn "optim not defined"
        end
    else
        if data.burstsinglet == :with
            if AutoDiff
                err_func = ini_optimAD(optim_struct, optim_struct.data.datagroup)
            else
                err_func = ini_optim(optim_struct, optim_struct.data.datagroup)
            end
        elseif data.burstsinglet == :without
            if AutoDiff
                err_func = ini_optim_withoutsingletAD(optim_struct, optim_struct.data.datagroup)
            else
                err_func = ini_optim_withoutsinglet(optim_struct, optim_struct.data.datagroup)
            end
        end

        data_fit = ini_data(optim_struct, FRange)

        #allocate memory for the utiles matrices
        if data.burstsinglet == :with
            if AutoDiff
                (stateTr, stateTr_on, stateAbs_on) = StoThyLiveCell.mo_basics_nonparam(model, maxrnaLC, data.detectionLimitLC, data.detectionLimitNS) 
                utileMat = (stateTr=stateTr, stateTr_on=stateTr_on, stateAbs_on=stateAbs_on,)
                optim_struct_wrapper = OptimStructWrapper{typeof(optim_struct.data),typeof(optim_struct.dist), typeof(optim_struct.model),typeof(err_func), typeof(utileMat)}(optim_struct.data, data_fit, optim_struct.dist, optim_struct.model, SRange, maxrnaLC, maxrnaFC, freeparametersidx,fixedparameters, utileMat, err_func)
            else
                (P,ssp, stateTr, stateTr_on, stateAbs_on, weightsTr_off,PabsOff, sspTr_Off, Pabs ) = StoThyLiveCell.mo_basics(model, ones(model.nbparameters+model.nbkini+1), maxrnaLC, data.detectionLimitLC, data.detectionLimitNS) 
                Qrna = zeros(model.nbstate*(maxrnaFC+1),model.nbstate*(maxrnaFC+1))
                utileMat = (stateTr=stateTr, stateTr_on=stateTr_on, stateAbs_on=stateAbs_on, weightsTr_off=weightsTr_off, P=P, ssp=ssp, PabsOff=PabsOff, sspTr_Off=sspTr_Off, Pabs=Pabs, Qrna=Qrna)
                optim_struct_wrapper = OptimStructWrapper{typeof(optim_struct.data),typeof(optim_struct.dist), typeof(optim_struct.model),typeof(err_func), typeof(utileMat)}(optim_struct.data, data_fit, optim_struct.dist, optim_struct.model, SRange, maxrnaLC, maxrnaFC, freeparametersidx,fixedparameters, utileMat, err_func)
            end
        elseif data.burstsinglet == :without
            if AutoDiff
                (nascentbin, stateTr, stateTr_on, stateAbs_on, totnbs, stateAbs_on_wos, statePre_on_wos, rnanbvec_on) = StoThyLiveCell.mo_basics_wosinglet_nonparam(model, maxrnaLC, data.detectionLimitLC, data.detectionLimitNS) 
                utileMat = (nascentbin=nascentbin, stateTr=stateTr, stateTr_on=stateTr_on, stateAbs_on=stateAbs_on, totnbs=totnbs, stateAbs_on_wos=stateAbs_on_wos, statePre_on_wos=statePre_on_wos,rnanbvec_on=rnanbvec_on)
                optim_struct_wrapper = OptimStructWrapper{typeof(optim_struct.data),typeof(optim_struct.dist), typeof(optim_struct.model),typeof(err_func), typeof(utileMat)}(optim_struct.data, data_fit, optim_struct.dist, optim_struct.model, SRange, maxrnaLC, maxrnaFC, freeparametersidx,fixedparameters, utileMat, err_func)
            else
                (nascentbin, P, ssp, stateTr, stateTr_on, stateAbs_on, totnbs, Pwos, stateAbs_on_wos, statePre_on_wos, weightsAbs_off_wos, sspTr_off_wos, weightsAbs_on, sspPreB, weightsTr_on, PabsOff, weightsTr_on_wos, weightsAbsorbed_off_wos, sspwos, weightsPre_on_and_on, Rn, NR, Nc, Qn, Nn, weightsTr_off_wos, Pabs_wos, weightsON_wos, rnanbvec_on, weightsPre_on_wos) = StoThyLiveCell.mo_basics_wosinglet(model, ones(model.nbparameters+model.nbkini+1), maxrnaLC, data.detectionLimitLC, data.detectionLimitNS) 
                Qrna = zeros(model.nbstate*(maxrnaFC+1),model.nbstate*(maxrnaFC+1))
                utileMat = (nascentbin=nascentbin, P=P, ssp=ssp, stateTr=stateTr, stateTr_on=stateTr_on, stateAbs_on=stateAbs_on, totnbs=totnbs, Pwos=Pwos, stateAbs_on_wos=stateAbs_on_wos, statePre_on_wos=statePre_on_wos, weightsAbs_off_wos=weightsAbs_off_wos, sspTr_off_wos=sspTr_off_wos, weightsAbs_on=weightsAbs_on, sspPreB=sspPreB, weightsTr_on=weightsTr_on, PabsOff=PabsOff, weightsTr_on_wos=weightsTr_on_wos, weightsAbsorbed_off_wos=weightsAbsorbed_off_wos, sspwos=sspwos, weightsPre_on_and_on=weightsPre_on_and_on, Rn=Rn, NR=NR, Nc=Nc, Qn=Qn, Nn=Nn, weightsTr_off_wos=weightsTr_off_wos, Pabs_wos=Pabs_wos, weightsON_wos=weightsON_wos, rnanbvec_on=rnanbvec_on, weightsPre_on_wos=weightsPre_on_wos, Qrna=Qrna)
                optim_struct_wrapper = OptimStructWrapper{typeof(optim_struct.data),typeof(optim_struct.dist), typeof(optim_struct.model),typeof(err_func), typeof(utileMat)}(optim_struct.data, data_fit, optim_struct.dist, optim_struct.model, SRange, maxrnaLC, maxrnaFC, freeparametersidx,fixedparameters, utileMat, err_func)
            end
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
                (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =    StoThyLiveCell.ModelOutput(optim_struct_wrapper.model, bfparameters, maxrnaLC, optim_struct_wrapper.data_fit.detectionLimitLC, optim_struct_wrapper.data_fit.detectionLimitNS, 10,400,400,10) 
                bfparameters_mrna = vcat(bfparameters[1:end-1],.01)
                mrna_distribution_model = StoThyLiveCell.distribution_mrna(optim_struct_wrapper.model, bfparameters_mrna, optim_struct_wrapper.maxrnaFC) 
            elseif data.datagroup == FixedCellData()
                bfparameters_lc = vcat(bfparameters[1:end-1],1.)
                (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =    StoThyLiveCell.ModelOutput(optim_struct_wrapper.model,bfparameters_lc, maxrnaLC, optim_struct_wrapper.data_fit.detectionLimitLC, optim_struct_wrapper.data_fit.detectionLimitNS,10,400,400,10)  
                mrna_distribution_model = StoThyLiveCell.distribution_mrna(optim_struct_wrapper.model, bfparameters, optim_struct_wrapper.maxrnaFC) 
            elseif data.datagroup == FixedAndLiveCellData()
                (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =    StoThyLiveCell.ModelOutput(optim_struct_wrapper.model, bfparameters[1:end-1], maxrnaLC, optim_struct_wrapper.data_fit.detectionLimitLC, optim_struct_wrapper.data_fit.detectionLimitNS, 10,400,400,10) 
                bfparameters_mrna = vcat(bfparameters[1:end-1],bfparameters[end])
                mrna_distribution_model = StoThyLiveCell.distribution_mrna(optim_struct_wrapper.model, bfparameters_mrna, optim_struct_wrapper.maxrnaFC) 
            end
        elseif data.burstsinglet == :without
            if data.datagroup == LiveCellData()
                (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =    StoThyLiveCell.ModelOutput_wosinglet(optim_struct_wrapper.model, bfparameters, maxrnaLC, optim_struct_wrapper.data_fit.detectionLimitLC, optim_struct_wrapper.data_fit.detectionLimitNS, 10,400,400,10) 
                bfparameters_mrna = vcat(bfparameters[1:end-1],.01)
                mrna_distribution_model = StoThyLiveCell.distribution_mrna(optim_struct_wrapper.model, bfparameters_mrna, optim_struct_wrapper.maxrnaFC) 
            elseif data.datagroup == FixedCellData()
                bfparameters_lc = vcat(bfparameters[1:end-1],1.)
                (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =    StoThyLiveCell.ModelOutput(optim_struct_wrapper.model,bfparameters_lc, maxrnaLC, optim_struct_wrapper.data_fit.detectionLimitLC, optim_struct_wrapper.data_fit.detectionLimitNS, 10,400,400,10) 
                mrna_distribution_model = StoThyLiveCell.distribution_mrna(optim_struct_wrapper.model, bfparameters, optim_struct_wrapper.maxrnaFC) 
            elseif data.datagroup == FixedAndLiveCellData()
                (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =    StoThyLiveCell.ModelOutput_wosinglet(optim_struct_wrapper.model, bfparameters[1:end-1], maxrnaLC, optim_struct_wrapper.data_fit.detectionLimitLC, optim_struct_wrapper.data_fit.detectionLimitNS, 10,400,400,10) 
                bfparameters_mrna = vcat(bfparameters[1:end-1],bfparameters[end])
                mrna_distribution_model = StoThyLiveCell.distribution_mrna(optim_struct_wrapper.model, bfparameters_mrna, optim_struct_wrapper.maxrnaFC) 
            end
        end
    end
    estimate_signal = (survival_burst = survival_burst, survival_interburst = survival_interburst, survival_nextburst = survival_nextburst, prob_burst = prob_burst, mean_nascentrna = mean_nascentrna, correlation_interburst = correlation_interburst, intensity_burst = intensity_burst, mrna_distribution = mrna_distribution_model)    
    return sol, bfparameters, minval, minidx, estimate_signal
end


"""
    optim_function(SRange, FRange, optim_struct::OptimStruct, args...; maxrnaLC = 10, maxrnaFC = 60, freeparametersidx =[x for x in eachindex(SRange)], fixedparameters =[-1],  kwargs...)

This version of the function required a function that reparameterize the model
"""
function optim_function(SRange, FRange, optim_struct::OptimStruct, f_parameter, args...; maxrnaLC = 10, maxrnaFC = 60, freeparametersidx =[x for x in eachindex(SRange)], fixedparameters =[-1],  kwargs...)
    @unpack data, dist, model, AutoDiff = optim_struct

    if data.burstsinglet == :with
        if AutoDiff
            err_func = ini_optimAD(optim_struct, optim_struct.data.datagroup, f_parameter)
        else
            err_func = ini_optim(optim_struct, optim_struct.data.datagroup, f_parameter)
        end
    elseif data.burstsinglet == :without
        if AutoDiff
            err_func = ini_optim_withoutsingletAD(optim_struct, optim_struct.data.datagroup, f_parameter)
        else
            err_func = ini_optim_withoutsinglet(optim_struct, optim_struct.data.datagroup, f_parameter)
        end
    end

    data_fit = ini_data(optim_struct, FRange)

    #allocate memory for the utiles matrices
    if data.burstsinglet == :with
        if AutoDiff
            (stateTr, stateTr_on, stateAbs_on) = StoThyLiveCell.mo_basics_nonparam(model, maxrnaLC, data.detectionLimitLC, data.detectionLimitNS) 
            utileMat = (stateTr=stateTr, stateTr_on=stateTr_on, stateAbs_on=stateAbs_on,)
            optim_struct_wrapper = OptimStructWrapper{typeof(optim_struct.data),typeof(optim_struct.dist), typeof(optim_struct.model),typeof(err_func), typeof(utileMat)}(optim_struct.data, data_fit, optim_struct.dist, optim_struct.model, SRange, maxrnaLC, maxrnaFC, freeparametersidx,fixedparameters, utileMat, err_func)
        else
            (P,ssp, stateTr, stateTr_on, stateAbs_on, weightsTr_off,PabsOff, sspTr_Off, Pabs ) = StoThyLiveCell.mo_basics(model, ones(model.nbparameters+model.nbkini+1), maxrnaLC, data.detectionLimitLC, data.detectionLimitNS) 
            Qrna = zeros(model.nbstate*(maxrnaFC+1),model.nbstate*(maxrnaFC+1))
            utileMat = (stateTr=stateTr, stateTr_on=stateTr_on, stateAbs_on=stateAbs_on, weightsTr_off=weightsTr_off, P=P, ssp=ssp, PabsOff=PabsOff, sspTr_Off=sspTr_Off, Pabs=Pabs, Qrna=Qrna)
            optim_struct_wrapper = OptimStructWrapper{typeof(optim_struct.data),typeof(optim_struct.dist), typeof(optim_struct.model),typeof(err_func), typeof(utileMat)}(optim_struct.data, data_fit, optim_struct.dist, optim_struct.model, SRange, maxrnaLC, maxrnaFC, freeparametersidx,fixedparameters, utileMat, err_func)
        end
    elseif data.burstsinglet == :without
        if AutoDiff
            (nascentbin, stateTr, stateTr_on, stateAbs_on, totnbs, stateAbs_on_wos, statePre_on_wos, rnanbvec_on) = StoThyLiveCell.mo_basics_wosinglet_nonparam(model, maxrnaLC, data.detectionLimitLC, data.detectionLimitNS) 
            utileMat = (nascentbin=nascentbin, stateTr=stateTr, stateTr_on=stateTr_on, stateAbs_on=stateAbs_on, totnbs=totnbs, stateAbs_on_wos=stateAbs_on_wos, statePre_on_wos=statePre_on_wos,rnanbvec_on=rnanbvec_on)
            optim_struct_wrapper = OptimStructWrapper{typeof(optim_struct.data),typeof(optim_struct.dist), typeof(optim_struct.model),typeof(err_func), typeof(utileMat)}(optim_struct.data, data_fit, optim_struct.dist, optim_struct.model, SRange, maxrnaLC, maxrnaFC, freeparametersidx,fixedparameters, utileMat, err_func)
        else
            (nascentbin, P, ssp, stateTr, stateTr_on, stateAbs_on, totnbs, Pwos, stateAbs_on_wos, statePre_on_wos, weightsAbs_off_wos, sspTr_off_wos, weightsAbs_on, sspPreB, weightsTr_on, PabsOff, weightsTr_on_wos, weightsAbsorbed_off_wos, sspwos, weightsPre_on_and_on, Rn, NR, Nc, Qn, Nn, weightsTr_off_wos, Pabs_wos, weightsON_wos, rnanbvec_on, weightsPre_on_wos) = StoThyLiveCell.mo_basics_wosinglet(model, ones(model.nbparameters+model.nbkini+1), maxrnaLC, data.detectionLimitLC, data.detectionLimitNS) 
            Qrna = zeros(model.nbstate*(maxrnaFC+1),model.nbstate*(maxrnaFC+1))
            utileMat = (nascentbin=nascentbin, P=P, ssp=ssp, stateTr=stateTr, stateTr_on=stateTr_on, stateAbs_on=stateAbs_on, totnbs=totnbs, Pwos=Pwos, stateAbs_on_wos=stateAbs_on_wos, statePre_on_wos=statePre_on_wos, weightsAbs_off_wos=weightsAbs_off_wos, sspTr_off_wos=sspTr_off_wos, weightsAbs_on=weightsAbs_on, sspPreB=sspPreB, weightsTr_on=weightsTr_on, PabsOff=PabsOff, weightsTr_on_wos=weightsTr_on_wos, weightsAbsorbed_off_wos=weightsAbsorbed_off_wos, sspwos=sspwos, weightsPre_on_and_on=weightsPre_on_and_on, Rn=Rn, NR=NR, Nc=Nc, Qn=Qn, Nn=Nn, weightsTr_off_wos=weightsTr_off_wos, Pabs_wos=Pabs_wos, weightsON_wos=weightsON_wos, rnanbvec_on=rnanbvec_on, weightsPre_on_wos=weightsPre_on_wos, Qrna=Qrna)
            optim_struct_wrapper = OptimStructWrapper{typeof(optim_struct.data),typeof(optim_struct.dist), typeof(optim_struct.model),typeof(err_func), typeof(utileMat)}(optim_struct.data, data_fit, optim_struct.dist, optim_struct.model, SRange, maxrnaLC, maxrnaFC, freeparametersidx,fixedparameters, utileMat, err_func)
        end
    end
    #run the optimization
    sol = start_optim(optim_struct_wrapper, args...; kwargs...)

    #bestfit parameters
    fvals = [sol[i].objective for i in eachindex(sol)] #collect all the optimization objective fct values
    minval, minidx = findmin(fvals)
    bfparameters= f_parameter(utiles.mergeparameter_base(fixedparameters, sol[minidx].u, freeparametersidx))
    bfreparameters = utiles.mergeparameter_base(fixedparameters, sol[minidx].u, freeparametersidx)
    #bestfit signal
    if data.burstsinglet == :with
        if data.datagroup == LiveCellData()
            (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =    StoThyLiveCell.ModelOutput(optim_struct_wrapper.model, bfparameters, maxrnaLC, optim_struct_wrapper.data_fit.detectionLimitLC, optim_struct_wrapper.data_fit.detectionLimitNS, 10,400,400,10) 
            bfparameters_mrna = vcat(bfparameters[1:end-1],.01)
            mrna_distribution_model = StoThyLiveCell.distribution_mrna(optim_struct_wrapper.model, bfparameters_mrna, optim_struct_wrapper.maxrnaFC) 
        elseif data.datagroup == FixedCellData()
            bfparameters_lc = vcat(bfparameters[1:end-1],1.)
            (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =    StoThyLiveCell.ModelOutput(optim_struct_wrapper.model,bfparameters_lc, maxrnaLC, optim_struct_wrapper.data_fit.detectionLimitLC, optim_struct_wrapper.data_fit.detectionLimitNS,10,400,400,10)  
            mrna_distribution_model = StoThyLiveCell.distribution_mrna(optim_struct_wrapper.model, bfparameters, optim_struct_wrapper.maxrnaFC) 
        elseif data.datagroup == FixedAndLiveCellData()
            (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =    StoThyLiveCell.ModelOutput(optim_struct_wrapper.model, bfparameters[1:end-1], maxrnaLC, optim_struct_wrapper.data_fit.detectionLimitLC, optim_struct_wrapper.data_fit.detectionLimitNS, 10,400,400,10) 
            bfparameters_mrna = vcat(bfparameters[1:end-1],bfparameters[end])
            mrna_distribution_model = StoThyLiveCell.distribution_mrna(optim_struct_wrapper.model, bfparameters_mrna, optim_struct_wrapper.maxrnaFC) 
        end
    elseif data.burstsinglet == :without
        if data.datagroup == LiveCellData()
            (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =    StoThyLiveCell.ModelOutput_wosinglet(optim_struct_wrapper.model, bfparameters, maxrnaLC, optim_struct_wrapper.data_fit.detectionLimitLC, optim_struct_wrapper.data_fit.detectionLimitNS, 10,400,400,10) 
            bfparameters_mrna = vcat(bfparameters[1:end-1],.01)
            mrna_distribution_model = StoThyLiveCell.distribution_mrna(optim_struct_wrapper.model, bfparameters_mrna, optim_struct_wrapper.maxrnaFC) 
        elseif data.datagroup == FixedCellData()
            bfparameters_lc = vcat(bfparameters[1:end-1],1.)
            (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =    StoThyLiveCell.ModelOutput(optim_struct_wrapper.model,bfparameters_lc, maxrnaLC, optim_struct_wrapper.data_fit.detectionLimitLC, optim_struct_wrapper.data_fit.detectionLimitNS, 10,400,400,10) 
            mrna_distribution_model = StoThyLiveCell.distribution_mrna(optim_struct_wrapper.model, bfparameters, optim_struct_wrapper.maxrnaFC) 
        elseif data.datagroup == FixedAndLiveCellData()
            (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =    StoThyLiveCell.ModelOutput_wosinglet(optim_struct_wrapper.model, bfparameters[1:end-1], maxrnaLC, optim_struct_wrapper.data_fit.detectionLimitLC, optim_struct_wrapper.data_fit.detectionLimitNS, 10,400,400,10) 
            bfparameters_mrna = vcat(bfparameters[1:end-1],bfparameters[end])
            mrna_distribution_model = StoThyLiveCell.distribution_mrna(optim_struct_wrapper.model, bfparameters_mrna, optim_struct_wrapper.maxrnaFC) 
        end
    end
    estimate_signal = (survival_burst = survival_burst, survival_interburst = survival_interburst, survival_nextburst = survival_nextburst, prob_burst = prob_burst, mean_nascentrna = mean_nascentrna, correlation_interburst = correlation_interburst, intensity_burst = intensity_burst, mrna_distribution = mrna_distribution_model)    
    return sol, bfreparameters, minval, minidx, estimate_signal
end

function optim_function_multipleModels(SRange, FRange_model1, FRange_model2, optim_struct_model1::OptimStruct,optim_struct_model2::OptimStruct, args...; maxrnaLC = 10, maxrnaFC = 60, freeparametersidx_m1 =[-1], freeparametersidx_m2 =[-1], fixedparameters_m1 =[-1], fixedparameters_m2 =[-1], paramidx_m1=[-1], paramidx_m2=[-1], pm1_infreeparam=[-1], pm2_infreeparam=[-1],  kwargs...)
    @unpack AutoDiff = optim_struct_model1


    if optim_struct_model1.data.burstsinglet == :with
        if AutoDiff
            err_func_m1 = ini_optimAD(optim_struct_model1, optim_struct_model1.data.datagroup)
            err_func_m2 = ini_optimAD(optim_struct_model2, optim_struct_model2.data.datagroup)
        else
            err_func_m1 = ini_optim(optim_struct_model1, optim_struct_model1.data.datagroup)
            err_func_m2 = ini_optim(optim_struct_model2, optim_struct_model2.data.datagroup)
        end
    elseif optim_struct_model1.data.burstsinglet == :without
        if AutoDiff
            err_func_m1= ini_optim_withoutsingletAD(optim_struct_model1, optim_struct_model1.data.datagroup)
            err_func_m2= ini_optim_withoutsingletAD(optim_struct_model2, optim_struct_model2.data.datagroup)
        else
            err_func_m1 = ini_optim_withoutsinglet(optim_struct_model1, optim_struct_model1.data.datagroup)
            err_func_m2 = ini_optim_withoutsinglet(optim_struct_model2, optim_struct_model2.data.datagroup)
        end
    end

    data_fit_m1 = ini_data(optim_struct_model1, FRange_model1)
    data_fit_m2 = ini_data(optim_struct_model2, FRange_model2)

    freeparametersidx_model1 = freeparametersidx_m1
    fixedparameters_model1 = fixedparameters_m1
    freeparametersidx_model2 = freeparametersidx_m2
    fixedparameters_model2 = fixedparameters_m2
    #allocate memory for the utiles matrices
    if optim_struct_model1.data.burstsinglet == :with
        if AutoDiff
            (stateTr, stateTr_on, stateAbs_on) = StoThyLiveCell.mo_basics_nonparam(optim_struct_model1.model, maxrnaLC, optim_struct_model1.data.detectionLimitLC, optim_struct_model1.data.detectionLimitNS) 
            utileMat = (stateTr=stateTr, stateTr_on=stateTr_on, stateAbs_on=stateAbs_on,)
            optim_struct_wrapper_m1 = OptimStructWrapper{typeof(optim_struct_model1.data),typeof(optim_struct_model1.dist), typeof(optim_struct_model1.model),typeof(err_func_m1), typeof(utileMat)}(optim_struct_model1.data, data_fit_m1, optim_struct_model1.dist, optim_struct_model1.model, SRange, maxrnaLC, maxrnaFC, freeparametersidx_model1,fixedparameters_model1, utileMat, err_func_m1)
            (stateTr, stateTr_on, stateAbs_on) = StoThyLiveCell.mo_basics_nonparam(optim_struct_model2.model, maxrnaLC, optim_struct_model2.data.detectionLimitLC, optim_struct_model2.data.detectionLimitNS) 
            utileMat = (stateTr=stateTr, stateTr_on=stateTr_on, stateAbs_on=stateAbs_on,)
            optim_struct_wrapper_m2 = OptimStructWrapper{typeof(optim_struct_model2.data),typeof(optim_struct_model2.dist), typeof(optim_struct_model2.model),typeof(err_func_m2), typeof(utileMat)}(optim_struct_model2.data, data_fit_m2, optim_struct_model2.dist, optim_struct_model2.model, SRange, maxrnaLC, maxrnaFC, freeparametersidx_model2,fixedparameters_model2, utileMat, err_func_m2)
        else
            (P,ssp, stateTr, stateTr_on, stateAbs_on, weightsTr_off,PabsOff, sspTr_Off, Pabs ) = StoThyLiveCell.mo_basics(optim_struct_model1.model, ones(optim_struct_model1.model.nbparameters+model.nbkini+1), maxrnaLC, optim_struct_model1.data.detectionLimitLC, optim_struct_model1.data.detectionLimitNS) 
            Qrna = zeros(optim_struct_model1.model.nbstate*(maxrnaFC+1),optim_struct_model1.model.nbstate*(maxrnaFC+1))
            utileMat = (stateTr=stateTr, stateTr_on=stateTr_on, stateAbs_on=stateAbs_on, weightsTr_off=weightsTr_off, P=P, ssp=ssp, PabsOff=PabsOff, sspTr_Off=sspTr_Off, Pabs=Pabs, Qrna=Qrna)
            optim_struct_wrapper_m1 = OptimStructWrapper{typeof(optim_struct_model1.data),typeof(optim_struct_model1.dist), typeof(optim_struct_model1.model),typeof(err_func_m1), typeof(utileMat)}(optim_struct_model1.data, data_fit_m1, optim_struct_model1.dist, optim_struct_model1.model, SRange, maxrnaLC, maxrnaFC, freeparametersidx_model1,fixedparameters_model1, utileMat, err_func_m1)
            (P,ssp, stateTr, stateTr_on, stateAbs_on, weightsTr_off,PabsOff, sspTr_Off, Pabs ) = StoThyLiveCell.mo_basics(optim_struct_model2.model, ones(optim_struct_model2.model.nbparameters+model.nbkini+1), maxrnaLC, optim_struct_model2.data.detectionLimitLC, optim_struct_model2.data.detectionLimitNS) 
            Qrna = zeros(optim_struct_model2.model.nbstate*(maxrnaFC+1),optim_struct_model2.model.nbstate*(maxrnaFC+1))
            utileMat = (stateTr=stateTr, stateTr_on=stateTr_on, stateAbs_on=stateAbs_on, weightsTr_off=weightsTr_off, P=P, ssp=ssp, PabsOff=PabsOff, sspTr_Off=sspTr_Off, Pabs=Pabs, Qrna=Qrna)
            optim_struct_wrapper_m2 = OptimStructWrapper{typeof(optim_struct_model2.data),typeof(optim_struct_model2.dist), typeof(optim_struct_model2.model),typeof(err_func_m2), typeof(utileMat)}(optim_struct_model2.data, data_fit_m2, optim_struct_model2.dist, optim_struct_model2.model, SRange, maxrnaLC, maxrnaFC, freeparametersidx_model2,fixedparameters_model2, utileMat, err_func_m2)
        end
    elseif optim_struct_model1.data.burstsinglet == :without
        if AutoDiff
            (nascentbin, stateTr, stateTr_on, stateAbs_on, totnbs, stateAbs_on_wos, statePre_on_wos, rnanbvec_on) = StoThyLiveCell.mo_basics_wosinglet_nonparam(optim_struct_model1.model, maxrnaLC, optim_struct_model1.data.detectionLimitLC, optim_struct_model1.data.detectionLimitNS) 
            utileMat = (nascentbin=nascentbin, stateTr=stateTr, stateTr_on=stateTr_on, stateAbs_on=stateAbs_on, totnbs=totnbs, stateAbs_on_wos=stateAbs_on_wos, statePre_on_wos=statePre_on_wos,rnanbvec_on=rnanbvec_on)
            optim_struct_wrapper_m1 = OptimStructWrapper{typeof(optim_struct_model1.data),typeof(optim_struct_model1.dist), typeof(optim_struct_model1.model),typeof(err_func_m1), typeof(utileMat)}(optim_struct_model1.data, data_fit_m1, optim_struct_model1.dist, optim_struct_model1.model, SRange, maxrnaLC, maxrnaFC, freeparametersidx_model1,fixedparameters_model1, utileMat, err_func_m1)
            (nascentbin, stateTr, stateTr_on, stateAbs_on, totnbs, stateAbs_on_wos, statePre_on_wos, rnanbvec_on) = StoThyLiveCell.mo_basics_wosinglet_nonparam(optim_struct_model2.model, maxrnaLC, optim_struct_model2.data.detectionLimitLC, optim_struct_model2.data.detectionLimitNS) 
            utileMat = (nascentbin=nascentbin, stateTr=stateTr, stateTr_on=stateTr_on, stateAbs_on=stateAbs_on, totnbs=totnbs, stateAbs_on_wos=stateAbs_on_wos, statePre_on_wos=statePre_on_wos,rnanbvec_on=rnanbvec_on)
            optim_struct_wrapper_m2 = OptimStructWrapper{typeof(optim_struct_model2.data),typeof(optim_struct_model2.dist), typeof(optim_struct_model2.model),typeof(err_func_m2), typeof(utileMat)}(optim_struct_model2.data, data_fit_m2, optim_struct_model2.dist, optim_struct_model2.model, SRange, maxrnaLC, maxrnaFC, freeparametersidx_model2,fixedparameters_model2, utileMat, err_func_m2)

        else
            (nascentbin, P, ssp, stateTr, stateTr_on, stateAbs_on, totnbs, Pwos, stateAbs_on_wos, statePre_on_wos, weightsAbs_off_wos, sspTr_off_wos, weightsAbs_on, sspPreB, weightsTr_on, PabsOff, weightsTr_on_wos, weightsAbsorbed_off_wos, sspwos, weightsPre_on_and_on, Rn, NR, Nc, Qn, Nn, weightsTr_off_wos, Pabs_wos, weightsON_wos, rnanbvec_on, weightsPre_on_wos) = StoThyLiveCell.mo_basics_wosinglet(optim_struct_model1.model, ones(optim_struct_model1.model.nbparameters+optim_struct_model1.model.nbkini+1), maxrnaLC, optim_struct_model1.data.detectionLimitLC, optim_struct_model1.data.detectionLimitNS) 
            Qrna = zeros(optim_struct_model1.model.nbstate*(maxrnaFC+1),optim_struct_model1.model.nbstate*(maxrnaFC+1))
            utileMat = (nascentbin=nascentbin, P=P, ssp=ssp, stateTr=stateTr, stateTr_on=stateTr_on, stateAbs_on=stateAbs_on, totnbs=totnbs, Pwos=Pwos, stateAbs_on_wos=stateAbs_on_wos, statePre_on_wos=statePre_on_wos, weightsAbs_off_wos=weightsAbs_off_wos, sspTr_off_wos=sspTr_off_wos, weightsAbs_on=weightsAbs_on, sspPreB=sspPreB, weightsTr_on=weightsTr_on, PabsOff=PabsOff, weightsTr_on_wos=weightsTr_on_wos, weightsAbsorbed_off_wos=weightsAbsorbed_off_wos, sspwos=sspwos, weightsPre_on_and_on=weightsPre_on_and_on, Rn=Rn, NR=NR, Nc=Nc, Qn=Qn, Nn=Nn, weightsTr_off_wos=weightsTr_off_wos, Pabs_wos=Pabs_wos, weightsON_wos=weightsON_wos, rnanbvec_on=rnanbvec_on, weightsPre_on_wos=weightsPre_on_wos, Qrna=Qrna)
            optim_struct_wrapper_m1 = OptimStructWrapper{typeof(optim_struct_model1.data),typeof(optim_struct_model1.dist), typeof(optim_struct_model1.model),typeof(err_func_m1), typeof(utileMat)}(optim_struct_model1.data, data_fit_m1, optim_struct_model1.dist, optim_struct_model1.model, SRange, maxrnaLC, maxrnaFC, freeparametersidx_model1,fixedparameters_model1, utileMat, err_func_m1)
            (nascentbin, P, ssp, stateTr, stateTr_on, stateAbs_on, totnbs, Pwos, stateAbs_on_wos, statePre_on_wos, weightsAbs_off_wos, sspTr_off_wos, weightsAbs_on, sspPreB, weightsTr_on, PabsOff, weightsTr_on_wos, weightsAbsorbed_off_wos, sspwos, weightsPre_on_and_on, Rn, NR, Nc, Qn, Nn, weightsTr_off_wos, Pabs_wos, weightsON_wos, rnanbvec_on, weightsPre_on_wos) = StoThyLiveCell.mo_basics_wosinglet(optim_struct_model2.model, ones(optim_struct_model2.model.nbparameters+optim_struct_model2.model.nbkini+1), maxrnaLC, optim_struct_model2.data.detectionLimitLC, optim_struct_model2.data.detectionLimitNS) 
            Qrna = zeros(optim_struct_model2.model.nbstate*(maxrnaFC+1),optim_struct_model2.model.nbstate*(maxrnaFC+1))
            utileMat = (nascentbin=nascentbin, P=P, ssp=ssp, stateTr=stateTr, stateTr_on=stateTr_on, stateAbs_on=stateAbs_on, totnbs=totnbs, Pwos=Pwos, stateAbs_on_wos=stateAbs_on_wos, statePre_on_wos=statePre_on_wos, weightsAbs_off_wos=weightsAbs_off_wos, sspTr_off_wos=sspTr_off_wos, weightsAbs_on=weightsAbs_on, sspPreB=sspPreB, weightsTr_on=weightsTr_on, PabsOff=PabsOff, weightsTr_on_wos=weightsTr_on_wos, weightsAbsorbed_off_wos=weightsAbsorbed_off_wos, sspwos=sspwos, weightsPre_on_and_on=weightsPre_on_and_on, Rn=Rn, NR=NR, Nc=Nc, Qn=Qn, Nn=Nn, weightsTr_off_wos=weightsTr_off_wos, Pabs_wos=Pabs_wos, weightsON_wos=weightsON_wos, rnanbvec_on=rnanbvec_on, weightsPre_on_wos=weightsPre_on_wos, Qrna=Qrna)
            optim_struct_wrapper_m2 = OptimStructWrapper{typeof(optim_struct_model2.data),typeof(optim_struct_model2.dist), typeof(optim_struct_model2.model),typeof(err_func_m2), typeof(utileMat)}(optim_struct_model2.data, data_fit_m2, optim_struct_model2.dist, optim_struct_model2.model, SRange, maxrnaLC, maxrnaFC, freeparametersidx_model2,fixedparameters_model2, utileMat, err_func_m2)
        end
    end
    #run the optimization
    sol = start_optim((optim_struct_wrapper_m1,optim_struct_wrapper_m2,),args...; paramidx_m1=paramidx_m1, paramidx_m2=paramidx_m2,pm1_infreeparam=pm1_infreeparam, pm2_infreeparam=pm2_infreeparam, kwargs...)

    #bestfit parameters
    fvals = [sol[i].objective for i in eachindex(sol)] #collect all the optimization objective fct values
    minval, minidx = findmin(fvals)
    bfparameters_m1= utiles.mergeparameter_base(fixedparameters_model1, sol[minidx].u[pm1_infreeparam], freeparametersidx_m1)
    bfparameters_m2= utiles.mergeparameter_base(fixedparameters_model2, sol[minidx].u[pm2_infreeparam], freeparametersidx_m2)

    #bestfit signal
    if optim_struct_model1.data.burstsinglet == :with
        if optim_struct_model1.data.datagroup == LiveCellData()
            (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =    StoThyLiveCell.ModelOutput(optim_struct_wrapper_m1.model, bfparameters_m1, maxrnaLC, optim_struct_wrapper_m1.data_fit.detectionLimitLC, optim_struct_wrapper_m1.data_fit.detectionLimitNS, 10,400,400,10) 
            bfparameters_mrna = vcat(bfparameters_m1[1:end-1],.01)
            mrna_distribution_model = StoThyLiveCell.distribution_mrna(optim_struct_wrapper_m1.model, bfparameters_mrna, optim_struct_wrapper_m1.maxrnaFC) 
        elseif optim_struct_model1.data.datagroup == FixedCellData()
            bfparameters_lc = vcat(bfparameters_m1[1:end-1],1.)
            (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =    StoThyLiveCell.ModelOutput(optim_struct_wrapper_m1.model,bfparameters_lc, maxrnaLC, optim_struct_wrapper_m1.data_fit.detectionLimitLC, optim_struct_wrapper_m1.data_fit.detectionLimitNS, 10,400,400,10) 
            mrna_distribution_model = StoThyLiveCell.distribution_mrna(optim_struct_wrapper_m1.model, bfparameters, optim_struct_wrapper_m1.maxrnaFC) 
        elseif optim_struct_model1.data.datagroup == FixedAndLiveCellData()
            (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =    StoThyLiveCell.ModelOutput(optim_struct_wrapper_m1.model, bfparameters_m1[1:end-1], maxrnaLC, optim_struct_wrapper_m1.data_fit.detectionLimitLC, optim_struct_wrapper_m1.data_fit.detectionLimitNS, 10,400,400,10)  
            bfparameters_mrna = vcat(bfparameters_m1[1:end-1],bfparameters_m1[end])
            mrna_distribution_model = StoThyLiveCell.distribution_mrna(optim_struct_wrapper_m1.model, bfparameters_mrna, optim_struct_wrapper_m1.maxrnaFC) 
        end
    elseif optim_struct_model1.data.burstsinglet == :without
        if optim_struct_model1.data.datagroup == LiveCellData()
            (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =    StoThyLiveCell.ModelOutput_wosinglet(optim_struct_wrapper_m1.model, bfparameters_m1, maxrnaLC, optim_struct_wrapper_m1.data_fit.detectionLimitLC, optim_struct_wrapper_m1.data_fit.detectionLimitNS, 10,400,400,10) 
            bfparameters_mrna = vcat(bfparameters_m1[1:end-1],.01)
            mrna_distribution_model = StoThyLiveCell.distribution_mrna(optim_struct_wrapper_m1.model, bfparameters_mrna, optim_struct_wrapper_m1.maxrnaFC) 
        elseif optim_struct_model1.data.datagroup == FixedCellData()
            bfparameters_lc = vcat(bfparameters_m1[1:end-1],1.)
            (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =    StoThyLiveCell.ModelOutput(optim_struct_wrapper_m1.model,bfparameters_lc, maxrnaLC, optim_struct_wrapper_m1.data_fit.detectionLimitLC, optim_struct_wrapper_m1.data_fit.detectionLimitNS, 10,400,400,10)  
            mrna_distribution_model = StoThyLiveCell.distribution_mrna(optim_struct_wrapper_m1.model, bfparameters, optim_struct_wrapper_m1.maxrnaFC) 
        elseif optim_struct_model1.data.datagroup == FixedAndLiveCellData()
            (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =    StoThyLiveCell.ModelOutput_wosinglet(optim_struct_wrapper_m1.model, bfparameters_m1[1:end-1], maxrnaLC, optim_struct_wrapper_m1data_fit.detectionLimitLC, optim_struct_wrapper_m1.data_fit.detectionLimitNS, 10,400,400,10) 
            bfparameters_mrna = vcat(bfparameters_m1[1:end-1],bfparameters_m1[end])
            mrna_distribution_model = StoThyLiveCell.distribution_mrna(optim_struct_wrapper_m1.model, bfparameters_mrna, optim_struct_wrapper_m1.maxrnaFC) 
        end
    end
    estimate_signal = (survival_burst = survival_burst, survival_interburst = survival_interburst, survival_nextburst = survival_nextburst, prob_burst = prob_burst, mean_nascentrna = mean_nascentrna, correlation_interburst = correlation_interburst, intensity_burst = intensity_burst, mrna_distribution = mrna_distribution_model)    
    return sol, (bfparameters_m1,bfparameters_m2), minval, minidx, estimate_signal
end


function start_optim(optim_struct_wrapper::OptimStructWrapper, args...; NbOptim::Int=1, maxtime::Int=1, maxiters::Int=1 , initialparameters = [],Method=BBO_adaptive_de_rand_1_bin_radiuslimited(), ADmethod=AutoForwardDiff(), pathToLog="", kwargs...)
    @unpack SRange, err_func = optim_struct_wrapper

    lbfull = [SRange[i][1] for i in eachindex(SRange)]
    ubfull = [SRange[i][2] for i in eachindex(SRange)]
    lb = lbfull[optim_struct_wrapper.freeparametersidx]
    ub = ubfull[optim_struct_wrapper.freeparametersidx]
    db = ub - lb
    sol = []
    if isempty(initialparameters)
        for i = 1: NbOptim
            u0 = lb .+ rand(length(lb)).*db
            optprob = OptimizationFunction(err_func, ADmethod);
            prob = OptimizationProblem(optprob, u0, optim_struct_wrapper, lb = lb, ub = ub)
            # Import a solver package and solve the optimization problem
            push!(sol, solve(prob, Method; maxtime = maxtime, maxiters = maxiters));
            open("$(pathToLog)log_bestfits.txt", "a") do io
                writedlm(io, sol[end].u')
            end
            open("$(pathToLog)log_fval.txt", "a") do io
                writedlm(io, sol[end].objective)    
            end
            open("$(pathToLog)log_iniparam.txt", "a") do io
                writedlm(io, u0')    
            end
        end
    else
        for row in eachrow(initialparameters)
            u0 = Vector(row)
            optprob = OptimizationFunction(err_func, ADmethod);
            prob = OptimizationProblem(optprob, u0, optim_struct_wrapper, lb = lb, ub = ub)
            # Import a solver package and solve the optimization problem
            push!(sol, solve(prob, Method; maxtime = maxtime, maxiters = maxiters));
            open("$(pathToLog)log_bestfits.txt", "a") do io
                writedlm(io, sol[end].u')
            end
            open("$(pathToLog)log_fval.txt", "a") do io
                writedlm(io, sol[end].objective)    
            end
            open("$(pathToLog)log_iniparam.txt", "a") do io
                writedlm(io, u0')    
            end
        end
    end
    return sol
end


function start_optim(optim_struct_wrapper::Tuple{OptimStructWrapper,OptimStructWrapper}, args...; paramidx_m1=[-1], paramidx_m2=[-1],pm1_infreeparam=[-1], pm2_infreeparam=[-1], NbOptim::Int=1, maxtime::Int=1, maxiters::Int=1 , initialparameters = [],Method=BBO_adaptive_de_rand_1_bin_radiuslimited(), ADmethod=AutoForwardDiff(), pathToLog="", kwargs...)
    @unpack SRange = optim_struct_wrapper[1]

    function err_func(params,optim_struct_wrapper::Tuple{OptimStructWrapper,OptimStructWrapper})
        return optim_struct_wrapper[1].err_func(params[pm1_infreeparam], optim_struct_wrapper[1]) + optim_struct_wrapper[2].err_func(params[pm2_infreeparam], optim_struct_wrapper[2])
    end

    lbfull = [SRange[i][1] for i in eachindex(SRange)]
    ubfull = [SRange[i][2] for i in eachindex(SRange)]
    lb = lbfull[union(paramidx_m1[optim_struct_wrapper[1].freeparametersidx],paramidx_m2[optim_struct_wrapper[2].freeparametersidx])]
    ub = ubfull[union(paramidx_m1[optim_struct_wrapper[1].freeparametersidx],paramidx_m2[optim_struct_wrapper[2].freeparametersidx])]
    db = ub - lb
    sol = []
    if isempty(initialparameters)
        for i = 1: NbOptim
            u0 = lb .+ rand(length(lb)).*db
            optprob = OptimizationFunction(err_func, ADmethod);
            prob = OptimizationProblem(optprob, u0, optim_struct_wrapper, lb = lb, ub = ub)
            # Import a solver package and solve the optimization problem
            push!(sol, solve(prob, Method; maxtime = maxtime, maxiters = maxiters));
            open("$(pathToLog)log_bestfits.txt", "a") do io
                writedlm(io, sol[end].u')
            end
            open("$(pathToLog)log_fval.txt", "a") do io
                writedlm(io, sol[end].objective)    
            end
            open("$(pathToLog)log_iniparam.txt", "a") do io
                writedlm(io, u0')    
            end
        end
    else
        for row in eachrow(initialparameters)
            u0 = Vector(row)
            optprob = OptimizationFunction(err_func, ADmethod);
            prob = OptimizationProblem(optprob, u0, optim_struct_wrapper, lb = lb, ub = ub)
            # Import a solver package and solve the optimization problem
            push!(sol, solve(prob, Method; maxtime = maxtime, maxiters = maxiters));
            open("$(pathToLog)log_bestfits.txt", "a") do io
                writedlm(io, sol[end].u')
            end
            open("$(pathToLog)log_fval.txt", "a") do io
                writedlm(io, sol[end].objective)    
            end
            open("$(pathToLog)log_iniparam.txt", "a") do io
                writedlm(io, u0')    
            end
        end
    end
    return sol
end


function ini_data(optim_struct::OptimStruct, FRange; kwargs...)
    @unpack data = optim_struct
    datafit = ()
    for i in eachindex(FRange)
        if (FRange[i][2]>0) & !(data.datatypes[i] == StoThyLiveCell.Distribution_RNA()) & !(data.datatypes[i] == StoThyLiveCell.ConvolutedDistribution_RNA())
            llimit = findfirst(data.data[i][1] .>=FRange[i][1])
            ulimit = findlast(data.data[i][1] .<=FRange[i][2])
            datatemp = (data.data[i][1][llimit:ulimit], data.data[i][2][llimit:ulimit],) 
            if minimum(datatemp[2])<=0
                @warn "There are survival probabilities equal or below zero"
            end
        elseif data.datatypes[i] == StoThyLiveCell.Distribution_RNA()
            datatemp = data.data[i][(data.data[i] .>=FRange[i][1]) .& (data.data[i] .<=FRange[i][2])]
        elseif data.datatypes[i] == StoThyLiveCell.ConvolutedDistribution_RNA()
            datatemp = data.data[i][(data.data[i] .>=FRange[i][1]) .& (data.data[i] .<=FRange[i][2])]
        else
            datatemp = data.data[i]
        end
        datafit = (datafit..., datatemp)
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


function ini_optim_multipleDistribution(optim_struct::OptimStruct, datagroup::FixedCellData; kwargs...)
    function err_func(params,optim_struct_wrapper::OptimStructWrapper)
        parameters = utiles.mergeparameter_base(optim_struct_wrapper.fixedparam, params, optim_struct_wrapper.freeparametersidx)
        estimate_signal = parameters[end-1]*parameters[end].*utiles.convolution(distribution_mrna(optim_struct_wrapper.model, parameters[1:4], optim_struct_wrapper.maxrnaFC)) .+ parameters[end-1]*(1-parameters[end]).*utiles.convolution(distribution_mrna(optim_struct_wrapper.model, vcat(parameters[5:7],parameters[4]), optim_struct_wrapper.maxrnaFC))
        estimate_signal[1] = estimate_signal[1] + (1-parameters[end-1])
        return optim_struct_wrapper.dist[1](estimate_signal, optim_struct_wrapper.data_fit.data[1])    
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

function ini_optimAD(optim_struct::OptimStruct, datagroup::LiveCellData; kwargs...)
    function err_func(params::AbstractVector{T},optim_struct_wrapper::OptimStructWrapper) where T
        @unpack utileMat = optim_struct_wrapper
        parameters = utiles.mergeparameter_base(optim_struct_wrapper.fixedparam, params, optim_struct_wrapper.freeparametersidx)
        (P, ssp, weightsTr_off, PabsOff, sspTr_Off, Pabs) = StoThyLiveCell.mo_basics!(optim_struct_wrapper.model, parameters, optim_struct_wrapper.maxrnaLC, optim_struct_wrapper.utileMat.stateTr_on, optim_struct_wrapper.utileMat.stateAbs_on)
        utileMat_tmp = (stateTr=utileMat.stateTr, stateTr_on=utileMat.stateTr_on, stateAbs_on=utileMat.stateAbs_on, weightsTr_off=weightsTr_off, P=P, ssp=ssp, PabsOff=PabsOff, sspTr_Off=sspTr_Off, Pabs=Pabs)
        optim_struct_wrapper_tmp = OptimStructWrapper{typeof(optim_struct_wrapper.data_fit),typeof(optim_struct_wrapper.dist), typeof(optim_struct_wrapper.model),typeof(optim_struct_wrapper.err_func), typeof(utileMat_tmp)}(optim_struct_wrapper.data_ref, optim_struct_wrapper.data_fit, optim_struct_wrapper.dist, optim_struct_wrapper.model, optim_struct_wrapper.SRange, optim_struct_wrapper.maxrnaLC, optim_struct_wrapper.maxrnaFC, optim_struct_wrapper.freeparametersidx, optim_struct_wrapper.fixedparam, utileMat_tmp, optim_struct_wrapper.err_func)
        error = 0.  
        for i in eachindex(optim_struct_wrapper_tmp.data_fit.datatypes)
            estimate_signal = optim_struct_wrapper_tmp.data_fit.datatypes[i](i, optim_struct_wrapper_tmp)
            error = error + optim_struct_wrapper_tmp.dist[i](estimate_signal, optim_struct_wrapper_tmp.data_fit.data[i])
        end
        return error
    end
    return err_func
end

function ini_optimAD(optim_struct::OptimStruct, datagroup::FixedCellData; kwargs...)
    function err_func(params::AbstractVector{T},optim_struct_wrapper::OptimStructWrapper) where T
        parameters = utiles.mergeparameter_base(optim_struct_wrapper.fixedparam, params, optim_struct_wrapper.freeparametersidx)
        Qrna = StoThyLiveCell.distrna_basic(optim_struct_wrapper.model, parameters, optim_struct_wrapper.maxrnaFC) 
        utileMat_tmp = (Qrna = Qrna,)
        optim_struct_wrapper_tmp = OptimStructWrapper{typeof(optim_struct_wrapper.data_fit),typeof(optim_struct_wrapper.dist), typeof(optim_struct_wrapper.model),typeof(optim_struct_wrapper.err_func), typeof(utileMat_tmp)}(optim_struct_wrapper.data_ref, optim_struct_wrapper.data_fit, optim_struct_wrapper.dist, optim_struct_wrapper.model, optim_struct_wrapper.SRange, optim_struct_wrapper.maxrnaLC, optim_struct_wrapper.maxrnaFC, optim_struct_wrapper.freeparametersidx, optim_struct_wrapper.fixedparam, utileMat_tmp, optim_struct_wrapper.err_func)
        estimate_signal = optim_struct_wrapper_tmp.data_fit.datatypes[1](1, optim_struct_wrapper_tmp)
        error = optim_struct_wrapper_tmp.dist[1](estimate_signal,optim_struct_wrapper_tmp.data_fit.data[1])
        return error
    end
    return err_func
end

function ini_optimAD(optim_struct::OptimStruct, datagroup::FixedAndLiveCellData; kwargs...)
    @warn "A mixture of live cell and fixed cell data is used in the error function"
    @warn " The  last parameter is interpreted as the degradation rate in the calculations of the mRNA number distribution"
    function err_func(params::AbstractVector{T},optim_struct_wrapper::OptimStructWrapper) where T
        @unpack utileMat = optim_struct_wrapper
        parameters = utiles.mergeparameter_base(optim_struct_wrapper.fixedparam, params, optim_struct_wrapper.freeparametersidx)
        (P, ssp, weightsTr_off, PabsOff, sspTr_Off, Pabs) = StoThyLiveCell.mo_basics!(optim_struct_wrapper.model, parameters, optim_struct_wrapper.maxrnaLC, optim_struct_wrapper.utileMat.stateTr_on, optim_struct_wrapper.utileMat.stateAbs_on)
        Qrna = StoThyLiveCell.distrna_basic(optim_struct_wrapper.model, parameters, optim_struct_wrapper.maxrnaFC) 
        utileMat_tmp = (stateTr=utileMat.stateTr, stateTr_on=utileMat.stateTr_on, stateAbs_on=utileMat.stateAbs_on, weightsTr_off=weightsTr_off, P=P, ssp=ssp, PabsOff=PabsOff, sspTr_Off=sspTr_Off, Pabs=Pabs, Qrna=Qrna)
        optim_struct_wrapper_tmp = OptimStructWrapper{typeof(optim_struct_wrapper.data_fit),typeof(optim_struct_wrapper.dist), typeof(optim_struct_wrapper.model),typeof(optim_struct_wrapper.err_func), typeof(utileMat_tmp)}(optim_struct_wrapper.data_ref, optim_struct_wrapper.data_fit, optim_struct_wrapper.dist, optim_struct_wrapper.model, optim_struct_wrapper.SRange, optim_struct_wrapper.maxrnaLC, optim_struct_wrapper.maxrnaFC, optim_struct_wrapper.freeparametersidx, optim_struct_wrapper.fixedparam, utileMat_tmp, optim_struct_wrapper.err_func)
        
        error = 0.  
        for i in eachindex(optim_struct_wrapper_tmp.data_fit.datatypes)
            estimate_signal = optim_struct_wrapper_tmp.data_fit.datatypes[i](i, optim_struct_wrapper_tmp)
            error = error + optim_struct_wrapper_tmp.dist[i](estimate_signal,optim_struct_wrapper_tmp.data_fit.data[i])
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
 
function ini_optim_withoutsingletAD(optim_struct::OptimStruct, datagroup::LiveCellData; kwargs...)
    function err_func(params::AbstractVector{T},optim_struct_wrapper::OptimStructWrapper) where T
        @unpack utileMat = optim_struct_wrapper
        parameters = utiles.mergeparameter_base(optim_struct_wrapper.fixedparam, params, optim_struct_wrapper.freeparametersidx)
        (P, ssp, Pwos,  weightsAbs_off_wos, sspTr_off_wos, weightsAbs_on, sspPreB, weightsTr_on, PabsOff, weightsTr_on_wos, weightsAbsorbed_off_wos, sspwos, weightsPre_on_and_on, Rn, NR, Nc, Qn, Nn, weightsTr_off_wos, Pabs_wos, weightsON_wos, weightsPre_on_wos) = StoThyLiveCell.mo_basics_wosingletAD(optim_struct_wrapper.model, parameters, optim_struct_wrapper.maxrnaLC, utileMat.stateTr_on, utileMat.stateAbs_on, utileMat.totnbs, utileMat.stateAbs_on_wos, utileMat.statePre_on_wos)
        utileMat_tmp = (nascentbin=utileMat.nascentbin, P=P, ssp=ssp, stateTr=utileMat.stateTr, stateTr_on=utileMat.stateTr_on, stateAbs_on=utileMat.stateAbs_on, totnbs=utileMat.totnbs, Pwos=Pwos, stateAbs_on_wos=utileMat.stateAbs_on_wos, statePre_on_wos=utileMat.statePre_on_wos, weightsAbs_off_wos=weightsAbs_off_wos, sspTr_off_wos=sspTr_off_wos, weightsAbs_on=weightsAbs_on, sspPreB=sspPreB, weightsTr_on=weightsTr_on, PabsOff=PabsOff, weightsTr_on_wos=weightsTr_on_wos, weightsAbsorbed_off_wos=weightsAbsorbed_off_wos, sspwos=sspwos, weightsPre_on_and_on=weightsPre_on_and_on, Rn=Rn, NR=NR, Nc=Nc, Qn=Qn, Nn=Nn, weightsTr_off_wos=weightsTr_off_wos, Pabs_wos=Pabs_wos, weightsON_wos=weightsON_wos, rnanbvec_on=utileMat.rnanbvec_on, weightsPre_on_wos=weightsPre_on_wos)
        optim_struct_wrapper_tmp = OptimStructWrapper{typeof(optim_struct_wrapper.data_fit),typeof(optim_struct_wrapper.dist), typeof(optim_struct_wrapper.model),typeof(optim_struct_wrapper.err_func), typeof(utileMat_tmp)}(optim_struct_wrapper.data_ref, optim_struct_wrapper.data_fit, optim_struct_wrapper.dist, optim_struct_wrapper.model, optim_struct_wrapper.SRange, optim_struct_wrapper.maxrnaLC, optim_struct_wrapper.maxrnaFC, optim_struct_wrapper.freeparametersidx, optim_struct_wrapper.fixedparam, utileMat_tmp, optim_struct_wrapper.err_func)
        error = 0.  
        for i in eachindex(optim_struct_wrapper_tmp.data_fit.datatypes)
            estimate_signal = optim_struct_wrapper_tmp.data_fit.datatypes[i](i, optim_struct_wrapper_tmp, :without)
            error = error + optim_struct_wrapper_tmp.dist[i](estimate_signal, optim_struct_wrapper_tmp.data_fit.data[i])
        end
        return error
    end
    return err_func
end

function ini_optim_withoutsingletAD(optim_struct::OptimStruct, datagroup::FixedAndLiveCellData; kwargs...)
    @warn "A mixture of live cell and fixed cell data is used in the error function"
    @warn " The  last parameter is interpreted as the degradation rate in the calculations of the mRNA number distribution"
    function err_func(params::AbstractVector{T},optim_struct_wrapper::OptimStructWrapper) where T
        @unpack utileMat = optim_struct_wrapper
        parameters = utiles.mergeparameter_base(optim_struct_wrapper.fixedparam, params, optim_struct_wrapper.freeparametersidx)
        (P, ssp, Pwos,  weightsAbs_off_wos, sspTr_off_wos, weightsAbs_on, sspPreB, weightsTr_on, PabsOff, weightsTr_on_wos, weightsAbsorbed_off_wos, sspwos, weightsPre_on_and_on, Rn, NR, Nc, Qn, Nn, weightsTr_off_wos, Pabs_wos, weightsON_wos, weightsPre_on_wos) = StoThyLiveCell.mo_basics_wosingletAD(optim_struct_wrapper.model, parameters, optim_struct_wrapper.maxrnaLC, utileMat.stateTr_on, utileMat.stateAbs_on, utileMat.totnbs, utileMat.stateAbs_on_wos, utileMat.statePre_on_wos)
        Qrna = StoThyLiveCell.distrna_basic(optim_struct_wrapper.model, parameters, optim_struct_wrapper.maxrnaFC) 
        utileMat_tmp = (nascentbin=utileMat.nascentbin, P=P, ssp=ssp, stateTr=utileMat.stateTr, stateTr_on=utileMat.stateTr_on, stateAbs_on=utileMat.stateAbs_on, totnbs=utileMat.totnbs, Pwos=Pwos, stateAbs_on_wos=utileMat.stateAbs_on_wos, statePre_on_wos=utileMat.statePre_on_wos, weightsAbs_off_wos=weightsAbs_off_wos, sspTr_off_wos=sspTr_off_wos, weightsAbs_on=weightsAbs_on, sspPreB=sspPreB, weightsTr_on=weightsTr_on, PabsOff=PabsOff, weightsTr_on_wos=weightsTr_on_wos, weightsAbsorbed_off_wos=weightsAbsorbed_off_wos, sspwos=sspwos, weightsPre_on_and_on=weightsPre_on_and_on, Rn=Rn, NR=NR, Nc=Nc, Qn=Qn, Nn=Nn, weightsTr_off_wos=weightsTr_off_wos, Pabs_wos=Pabs_wos, weightsON_wos=weightsON_wos, rnanbvec_on=utileMat.rnanbvec_on, weightsPre_on_wos=weightsPre_on_wos, Qrna=Qrna)
        optim_struct_wrapper_tmp = OptimStructWrapper{typeof(optim_struct_wrapper.data_fit),typeof(optim_struct_wrapper.dist), typeof(optim_struct_wrapper.model),typeof(optim_struct_wrapper.err_func), typeof(utileMat_tmp)}(optim_struct_wrapper.data_ref, optim_struct_wrapper.data_fit, optim_struct_wrapper.dist, optim_struct_wrapper.model, optim_struct_wrapper.SRange, optim_struct_wrapper.maxrnaLC, optim_struct_wrapper.maxrnaFC, optim_struct_wrapper.freeparametersidx, optim_struct_wrapper.fixedparam, utileMat_tmp, optim_struct_wrapper.err_func)
        error = 0.  
        for i in eachindex(optim_struct_wrapper_tmp.data_fit.datatypes)
            estimate_signal = optim_struct_wrapper_tmp.data_fit.datatypes[i](i, optim_struct_wrapper_tmp, :without)
            error = error + optim_struct_wrapper_tmp.dist[i](estimate_signal,optim_struct_wrapper_tmp.data_fit.data[i])
        end
        return error
    end 
    return err_func
end

#initiation with reparameterization

function ini_optim(optim_struct::OptimStruct, datagroup::LiveCellData, f_parameter; kwargs...)
    function err_func(params,optim_struct_wrapper::OptimStructWrapper)
        parameters = f_parameter(utiles.mergeparameter_base(optim_struct_wrapper.fixedparam, params, optim_struct_wrapper.freeparametersidx))
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

function ini_optim(optim_struct::OptimStruct, datagroup::FixedCellData, f_parameter; kwargs...)
    function err_func(params,optim_struct_wrapper::OptimStructWrapper)
        parameters = f_parameter(utiles.mergeparameter_base(optim_struct_wrapper.fixedparam, params, optim_struct_wrapper.freeparametersidx))
        #@unpack utileMat = optim_struct_wrapper
        #model outputs
        StoThyLiveCell.distrna_basic!(optim_struct_wrapper.model, parameters, optim_struct_wrapper.maxrnaFC, optim_struct_wrapper.utileMat.Qrna) 
        estimate_signal = optim_struct_wrapper.data_fit.datatypes[1](1,optim_struct_wrapper)
        error = optim_struct_wrapper.dist[1](estimate_signal,optim_struct_wrapper.data_fit.data[1])
        return error
    end
    return err_func
end

function ini_optim(optim_struct::OptimStruct, datagroup::FixedAndLiveCellData, f_parameter; kwargs...)
    @warn "A mixture of live cell and fixed cell data is used in the error function"
    @warn " The  last parameter is interpreted as the degradation rate in the calculations of the mRNA number distribution"
    function err_func(params,optim_struct_wrapper::OptimStructWrapper)
        parameters = f_parameter(utiles.mergeparameter_base(optim_struct_wrapper.fixedparam, params, optim_struct_wrapper.freeparametersidx))
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

function ini_optimAD(optim_struct::OptimStruct, datagroup::LiveCellData, f_parameter; kwargs...)
    function err_func(params::AbstractVector{T},optim_struct_wrapper::OptimStructWrapper) where T
        @unpack utileMat = optim_struct_wrapper
        parameters = f_parameter(utiles.mergeparameter_base(optim_struct_wrapper.fixedparam, params, optim_struct_wrapper.freeparametersidx))
        (P, ssp, weightsTr_off, PabsOff, sspTr_Off, Pabs) = StoThyLiveCell.mo_basics!(optim_struct_wrapper.model, parameters, optim_struct_wrapper.maxrnaLC, optim_struct_wrapper.utileMat.stateTr_on, optim_struct_wrapper.utileMat.stateAbs_on)
        utileMat_tmp = (stateTr=utileMat.stateTr, stateTr_on=utileMat.stateTr_on, stateAbs_on=utileMat.stateAbs_on, weightsTr_off=weightsTr_off, P=P, ssp=ssp, PabsOff=PabsOff, sspTr_Off=sspTr_Off, Pabs=Pabs)
        optim_struct_wrapper_tmp = OptimStructWrapper{typeof(optim_struct_wrapper.data_fit),typeof(optim_struct_wrapper.dist), typeof(optim_struct_wrapper.model),typeof(optim_struct_wrapper.err_func), typeof(utileMat_tmp)}(optim_struct_wrapper.data_ref, optim_struct_wrapper.data_fit, optim_struct_wrapper.dist, optim_struct_wrapper.model, optim_struct_wrapper.SRange, optim_struct_wrapper.maxrnaLC, optim_struct_wrapper.maxrnaFC, optim_struct_wrapper.freeparametersidx, optim_struct_wrapper.fixedparam, utileMat_tmp, optim_struct_wrapper.err_func)
        error = 0.  
        for i in eachindex(optim_struct_wrapper_tmp.data_fit.datatypes)
            estimate_signal = optim_struct_wrapper_tmp.data_fit.datatypes[i](i, optim_struct_wrapper_tmp)
            error = error + optim_struct_wrapper_tmp.dist[i](estimate_signal, optim_struct_wrapper_tmp.data_fit.data[i])
        end
        return error
    end
    return err_func
end

function ini_optimAD(optim_struct::OptimStruct, datagroup::FixedCellData, f_parameter; kwargs...)
    function err_func(params::AbstractVector{T},optim_struct_wrapper::OptimStructWrapper) where T
        parameters = f_parameter(utiles.mergeparameter_base(optim_struct_wrapper.fixedparam, params, optim_struct_wrapper.freeparametersidx))
        Qrna = StoThyLiveCell.distrna_basic(optim_struct_wrapper.model, parameters, optim_struct_wrapper.maxrnaFC) 
        utileMat_tmp = (Qrna = Qrna,)
        optim_struct_wrapper_tmp = OptimStructWrapper{typeof(optim_struct_wrapper.data_fit),typeof(optim_struct_wrapper.dist), typeof(optim_struct_wrapper.model),typeof(optim_struct_wrapper.err_func), typeof(utileMat_tmp)}(optim_struct_wrapper.data_ref, optim_struct_wrapper.data_fit, optim_struct_wrapper.dist, optim_struct_wrapper.model, optim_struct_wrapper.SRange, optim_struct_wrapper.maxrnaLC, optim_struct_wrapper.maxrnaFC, optim_struct_wrapper.freeparametersidx, optim_struct_wrapper.fixedparam, utileMat_tmp, optim_struct_wrapper.err_func)
        estimate_signal = optim_struct_wrapper_tmp.data_fit.datatypes[1](1, optim_struct_wrapper_tmp)
        error = optim_struct_wrapper_tmp.dist[1](estimate_signal,optim_struct_wrapper_tmp.data_fit.data[1])
        return error
    end
    return err_func
end

function ini_optimAD(optim_struct::OptimStruct, datagroup::FixedAndLiveCellData, f_parameter; kwargs...)
    @warn "A mixture of live cell and fixed cell data is used in the error function"
    @warn " The  last parameter is interpreted as the degradation rate in the calculations of the mRNA number distribution"
    function err_func(params::AbstractVector{T},optim_struct_wrapper::OptimStructWrapper) where T
        @unpack utileMat = optim_struct_wrapper
        parameters = f_parameter(utiles.mergeparameter_base(optim_struct_wrapper.fixedparam, params, optim_struct_wrapper.freeparametersidx))
        (P, ssp, weightsTr_off, PabsOff, sspTr_Off, Pabs) = StoThyLiveCell.mo_basics!(optim_struct_wrapper.model, parameters, optim_struct_wrapper.maxrnaLC, optim_struct_wrapper.utileMat.stateTr_on, optim_struct_wrapper.utileMat.stateAbs_on)
        Qrna = StoThyLiveCell.distrna_basic(optim_struct_wrapper.model, parameters, optim_struct_wrapper.maxrnaFC) 
        utileMat_tmp = (stateTr=utileMat.stateTr, stateTr_on=utileMat.stateTr_on, stateAbs_on=utileMat.stateAbs_on, weightsTr_off=weightsTr_off, P=P, ssp=ssp, PabsOff=PabsOff, sspTr_Off=sspTr_Off, Pabs=Pabs, Qrna=Qrna)
        optim_struct_wrapper_tmp = OptimStructWrapper{typeof(optim_struct_wrapper.data_fit),typeof(optim_struct_wrapper.dist), typeof(optim_struct_wrapper.model),typeof(optim_struct_wrapper.err_func), typeof(utileMat_tmp)}(optim_struct_wrapper.data_ref, optim_struct_wrapper.data_fit, optim_struct_wrapper.dist, optim_struct_wrapper.model, optim_struct_wrapper.SRange, optim_struct_wrapper.maxrnaLC, optim_struct_wrapper.maxrnaFC, optim_struct_wrapper.freeparametersidx, optim_struct_wrapper.fixedparam, utileMat_tmp, optim_struct_wrapper.err_func)
        
        error = 0.  
        for i in eachindex(optim_struct_wrapper_tmp.data_fit.datatypes)
            estimate_signal = optim_struct_wrapper_tmp.data_fit.datatypes[i](i, optim_struct_wrapper_tmp)
            error = error + optim_struct_wrapper_tmp.dist[i](estimate_signal,optim_struct_wrapper_tmp.data_fit.data[i])
        end
        return error
    end 
    return err_func
end

function ini_optim_withoutsinglet(optim_struct::OptimStruct, datagroup::LiveCellData, f_parameter; kwargs...)
    function err_func(params,optim_struct_wrapper::OptimStructWrapper)
        parameters = f_parameter(utiles.mergeparameter_base(optim_struct_wrapper.fixedparam, params, optim_struct_wrapper.freeparametersidx))
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

function ini_optim_withoutsinglet(optim_struct::OptimStruct, datagroup::FixedAndLiveCellData, f_parameter; kwargs...)
    @warn "A mixture of live cell and fixed cell data is used in the error function"
    @warn " The  last parameter is interpreted as the degradation rate in the calculations of the mRNA number distribution"
    function err_func(params,optim_struct_wrapper::OptimStructWrapper)
        parameters = f_parameter(utiles.mergeparameter_base(optim_struct_wrapper.fixedparam, params, optim_struct_wrapper.freeparametersidx))
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
 
function ini_optim_withoutsingletAD(optim_struct::OptimStruct, datagroup::LiveCellData, f_parameter; kwargs...)
    function err_func(params::AbstractVector{T},optim_struct_wrapper::OptimStructWrapper) where T
        @unpack utileMat = optim_struct_wrapper
        parameters =f_parameter(utiles.mergeparameter_base(optim_struct_wrapper.fixedparam, params, optim_struct_wrapper.freeparametersidx))
        (P, ssp, Pwos,  weightsAbs_off_wos, sspTr_off_wos, weightsAbs_on, sspPreB, weightsTr_on, PabsOff, weightsTr_on_wos, weightsAbsorbed_off_wos, sspwos, weightsPre_on_and_on, Rn, NR, Nc, Qn, Nn, weightsTr_off_wos, Pabs_wos, weightsON_wos, weightsPre_on_wos) = StoThyLiveCell.mo_basics_wosingletAD(optim_struct_wrapper.model, parameters, optim_struct_wrapper.maxrnaLC, utileMat.stateTr_on, utileMat.stateAbs_on, utileMat.totnbs, utileMat.stateAbs_on_wos, utileMat.statePre_on_wos)
        utileMat_tmp = (nascentbin=utileMat.nascentbin, P=P, ssp=ssp, stateTr=utileMat.stateTr, stateTr_on=utileMat.stateTr_on, stateAbs_on=utileMat.stateAbs_on, totnbs=utileMat.totnbs, Pwos=Pwos, stateAbs_on_wos=utileMat.stateAbs_on_wos, statePre_on_wos=utileMat.statePre_on_wos, weightsAbs_off_wos=weightsAbs_off_wos, sspTr_off_wos=sspTr_off_wos, weightsAbs_on=weightsAbs_on, sspPreB=sspPreB, weightsTr_on=weightsTr_on, PabsOff=PabsOff, weightsTr_on_wos=weightsTr_on_wos, weightsAbsorbed_off_wos=weightsAbsorbed_off_wos, sspwos=sspwos, weightsPre_on_and_on=weightsPre_on_and_on, Rn=Rn, NR=NR, Nc=Nc, Qn=Qn, Nn=Nn, weightsTr_off_wos=weightsTr_off_wos, Pabs_wos=Pabs_wos, weightsON_wos=weightsON_wos, rnanbvec_on=utileMat.rnanbvec_on, weightsPre_on_wos=weightsPre_on_wos)
        optim_struct_wrapper_tmp = OptimStructWrapper{typeof(optim_struct_wrapper.data_fit),typeof(optim_struct_wrapper.dist), typeof(optim_struct_wrapper.model),typeof(optim_struct_wrapper.err_func), typeof(utileMat_tmp)}(optim_struct_wrapper.data_ref, optim_struct_wrapper.data_fit, optim_struct_wrapper.dist, optim_struct_wrapper.model, optim_struct_wrapper.SRange, optim_struct_wrapper.maxrnaLC, optim_struct_wrapper.maxrnaFC, optim_struct_wrapper.freeparametersidx, optim_struct_wrapper.fixedparam, utileMat_tmp, optim_struct_wrapper.err_func)
        error = 0.  
        for i in eachindex(optim_struct_wrapper_tmp.data_fit.datatypes)
            estimate_signal = optim_struct_wrapper_tmp.data_fit.datatypes[i](i, optim_struct_wrapper_tmp, :without)
            error = error + optim_struct_wrapper_tmp.dist[i](estimate_signal, optim_struct_wrapper_tmp.data_fit.data[i])
        end
        return error
    end
    return err_func
end

function ini_optim_withoutsingletAD(optim_struct::OptimStruct, datagroup::FixedAndLiveCellData, f_parameter; kwargs...)
    @warn "A mixture of live cell and fixed cell data is used in the error function"
    @warn " The  last parameter is interpreted as the degradation rate in the calculations of the mRNA number distribution"
    function err_func(params::AbstractVector{T},optim_struct_wrapper::OptimStructWrapper) where T
        @unpack utileMat = optim_struct_wrapper
        parameters = f_parameter(utiles.mergeparameter_base(optim_struct_wrapper.fixedparam, params, optim_struct_wrapper.freeparametersidx))
        (P, ssp, Pwos,  weightsAbs_off_wos, sspTr_off_wos, weightsAbs_on, sspPreB, weightsTr_on, PabsOff, weightsTr_on_wos, weightsAbsorbed_off_wos, sspwos, weightsPre_on_and_on, Rn, NR, Nc, Qn, Nn, weightsTr_off_wos, Pabs_wos, weightsON_wos, weightsPre_on_wos) = StoThyLiveCell.mo_basics_wosingletAD(optim_struct_wrapper.model, parameters, optim_struct_wrapper.maxrnaLC, utileMat.stateTr_on, utileMat.stateAbs_on, utileMat.totnbs, utileMat.stateAbs_on_wos, utileMat.statePre_on_wos)
        Qrna = StoThyLiveCell.distrna_basic(optim_struct_wrapper.model, parameters, optim_struct_wrapper.maxrnaFC) 
        utileMat_tmp = (nascentbin=utileMat.nascentbin, P=P, ssp=ssp, stateTr=utileMat.stateTr, stateTr_on=utileMat.stateTr_on, stateAbs_on=utileMat.stateAbs_on, totnbs=utileMat.totnbs, Pwos=Pwos, stateAbs_on_wos=utileMat.stateAbs_on_wos, statePre_on_wos=utileMat.statePre_on_wos, weightsAbs_off_wos=weightsAbs_off_wos, sspTr_off_wos=sspTr_off_wos, weightsAbs_on=weightsAbs_on, sspPreB=sspPreB, weightsTr_on=weightsTr_on, PabsOff=PabsOff, weightsTr_on_wos=weightsTr_on_wos, weightsAbsorbed_off_wos=weightsAbsorbed_off_wos, sspwos=sspwos, weightsPre_on_and_on=weightsPre_on_and_on, Rn=Rn, NR=NR, Nc=Nc, Qn=Qn, Nn=Nn, weightsTr_off_wos=weightsTr_off_wos, Pabs_wos=Pabs_wos, weightsON_wos=weightsON_wos, rnanbvec_on=utileMat.rnanbvec_on, weightsPre_on_wos=weightsPre_on_wos, Qrna=Qrna)
        optim_struct_wrapper_tmp = OptimStructWrapper{typeof(optim_struct_wrapper.data_fit),typeof(optim_struct_wrapper.dist), typeof(optim_struct_wrapper.model),typeof(optim_struct_wrapper.err_func), typeof(utileMat_tmp)}(optim_struct_wrapper.data_ref, optim_struct_wrapper.data_fit, optim_struct_wrapper.dist, optim_struct_wrapper.model, optim_struct_wrapper.SRange, optim_struct_wrapper.maxrnaLC, optim_struct_wrapper.maxrnaFC, optim_struct_wrapper.freeparametersidx, optim_struct_wrapper.fixedparam, utileMat_tmp, optim_struct_wrapper.err_func)
        error = 0.  
        for i in eachindex(optim_struct_wrapper_tmp.data_fit.datatypes)
            estimate_signal = optim_struct_wrapper_tmp.data_fit.datatypes[i](i, optim_struct_wrapper_tmp, :without)
            error = error + optim_struct_wrapper_tmp.dist[i](estimate_signal,optim_struct_wrapper_tmp.data_fit.data[i])
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
#= 
function (f::Survival_Burst)(dataidx::Int,  parameters::AbstractVector{T},optimstruct::OptimStructWrapper)
    @unpack utileMat, data_fit,model, maxrnaLC = optimstruct
    @unpack stateTr_on, stateAbs_on = utileMat

    return survival_burst(model, parameters, maxrnaLC, stateTr_on, stateAbs_on, data_fit.data[dataidx][1]) 
   
end =#

function (f::Survival_InterBurst)(dataidx::Int, optimstruct::OptimStructWrapper)
    @unpack utileMat, data_fit = optimstruct
    @unpack weightsTr_off, PabsOff = utileMat
    return survival_interburst(PabsOff, weightsTr_off,data_fit.data[dataidx][1])
end
#= 
function (f::Survival_InterBurst)(dataidx::Int, parameters::AbstractVector{T}, optimstruct::OptimStructWrapper)
    @unpack data_fit,model, maxrnaLC = optimstruct
    return survival_interburst(model, parameters, maxrnaLC,data_fit.data[dataidx][1])
end =#

function (f::Survival_NextBurst)(dataidx::Int, optimstruct::OptimStructWrapper)
    @unpack utileMat, data_fit = optimstruct
    @unpack sspTr_off, PabsOff = utileMat

    return survival_nextburst(sspTr_off::Vector{Float64},PabsOff,data_fit.data[dataidx][1])
end
#= 
function (f::Survival_NextBurst)(dataidx::Int,parameters::AbstractVector{T}, optimstruct::OptimStructWrapper)
    @unpack data_fit, model, maxrnaLC = optimstruct

    return survival_nextburst(model, parameters, maxrnaLC, data_fit.data[dataidx][1])
end =#

function (f::Mean_Nascent)(dataidx::Int, optimstruct::OptimStructWrapper)
    @unpack utileMat, data_fit = optimstruct
    @unpack stateTr, ssp = utileMat

    return mean_nascentrna(ssp, optimstruct.maxrnaLC, stateTr, optimstruct.model.nbstate, data_fit.detectionLimitNS)
end
#= 
function (f::Mean_Nascent)(dataidx::Int,parameters::AbstractVector{T}, optimstruct::OptimStructWrapper)
    @unpack utileMat, data_fit, model, maxrnaLC = optimstruct
    @unpack stateTr = utileMat

    return mean_nascentrna(model, parameters, maxrnaLC, stateTr, data_fit.detectionLimitNS)
end =#


function (f::Prob_Burst)(dataidx::Int, optimstruct::OptimStructWrapper)
    @unpack utileMat = optimstruct
    @unpack stateTr_on, ssp = utileMat

    return prob_burst(ssp,stateTr_on) 
end
#= 
function (f::Prob_Burst)(dataidx::Int, parameters::AbstractVector{T}, optimstruct::OptimStructWrapper)
    @unpack utileMat, model, maxrnaLC = optimstruct
    @unpack stateTr_on = utileMat

    return prob_burst(model, maxrnaLC, stateTr_on) 
end
 =#
function (f::Intensity_Burst)(dataidx::Int, optimstruct::OptimStructWrapper)
    @unpack utileMat, data_fit, model, maxrnaLC = optimstruct
    @unpack stateTr_on, stateAbs_on, Pabs, P, sspTr_off = utileMat

    return intensity_burst(data_fit.detectionLimitLC, P, Pabs, sspTr_off,stateTr_on, stateAbs_on,data_fit.data[dataidx][1], model.nbstate, maxrnaLC)
end
#= 
function (f::Intensity_Burst)(dataidx::Int,parameters::AbstractVector{T}, optimstruct::OptimStructWrapper)
    @unpack utileMat, data_fit, model, maxrnaLC = optimstruct
    @unpack stateTr_on, stateAbs_on = utileMat

    return intensity_burst(model, parameters, data_fit.detectionLimitLC,stateTr_on, stateAbs_on, data_fit.data[dataidx][1], maxrnaLC)
end =#

function (f::Correlation_InterBurst)(dataidx::Int, optimstruct::OptimStructWrapper)
    @unpack utileMat = optimstruct
    
    @unpack stateTr_on, stateAbs_on, weightsTr_off, P = utileMat

  return correlation_interburst(P, weightsTr_off,stateAbs_on, stateTr_on, 15000)
end
#= 
function (f::Correlation_InterBurst)(dataidx::Int, parameters::AbstractVector{T}, optimstruct::OptimStructWrapper)
    @unpack utileMat, model, maxrnaLC = optimstruct
    @unpack stateTr_on, stateAbs_on = utileMat

  return correlation_interburst(model, parameters, maxrnaLC, stateAbs_on, stateTr_on, 15000)
end =#


function (f::Distribution_RNA)(dataidx::Int, optimstruct::OptimStructWrapper)
    @unpack utileMat, maxrnaFC, model = optimstruct
    @unpack Qrna = utileMat
    return distribution_mrna(Qrna, maxrnaFC, model.nbstate)
end

#= 
function (f::Distribution_RNA)(dataidx::Int, parameters::AbstractVector{T}, optimstruct::OptimStructWrapper) where T
    @unpack maxrnaFC, model = optimstruct
    return distribution_mrna(model,parameters,maxrnaFC,)
end
 =#

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