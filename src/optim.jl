#functions used in the optimization procedure

struct OptimStruct{DT,D,DI,M}
    data::DataFit{DT,D}
    dist::Vector{DI}
    model::M
end

function OptimStruct(data::Vector, dist::D, model::M) where {D,M}
    return OptimStruct{typeof(dist),typeof(model)}(data, dist, model)
end

struct OptimStructWrapper{DT,D,DI,M,SR,EF}
    data::DataFit{DT,D}
    FRange::Vector{Int}
    dist::Vector{DI}
    model::M
    SRange::SR
    freeparametersidx::Vector{Int}
    fixedparam::Vector{Float32}
    utileMat::NamedTuple{(:stateTr, :stateTr_on, :stateAbs_on, :weightsTr_off, :P, :ssp, :PabsOff),Tuple{Vector{Int64}, Vector{Int64}, Vector{Int64}, Vector{Float64}, Array{Float64,2}, Vector{Float64}, Array{Float64,2}}}
    err_func::EF
end


function optim_function(SRange, optim_struct::OptimStruct, args...; kwargs...)

    err_func = ini_optim(optim_struct; kwargs...)

    optim_struct_wrapper = OptimStructWrapper(optim_struct.data,FRange, optim_struct.dist, optim_struct.model, SRange, freeparameteridx,fixedparam, utileMat, err_func)

    thetax = start_optim(optim_struct_wrapper, args...; kwargs...)

    model, stage = optim_struct.model, optim_struct.stage
    params = thetax[1:end-1]
    estimate_data = compute_distribution(params, NT, model, stage, infer_counts, filter=filter_uniform)
    return thetax, hcat(estimate_data, reference_data)
end


function start_optim(optim_struct_wrapper::OptimStructWrapper, args...; maxtime::Int=100, maxiters::Int=100 , Method=:adaptive_de_rand_1_bin_radiuslimited, kwargs...)
    @unpack SRange, err_func = optim_struct_wrapper
    optprob = OptimizationFunction(lsq);
    prob = OptimizationProblem(optprob, u0, metaparam, lb = lbu, ub = ubu)
    # Import a solver package and solve the optimization problem
    sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxtime = maxtime, maxiters = maxiters);
    return thetax
end

function ini_optim(optim_struct::OptimStruct; kwargs...)
    @unpack data, dist, model = optim_struct

    utileMat = StoThyLiveCell.mo_basics(model, zeros(model.nbparameters+model.nbkini+1), data.maxrnaLC, data.detectionlimitLC) 
    
    optim_struct_wrapper = OptimStructWrapper(optim_struct.data,FRange, optim_struct.dist, optim_struct.model, SRange, freeparameteridx,fixedparam, utileMat, [])

    if :burst == data.datagroups 
        function err_func(params,optim_struct_wrapper::OptimStructWrapper)
            parameters = usefullfunctions.mergeparameter_base(moptim_struct_wrapper.fixedparam, params, optim_struct_wrapper.freeparametersidx)
            @unpack utileMat = optim_struct_wrapper
            #model outputs
            StoThyLiveCell.mo_basics!(model, parameters[1:end-1], optim_struct_wrapper.maxrna, utileMat.P, utileMat.ssp, utileMat.stateTr_on, utileMat.stateAbs_on, utileMat.weightsTr_off, utileMat.PabsOff) 
            error = 0  
            for i in eachindex(optim_struct_wrapper.data.datatype)
                estimate_signal = optim_struct_wrapper.data.datatype[i](tmax,optim_struct_wrapper)
                error = error + optim_struct_wrapper.data.dist[i](estimate_signal,optim_struct_wrapper.data.data[i])
            end
        end
    end
    return err_func
end


