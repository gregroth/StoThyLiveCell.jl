using StoThyLiveCell
using Test
using DataFrames, FileIO, JLD2
using BenchmarkTools
@testset "Test 2s3r model" begin
    #define the 2s 3r model
    #(r12,r21,r23,r32,k1on,k2on,k3on,koff)
    Qstate = [0    8    4    0    0    0;
            7    0    0    4    0    0;
            3    0    0    8    2    0;
            0    3    6    0    0    2;
            0    0    1    0    0    8;
            0    0    0    1    5    0]
    paramToRate_idx = findall(Qstate .>0)
    paramToRate_val = Qstate[findall(Qstate .>0)]
    model1 = StoThyLiveCell.StandardStoModel(6,8,1,paramToRate_idx,paramToRate_val,[1,3,5],[9,9,9],10)

    #creating an instance 
    parameters = [0.0178504,  0.0436684,  0.0543096,  0.427785,  0.023986,  0.308174,  2.24418,  1.28387, 3.80846,1.]
    maxrna = 25
    detectionLimitLC = 1
    detectionLimitNS = 2

    P1 = StoModel(model1, parameters, maxrna)

    timevec = 1:1:200
    timevec_on = 1:1:10
    timevec_int = 1:1:20

    (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =  ModelOutput(model1, parameters, maxrna,detectionLimitLC, detectionLimitNS, timevec_on[end],timevec[end],timevec[end],timevec_int[end])


    @test mean_nascentrna ≈ 3.0337754651143656
    @test prob_burst ≈ 0.16568033517579464
    @test correlation_interburst ≈ 0.10866194618122918
    @test intensity_burst[1] ≈ 1
    @test intensity_burst[20] ≈ 0.0021212008496347893/1.9075373067561898
    @test survival_nextburst[1] ≈ 0.9409352782206553
    @test survival_nextburst[200] ≈  0.0020411810555407066
    @test survival_interburst[1] ≈ 0.8346193347912076
    @test survival_interburst[200] ≈  0.0010301814686152551
    @test survival_burst[1] ≈ 0.742997360615419
    @test survival_burst[10] ≈  0.024598091184660095
end

@testset "Test single output fcts" begin
    Qstate = [0    8    4    0    0    0;
                7    0    0    4    0    0;
                3    0    0    8    2    0;
                0    3    6    0    0    2;
                0    0    1    0    0    8;
                0    0    0    1    5    0]
    paramToRate_idx = findall(Qstate .>0)
    paramToRate_val = Qstate[findall(Qstate .>0)]
    model1 = StoThyLiveCell.StandardStoModel(6,8,1,paramToRate_idx,paramToRate_val,[1,3,5],[9,9,9],10)

    #creating an instance 
    parameters = [0.0178504,  0.0436684,  0.0543096,  0.427785,  0.023986,  0.308174,  2.24418,  1.28387, 3.80846,1.]

    maxrna = 25
    detectionLimitLC = 1
    detectionLimitNS = 2

    timevec = 1:1:200
    timevec_on = 1:1:10
    timevec_intensity = 1:1:20
    
    (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =  ModelOutput(model1, parameters, maxrna,detectionLimitLC, detectionLimitNS, timevec_on[end],timevec[end],timevec[end],timevec_intensity[end])

    (P,ssp, stateTr, stateTr_on, stateAbs_on, weightsTr_off,PabsOff, sspTr_off,Pabs) = StoThyLiveCell.mo_basics(model1, parameters, maxrna, detectionLimitLC, detectionLimitNS,) 

    mnascent_s = StoThyLiveCell.mo_mnascent(ssp, maxrna, stateTr, model1.nbstate, detectionLimitNS) 
    survivalon_s = StoThyLiveCell.mo_ontime(P, ssp,stateTr_on, stateAbs_on,timevec_on)
    survivaloff_s = StoThyLiveCell.mo_offtime(PabsOff, weightsTr_off,timevec)
    survivalnb_s = StoThyLiveCell.mo_nextbursttime(sspTr_off,PabsOff,timevec)
    pburst_s = StoThyLiveCell.mo_pon(ssp,stateTr_on)
    corr_s = StoThyLiveCell.mo_interburstcorr(P, weightsTr_off,stateAbs_on, stateTr_on, 15000) 
    avgint_s =StoThyLiveCell.mo_avgintensity(detectionLimitLC, P,Pabs, sspTr_off,stateAbs_on, stateTr_on,timevec_intensity, model1.nbstate, maxrna)
    @test mean_nascentrna ≈ mnascent_s
    @test prob_burst ≈ pburst_s
    @test correlation_interburst ≈ corr_s
    @test intensity_burst[1] ≈ avgint_s[1]
    @test intensity_burst[20] ≈ avgint_s[20]
    @test survival_nextburst[1] ≈ survivalnb_s[1]
    @test survival_nextburst[200] ≈  survivalnb_s[200]
    @test survival_interburst[1] ≈  survivaloff_s[1]
    @test survival_interburst[200] ≈   survivaloff_s[200]
    @test survival_burst[1] ≈ survivalon_s[1]
    @test survival_burst[10] ≈  survivalon_s[10]
end

@testset "Test single for only on/off" begin
    Qstate = [0    8    4    0    0    0;
                7    0    0    4    0    0;
                3    0    0    8    2    0;
                0    3    6    0    0    2;
                0    0    1    0    0    8;
                0    0    0    1    5    0]
    paramToRate_idx = findall(Qstate .>0)
    paramToRate_val = Qstate[findall(Qstate .>0)]
    model1 = StoThyLiveCell.StandardStoModel(6,8,1,paramToRate_idx,paramToRate_val,[1,3,5],[9,9,9],10)

    #creating an instance 
    parameters = [0.0178504,  0.0436684,  0.0543096,  0.427785,  0.023986,  0.308174,  2.24418,  1.28387, 3.80846,1.]

    maxrna = 25
    detectionLimitLC = 1
    detectionLimitNS = 2

    timevec = 1:1:200
    timevec_on = 1:1:10
    timevec_intensity = 1:1:20
    
 
    (P,ssp, stateTr, stateTr_on, stateAbs_on, weightsTr_off,PabsOff, sspTr_off,Pabs) = StoThyLiveCell.mo_basics(model1, parameters, maxrna, detectionLimitLC, detectionLimitNS) 
    mnascent_s = StoThyLiveCell.mo_mnascent(ssp, maxrna, stateTr, model1.nbstate, detectionLimitNS) 
    survivalon_s = StoThyLiveCell.mo_ontime(P, ssp,stateTr_on, stateAbs_on,timevec_on)
    survivaloff_s = StoThyLiveCell.mo_offtime(PabsOff, weightsTr_off,timevec)
    pburst_s = StoThyLiveCell.mo_pon(ssp,stateTr_on)
    corr_s = StoThyLiveCell.mo_interburstcorr(P, weightsTr_off,stateAbs_on, stateTr_on, 15000) 

   StoThyLiveCell.mo_basics!(model1, parameters, maxrna, P,ssp, stateTr_on, stateAbs_on, weightsTr_off,PabsOff) 
    mnascent_s2 = StoThyLiveCell.mo_mnascent(ssp, maxrna, stateTr, model1.nbstate, detectionLimitNS) 
    survivalon_s2 = StoThyLiveCell.mo_ontime(P, ssp,stateTr_on, stateAbs_on,timevec_on)
    survivaloff_s2 = StoThyLiveCell.mo_offtime(PabsOff, weightsTr_off,timevec)
    pburst_s2 = StoThyLiveCell.mo_pon(ssp,stateTr_on)
    corr_s2 = StoThyLiveCell.mo_interburstcorr(P, weightsTr_off,stateAbs_on, stateTr_on, 15000) 

    @test mnascent_s ≈ mnascent_s2
    @test pburst_s≈ pburst_s2
    @test  corr_s ≈ corr_s2
    @test survivaloff_s[1] ≈  survivaloff_s2[1]
    @test survivaloff_s[200] ≈   survivaloff_s2[200]
    @test survivalon_s[1] ≈ survivalon_s2[1]
    @test survivalon_s[10] ≈  survivalon_s2[10]
end


@testset "Test single for only nextburst" begin
    Qstate = [0    8    4    0    0    0;
                7    0    0    4    0    0;
                3    0    0    8    2    0;
                0    3    6    0    0    2;
                0    0    1    0    0    8;
                0    0    0    1    5    0]
    paramToRate_idx = findall(Qstate .>0)
    paramToRate_val = Qstate[findall(Qstate .>0)]
    model1 = StoThyLiveCell.StandardStoModel(6,8,1,paramToRate_idx,paramToRate_val,[1,3,5],[9,9,9],10)

    #creating an instance 
    parameters = [0.0178504,  0.0436684,  0.0543096,  0.427785,  0.023986,  0.308174,  2.24418,  1.28387, 3.80846,1.]

    maxrna = 25
    detectionLimitLC = 1
    detectionLimitNS = 2

    timevec = 1:1:200
    timevec_on = 1:1:10
    timevec_intensity = 1:1:20
    
 
    (P,ssp, stateTr, stateTr_on, stateAbs_on, weightsTr_off,PabsOff, sspTr_off,Pabs) = StoThyLiveCell.mo_basics(model1, parameters, maxrna, detectionLimitLC, detectionLimitNS) 
    survivalnb_s = StoThyLiveCell.mo_nextbursttime(sspTr_off,PabsOff,timevec)

    StoThyLiveCell.mo_basics!(model1, parameters, maxrna, P,ssp, stateTr_on, stateAbs_on, weightsTr_off,PabsOff,sspTr_off) 
    survivalnb_s2 = StoThyLiveCell.mo_nextbursttime(sspTr_off,PabsOff,timevec)

    @test survivalnb_s[1] ≈  survivalnb_s2[1]

end



@testset "Test optim error function" begin

    datafile= load("./data_test.jld2") ;
    data_test = datafile["data_all"];

    datatype = (StoThyLiveCell.Survival_InterBurst(),StoThyLiveCell.Survival_Burst(),StoThyLiveCell.Mean_Nascent(), StoThyLiveCell.Prob_Burst(), StoThyLiveCell.Correlation_InterBurst(),)
    datalist = data_test[[2,1,5,6,7]]
    datagroup = StoThyLiveCell.LiveCellData()
    dist = (StoThyLiveCell.LsqSurvival(), StoThyLiveCell.LsqSurvival(), StoThyLiveCell.LsqNumber(), StoThyLiveCell.LsqProb(), StoThyLiveCell.LsqNumber(),)
    maxrnaLC = 10
    maxrnaFC = 40
    detectionLimitLC = 1
    detectionLimitNS = 2
    burstsinglet = :with

    data = StoThyLiveCell.DataFit{typeof(datatype),typeof(datalist)}(datatype, datagroup, datalist,detectionLimitLC, detectionLimitNS, burstsinglet)

    #model
    Qstate = [0    8    4    0    0    0;
    7    0    0    4    0    0;
    3    0    0    8    2    0;
    0    3    6    0    0    2;
    0    0    1    0    0    8;
    0    0    0    1    5    0]
    paramToRate_idx = findall(Qstate .>0)
    paramToRate_val = Qstate[findall(Qstate .>0)]
    model = StoThyLiveCell.StandardStoModel(6,8,1,paramToRate_idx,paramToRate_val,[1,3,5],[9,9,9],10)

    #setting up the optimiziation
    optimtest = StoThyLiveCell.OptimStruct{typeof(data), typeof(dist), typeof(model)}(data,dist,model)

    SRange = [(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),]

    FRange = [(0,200),(0,7),(0,0),(0,0),(0,0),]
    fixedparameters = [1.]
    #indices of the free parameters
    freeparametersidx = [1,2,3,4,5,6,7,8,9]


    err_func = StoThyLiveCell.ini_optim(optimtest, optimtest.data.datagroup)
    
    data_fit = StoThyLiveCell.ini_data(optimtest, FRange)


    freeparameters = [.01,.01,.01,.01,.1,.1,.1,.1,10]

    (P,ssp, stateTr, stateTr_on, stateAbs_on, weightsTr_off,PabsOff, sspTr_Off, Pabs ) = StoThyLiveCell.mo_basics(model, zeros(model.nbparameters+model.nbkini+1), maxrnaLC, data.detectionLimitLC, data.detectionLimitNS) 
    Qrna = zeros(model.nbstate*(maxrnaFC+1),model.nbstate*(maxrnaFC+1))
    utileMat = (stateTr=stateTr, stateTr_on=stateTr_on, stateAbs_on=stateAbs_on, weightsTr_off=weightsTr_off, P=P, ssp=ssp, PabsOff=PabsOff, sspTr_Off=sspTr_Off, Pabs=Pabs, Qrna=Qrna)
    

    optim_struct_wrapper = StoThyLiveCell.OptimStructWrapper{typeof(optimtest.data),typeof(optimtest.dist), typeof(optimtest.model),typeof(err_func), typeof(utileMat)}(optimtest.data, data_fit, optimtest.dist, optimtest.model, SRange, maxrnaLC, maxrnaFC, freeparametersidx,fixedparameters, utileMat, err_func)

   @btime $err_func($freeparameters,$optim_struct_wrapper )

    @test err_func(freeparameters,optim_struct_wrapper ) >= 0
end


@testset "Test optim for live cells" begin

    datafile= load("./data_test.jld2") ;
    data_test = datafile["data_all"];

    datatype = (StoThyLiveCell.Survival_InterBurst(),StoThyLiveCell.Survival_Burst(),StoThyLiveCell.Mean_Nascent(), StoThyLiveCell.Prob_Burst(), StoThyLiveCell.Correlation_InterBurst(),)
    datalist = data_test[[2,1,5,6,7]]
    datagroup = StoThyLiveCell.LiveCellData()
    dist = (StoThyLiveCell.LsqSurvival(), StoThyLiveCell.LsqSurvival(), StoThyLiveCell.LsqNumber(), StoThyLiveCell.LsqProb(), StoThyLiveCell.LsqNumber(),)
    maxrnaLC = 10
    maxrnaFC = 40
    detectionLimitLC = 1
    detectionLimitNS = 2
    burstsinglet = :with

    data = StoThyLiveCell.DataFit{typeof(datatype),typeof(datalist)}(datatype, datagroup, datalist,detectionLimitLC, detectionLimitNS, burstsinglet)

    #model
    Qstate = [0    8    4    0    0    0;
    7    0    0    4    0    0;
    3    0    0    8    2    0;
    0    3    6    0    0    2;
    0    0    1    0    0    8;
    0    0    0    1    5    0]
    paramToRate_idx = findall(Qstate .>0)
    paramToRate_val = Qstate[findall(Qstate .>0)]
    model = StoThyLiveCell.StandardStoModel(6,8,1,paramToRate_idx,paramToRate_val,[1,3,5],[9,9,9],10)

    #setting up the optimiziation
    optimtest = StoThyLiveCell.OptimStruct{typeof(data), typeof(dist), typeof(model)}(data,dist,model)

    SRange = [(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),]

    FRange = [(0,200),(0,7),(0,0),(0,0),(0,0),]
    fixedparameters = [1.]
    #indices of the free parameters
    freeparametersidx = [1,2,3,4,5,6,7,8,9]

    sol, bfparameters, minval, minidx, estimate_signal = StoThyLiveCell.optim_function(SRange, FRange, optimtest; NbOptim=2, fixedparameters=fixedparameters,  freeparametersidx=freeparametersidx, maxrnaLC=maxrnaLC, maxrnaFC=maxrnaFC)

    println(sol[1].objective)
    @test typeof(sol[2].u) <: Vector
end


@testset "Test optim for fixed cells" begin

    datafile= load("./data_test.jld2") ;
    data_test = datafile["data_all"];

    datatype = (StoThyLiveCell.Distribution_RNA(),)
    datalist = data_test[[8]]
    datagroup = StoThyLiveCell.FixedCellData()
    dist = (StoThyLiveCell.LikelihoodRNA(),)
    maxrnaLC = 10
    maxrnaFC = maximum(datalist[1])
    detectionLimitLC = 1
    detectionLimitNS = 2
    burstsinglet = :with

    data = StoThyLiveCell.DataFit{typeof(datatype),typeof(datalist)}(datatype, datagroup, datalist,detectionLimitLC, detectionLimitNS, burstsinglet)

    #model
    Qstate = [0    8    4    0    0    0;
    7    0    0    4    0    0;
    3    0    0    8    2    0;
    0    3    6    0    0    2;
    0    0    1    0    0    8;
    0    0    0    1    5    0]
    paramToRate_idx = findall(Qstate .>0)
    paramToRate_val = Qstate[findall(Qstate .>0)]
    model = StoThyLiveCell.StandardStoModel(6,8,1,paramToRate_idx,paramToRate_val,[1,3,5],[9,9,9],10)

    #setting up the optimiziation
    optimtest = StoThyLiveCell.OptimStruct{typeof(data), typeof(dist), typeof(model)}(data,dist,model)

    SRange = [(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),]

    FRange = [(0,55),]
    fixedparameters = [1.]
    #indices of the free parameters
    freeparametersidx = [1,2,3,4,5,6,7,8,9]

    sol, bfparameters, minval, minidx, estimate_signal = StoThyLiveCell.optim_function(SRange, FRange, optimtest; NbOptim=2, fixedparameters=fixedparameters,  freeparametersidx=freeparametersidx, maxrnaLC=maxrnaLC, maxrnaFC=maxrnaFC)

    @test typeof(sol[2].u) <: Vector
end

@testset "Test optim for mixture data" begin

    datafile= load("./data_test.jld2") ;
    data_test = datafile["data_all"];

    datatype = (StoThyLiveCell.Survival_InterBurst(),StoThyLiveCell.Survival_Burst(), StoThyLiveCell.Distribution_RNA(),)
    datalist = data_test[[2,1,8]]
    datagroup = StoThyLiveCell.FixedAndLiveCellData()
    dist = (StoThyLiveCell.LsqSurvival(), StoThyLiveCell.LsqSurvival(), StoThyLiveCell.LikelihoodRNA(),)
    maxrnaLC = 10
    maxrnaFC = 55
    detectionLimitLC = 1
    detectionLimitNS = 2
    burstsinglet = :with

    data = StoThyLiveCell.DataFit{typeof(datatype),typeof(datalist)}(datatype, datagroup, datalist,detectionLimitLC, detectionLimitNS, burstsinglet)

    #model
    Qstate = [0    8    4    0    0    0;
    7    0    0    4    0    0;
    3    0    0    8    2    0;
    0    3    6    0    0    2;
    0    0    1    0    0    8;
    0    0    0    1    5    0]
    paramToRate_idx = findall(Qstate .>0)
    paramToRate_val = Qstate[findall(Qstate .>0)]
    model = StoThyLiveCell.StandardStoModel(6,8,1,paramToRate_idx,paramToRate_val,[1,3,5],[9,9,9],10)

    #setting up the optimiziation
    optimtest = StoThyLiveCell.OptimStruct{typeof(data), typeof(dist), typeof(model)}(data,dist,model)

    SRange = [(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0, .1),]

    FRange = [(0,200),(0,7), (0,55),]
    fixedparameters = [1.]
    #indices of the free parameters
    freeparametersidx = [1,2,3,4,5,6,7,8,9,11]

    sol, bfparameters, minval, minidx, estimate_signal = StoThyLiveCell.optim_function(SRange, FRange, optimtest; NbOptim=2, fixedparameters=fixedparameters,  freeparametersidx=freeparametersidx, maxrnaLC=maxrnaLC, maxrnaFC=maxrnaFC)

    println(bfparameters)
    @test typeof(sol[2].u) <: Vector
end
