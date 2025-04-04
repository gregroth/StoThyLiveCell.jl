using StoThyLiveCell
using Test
using DataFrames, FileIO, JLD2
using BenchmarkTools
using Optimization
using OptimizationOptimJL
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

    timevec = collect(1:1:200)
    timevec_on = collect(1:1:10)
    timevec_int = collect(1:1:20)

    (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =  StoThyLiveCell.ModelOutput(model1, parameters, maxrna,detectionLimitLC, detectionLimitNS, timevec_on[end],timevec[end],timevec[end],timevec_int[end])


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

    timevec = collect(1:1:200)
    timevec_on = collect(1:1:10)
    timevec_int = collect(1:1:20)
    
    (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =  StoThyLiveCell.ModelOutput(model1, parameters, maxrna,detectionLimitLC, detectionLimitNS, timevec_on[end],timevec[end],timevec[end],timevec_int[end])



    (P,ssp, stateTr, stateTr_on, stateAbs_on, weightsTr_off,PabsOff, sspTr_off,Pabs) = StoThyLiveCell.mo_basics(model1, parameters, maxrna, detectionLimitLC, detectionLimitNS,) 

    mnascent_s = StoThyLiveCell.mean_nascentrna(ssp, maxrna, stateTr, model1.nbstate, detectionLimitNS) 
    survivalon_s = StoThyLiveCell.survival_burst(P, ssp,stateTr_on, stateAbs_on,timevec_on)
    survivaloff_s = StoThyLiveCell.survival_interburst(PabsOff, weightsTr_off,timevec)
    survivalnb_s = StoThyLiveCell.survival_nextburst(sspTr_off,PabsOff,timevec)
    pburst_s = StoThyLiveCell.prob_burst(ssp,stateTr_on)
    corr_s = StoThyLiveCell.correlation_interburst(P, weightsTr_off,stateAbs_on, stateTr_on, 15000) 
    avgint_s =StoThyLiveCell.intensity_burst(detectionLimitLC, P,Pabs, sspTr_off,stateTr_on, stateAbs_on,timevec_int, model1.nbstate, maxrna)
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

    timevec = collect(1:1:200)
    timevec_on = collect(1:1:10)
    timevec_intensity = collect(1:1:20)
    
 
    (P,ssp, stateTr, stateTr_on, stateAbs_on, weightsTr_off,PabsOff, sspTr_off,Pabs) = StoThyLiveCell.mo_basics(model1, parameters, maxrna, detectionLimitLC, detectionLimitNS) 
    mnascent_s = StoThyLiveCell.mean_nascentrna(ssp, maxrna, stateTr, model1.nbstate, detectionLimitNS) 
    survivalon_s = StoThyLiveCell.survival_burst(P, ssp,stateTr_on, stateAbs_on,timevec_on)
    survivaloff_s = StoThyLiveCell.survival_interburst(PabsOff, weightsTr_off,timevec)
    pburst_s = StoThyLiveCell.prob_burst(ssp,stateTr_on)
    corr_s = StoThyLiveCell.correlation_interburst(P, weightsTr_off,stateAbs_on, stateTr_on, 15000) 

   StoThyLiveCell.mo_basics!(model1, parameters, maxrna, P,ssp, stateTr_on, stateAbs_on, weightsTr_off,PabsOff) 
    mnascent_s2 = StoThyLiveCell.mean_nascentrna(ssp, maxrna, stateTr, model1.nbstate, detectionLimitNS) 
    survivalon_s2 = StoThyLiveCell.survival_burst(P, ssp,stateTr_on, stateAbs_on,timevec_on)
    survivaloff_s2 = StoThyLiveCell.survival_interburst(PabsOff, weightsTr_off,timevec)
    pburst_s2 = StoThyLiveCell.prob_burst(ssp,stateTr_on)
    corr_s2 = StoThyLiveCell.correlation_interburst(P, weightsTr_off,stateAbs_on, stateTr_on, 15000) 

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

    timevec = collect(1:1:200)
    timevec_on = collect(1:1:10)
    timevec_intensity = collect(1:1:20)
    
 
    (P,ssp, stateTr, stateTr_on, stateAbs_on, weightsTr_off,PabsOff, sspTr_off,Pabs) = StoThyLiveCell.mo_basics(model1, parameters, maxrna, detectionLimitLC, detectionLimitNS) 
    survivalnb_s = StoThyLiveCell.survival_nextburst(sspTr_off,PabsOff,timevec)

    StoThyLiveCell.mo_basics!(model1, parameters, maxrna, P,ssp, stateTr_on, stateAbs_on, weightsTr_off,PabsOff,sspTr_off) 
    survivalnb_s2 = StoThyLiveCell.survival_nextburst(sspTr_off,PabsOff,timevec)

    @test survivalnb_s[1] ≈  survivalnb_s2[1]

end

@testset "Test optim error function" begin

    datafile= load("./data_test.jld2") ;
    data_test = datafile["data_all"];

    datatype = (StoThyLiveCell.Survival_InterBurst(),StoThyLiveCell.Survival_Burst(),StoThyLiveCell.Mean_Nascent(), StoThyLiveCell.Prob_Burst(), StoThyLiveCell.Correlation_InterBurst(),)
    datalist = data_test[[2,1,5,6,7]]
    datagroup = StoThyLiveCell.LiveCellData()
    dist = (StoThyLiveCell.LsqSurvival(1.0), StoThyLiveCell.LsqSurvival(1.0), StoThyLiveCell.LsqNumber(1.0), StoThyLiveCell.LsqProb(1.0), StoThyLiveCell.LsqNumber(1.0),)
    maxrnaLC = 10
    maxrnaFC = 40
    detectionLimitLC = 1
    detectionLimitNS = 2
    burstsinglet = :without

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
    optimtest = StoThyLiveCell.OptimStruct{typeof(data), typeof(dist), typeof(model)}(data,dist,model,true)

    SRange = [(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),]

    FRange = [(0,200),(0,7),(0,0),(0,0),(0,0),]
    fixedparameters = [1.]
    #indices of the free parameters
    freeparametersidx = [1,2,3,4,5,6,7,8,9]

if burstsinglet==:without
    err_func = StoThyLiveCell.ini_optim_withoutsinglet(optimtest, optimtest.data.datagroup)
else
    err_func = StoThyLiveCell.ini_optim(optimtest, optimtest.data.datagroup)
end
    data_fit = StoThyLiveCell.ini_data(optimtest, FRange)


    freeparameters = [.01,.01,.01,.01,.1,.1,.1,.1,10]
    freeparameters = [0.01310826520418289, 0.6672789074705456, 0.2753473176958067, 0.5606025216433805, 0.02547352716458086, 0.37848638394059114, 3.056353368641461, 15.786212952072278, 32.628147442566735]
   
 if  burstsinglet==:without
    (nascentbin, P, ssp, stateTr, stateTr_on, stateAbs_on, totnbs, Pwos, stateAbs_on_wos, statePre_on_wos, weightsAbs_off_wos, sspTr_off_wos, weightsAbs_on, sspPreB, weightsTr_on, PabsOff, weightsTr_on_wos, weightsAbsorbed_off_wos, sspwos, weightsPre_on_and_on, Rn, NR, Nc, Qn, Nn, weightsTr_off_wos, Pabs_wos, weightsON_wos, rnanbvec_on, weightsPre_on_wos) = StoThyLiveCell.mo_basics_wosinglet(model, ones(model.nbparameters+model.nbkini+1), maxrnaLC, data.detectionLimitLC, data.detectionLimitNS) 
    Qrna = zeros(model.nbstate*(maxrnaFC+1),model.nbstate*(maxrnaFC+1))
    utileMat = (nascentbin=nascentbin, P=P, ssp=ssp, stateTr=stateTr, stateTr_on=stateTr_on, stateAbs_on=stateAbs_on, totnbs=totnbs, Pwos=Pwos, stateAbs_on_wos=stateAbs_on_wos, statePre_on_wos=statePre_on_wos, weightsAbs_off_wos=weightsAbs_off_wos, sspTr_off_wos=sspTr_off_wos, weightsAbs_on=weightsAbs_on, sspPreB=sspPreB, weightsTr_on=weightsTr_on, PabsOff=PabsOff, weightsTr_on_wos=weightsTr_on_wos, weightsAbsorbed_off_wos=weightsAbsorbed_off_wos, sspwos=sspwos, weightsPre_on_and_on=weightsPre_on_and_on, Rn=Rn, NR=NR, Nc=Nc, Qn=Qn, Nn=Nn, weightsTr_off_wos=weightsTr_off_wos, Pabs_wos=Pabs_wos, weightsON_wos=weightsON_wos, rnanbvec_on=rnanbvec_on, weightsPre_on_wos=weightsPre_on_wos, Qrna=Qrna)
 else
    (P,ssp, stateTr, stateTr_on, stateAbs_on, weightsTr_off,PabsOff, sspTr_Off, Pabs ) = StoThyLiveCell.mo_basics(model, ones(model.nbparameters+model.nbkini+1), maxrnaLC, data.detectionLimitLC, data.detectionLimitNS) 
    Qrna = zeros(model.nbstate*(maxrnaFC+1),model.nbstate*(maxrnaFC+1))
    utileMat = (stateTr=stateTr, stateTr_on=stateTr_on, stateAbs_on=stateAbs_on, weightsTr_off=weightsTr_off, P=P, ssp=ssp, PabsOff=PabsOff, sspTr_Off=sspTr_Off, Pabs=Pabs, Qrna=Qrna)
 end 

    optim_struct_wrapper = StoThyLiveCell.OptimStructWrapper{typeof(optimtest.data),typeof(optimtest.dist), typeof(optimtest.model),typeof(err_func), typeof(utileMat)}(optimtest.data, data_fit, optimtest.dist, optimtest.model, SRange, maxrnaLC, maxrnaFC, freeparametersidx,fixedparameters, utileMat, err_func)

   @btime $err_func($freeparameters,$optim_struct_wrapper )
println(err_func(freeparameters,optim_struct_wrapper ))
    @test err_func(freeparameters,optim_struct_wrapper ) >= 0
end

@testset "Test optim for live cells" begin

    datafile= load("./data_test.jld2") ;
    data_test = datafile["data_all"];

    datatype = (StoThyLiveCell.Survival_InterBurst(),StoThyLiveCell.Survival_Burst(),StoThyLiveCell.Mean_Nascent(), StoThyLiveCell.Prob_Burst(), StoThyLiveCell.Correlation_InterBurst(),)
    datalist = data_test[[2,1,5,6,7]]
    datagroup = StoThyLiveCell.LiveCellData()
    dist = (StoThyLiveCell.LsqSurvival(1.0), StoThyLiveCell.LsqSurvival(1.0), StoThyLiveCell.LsqNumber(1.0), StoThyLiveCell.LsqProb(1.0), StoThyLiveCell.LsqNumber(1.0),)
    maxrnaLC = 10
    maxrnaFC = 40
    detectionLimitLC = 1
    detectionLimitNS = 2
    burstsinglet = :without

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
    optimtest = StoThyLiveCell.OptimStruct{typeof(data), typeof(dist), typeof(model)}(data,dist,model,  true)

    SRange = [(0.0,1.0),(0.0,1.0),(0.0,1.0),(0.0,1.0),(0.0,.1),(0.0,1.0),(0.0,5.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),]

    FRange = [(0,200),(0,7),(0,0),(0,0),(0,0),]
    fixedparameters = [1.]
    #indices of the free parameters
    freeparametersidx = [1,2,3,4,5,6,7,8,9]

    sol, bfparameters, minval, minidx, estimate_signal = StoThyLiveCell.optim_function(SRange, FRange, optimtest; NbOptim=2, maxtime=10, maxiters=10, fixedparameters=fixedparameters,  freeparametersidx=freeparametersidx, maxrnaLC=maxrnaLC, maxrnaFC=maxrnaFC,Method=BFGS())

    println(sol[1].objective)
    println(sol[1].u)
    @test typeof(sol[1].u) <: Vector
end 

 
@testset "Test optim for fixed cells" begin

    datafile= load("./data_test.jld2") ;
    data_test = datafile["data_all"];

    datatype = (StoThyLiveCell.Distribution_RNA(),)
    datalist = data_test[[8]]
    datagroup = StoThyLiveCell.FixedCellData()
    dist = (StoThyLiveCell.LikelihoodRNA(1.0),)
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
    optimtest = StoThyLiveCell.OptimStruct{typeof(data), typeof(dist), typeof(model)}(data,dist,model, true)

    SRange = [(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),]

    FRange = [(0,55),]
    fixedparameters = [1.]
    #indices of the free parameters
    freeparametersidx = [1,2,3,4,5,6,7,8,9]

    sol, bfparameters, minval, minidx, estimate_signal = StoThyLiveCell.optim_function(SRange, FRange, optimtest; NbOptim=2, fixedparameters=fixedparameters,  freeparametersidx=freeparametersidx, maxrnaLC=maxrnaLC, maxrnaFC=maxrnaFC, Method=BFGS())

    @test typeof(sol[2].u) <: Vector
end 

@testset "Test optim for mixture data" begin

    datafile= load("./data_test.jld2") ;
    data_test = datafile["data_all"];

    datatype = (StoThyLiveCell.Survival_InterBurst(),StoThyLiveCell.Survival_Burst(), StoThyLiveCell.Distribution_RNA(),)
    datalist = data_test[[2,1,8]]
    datagroup = StoThyLiveCell.FixedAndLiveCellData()
    dist = (StoThyLiveCell.LsqSurvival(1.0), StoThyLiveCell.LsqSurvival(1.0), StoThyLiveCell.LikelihoodRNA(1.0),)
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
    optimtest = StoThyLiveCell.OptimStruct{typeof(data), typeof(dist), typeof(model)}(data,dist,model,true)

    SRange = [(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0, .1),]

    FRange = [(0,200),(0,7), (0,55),]
    fixedparameters = [1.]
    #indices of the free parameters
    freeparametersidx = [1,2,3,4,5,6,7,8,9,11]

    sol, bfparameters, minval, minidx, estimate_signal = StoThyLiveCell.optim_function(SRange, FRange, optimtest; NbOptim=2, fixedparameters=fixedparameters,  freeparametersidx=freeparametersidx, maxrnaLC=maxrnaLC, maxrnaFC=maxrnaFC)

    println(bfparameters)
    @test typeof(sol[2].u) <: Vector
end




@testset "Test optim for live cells without singlets" begin

    datafile= load("./data_test.jld2") ;
    data_test = datafile["data_all"];

    datatype = (StoThyLiveCell.Survival_InterBurst(),StoThyLiveCell.Survival_Burst(),StoThyLiveCell.Mean_Nascent(), StoThyLiveCell.Prob_Burst(), StoThyLiveCell.Correlation_InterBurst(),)
    datalist = data_test[[2,1,5,6,7]]
    datagroup = StoThyLiveCell.LiveCellData()
    dist = (StoThyLiveCell.LsqSurvival(1.0), StoThyLiveCell.LsqSurvival(1.0), StoThyLiveCell.LsqNumber(1.0), StoThyLiveCell.LsqProb(1.0), StoThyLiveCell.LsqNumber(1.0),)
    maxrnaLC = 10
    maxrnaFC = 40
    detectionLimitLC = 1
    detectionLimitNS = 2
    burstsinglet = :without

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
    optimtest = StoThyLiveCell.OptimStruct{typeof(data), typeof(dist), typeof(model)}(data,dist,model,false)

    SRange = [(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),]

    FRange = [(0,200),(0,7),(0,0),(0,0),(0,0),]
    fixedparameters = [1.]
    #indices of the free parameters
    freeparametersidx = [1,2,3,4,5,6,7,8,9]

    sol, bfparameters, minval, minidx, estimate_signal = StoThyLiveCell.optim_function(SRange, FRange, optimtest; NbOptim=2, maxtime=1, maxiters=1,fixedparameters=fixedparameters,  freeparametersidx=freeparametersidx, maxrnaLC=maxrnaLC, maxrnaFC=maxrnaFC,Method=SAMIN())

    println(sol[1].objective)
    @test typeof(sol[2].u) <: Vector
end

@testset "Test optim for mixture data without singlets" begin

    datafile= load("./data_test.jld2") ;
    data_test = datafile["data_all"];

    datatype = (StoThyLiveCell.Survival_InterBurst(),StoThyLiveCell.Survival_Burst(), StoThyLiveCell.Distribution_RNA(),)
    datalist = data_test[[2,1,8]]
    datagroup = StoThyLiveCell.FixedAndLiveCellData()
    dist = (StoThyLiveCell.LsqSurvival(1.0), StoThyLiveCell.LsqSurvival(1.0), StoThyLiveCell.LikelihoodRNA(1.0),)
    maxrnaLC = 10
    maxrnaFC = 55
    detectionLimitLC = 1
    detectionLimitNS = 2
    burstsinglet = :without

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
    optimtest = StoThyLiveCell.OptimStruct{typeof(data), typeof(dist), typeof(model)}(data,dist,model,true)

    SRange = [(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0, .1),]

    FRange = [(0,200),(0,7), (0,55),]
    fixedparameters = [1.]
    #indices of the free parameters
    freeparametersidx = [1,2,3,4,5,6,7,8,9,11]

    sol, bfparameters, minval, minidx, estimate_signal = StoThyLiveCell.optim_function(SRange, FRange, optimtest; NbOptim=2, fixedparameters=fixedparameters,  freeparametersidx=freeparametersidx, maxrnaLC=maxrnaLC, maxrnaFC=maxrnaFC)

    println(sol[1].objective)
    @test typeof(sol[2].u) <: Vector
end

@testset "Test 2s3r model without singlet" begin
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
    parameters = [0.0045506,  0.00421043,  0.001,  0.0208397,  0.000605331,  0.0315689,  0.42817,  0.637767,  9.35007,  2.4761]
    maxrna = 15
    detectionLimitLC = 1
    detectionLimitNS = 2

    P1 = StoModel(model1, parameters, maxrna)

    timevec = 1:1:200
    timevec_on = 1:1:10
    timevec_int = 1:1:20

    (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =  StoThyLiveCell.ModelOutput_wosinglet(model1, parameters, maxrna,detectionLimitLC, detectionLimitNS, timevec_on[end],timevec[end],timevec[end],timevec_int[end])


    @test mean_nascentrna ≈ 3.6401850701364196
    @test prob_burst ≈ 0.039305138190572735
    @test correlation_interburst ≈ 0.07048666712782842
    @test intensity_burst[1] ≈ 0.9261364173944553
    @test intensity_burst[5] ≈ 0.20014864762599732
    @test survival_nextburst[2] ≈ 0.9772514543254125
    @test survival_nextburst[200] ≈  0.3281989569400399
    @test survival_interburst[2] ≈ 0.8849072429323206
    @test survival_interburst[200] ≈  0.11172174493227556
    @test survival_burst[2] ≈ 0.601658802857394
    @test survival_burst[10] ≈  0.01088478468427921
end

@testset "Test 2s3r model without singlet, basic matrices" begin
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
    parameters = [0.0045506,  0.00421043,  0.001,  0.0208397,  0.000605331,  0.0315689,  0.42817,  0.637767,  9.35007,  2.4761]
    maxrna = 15
    detectionLimitLC = 1
    detectionLimitNS = 2

    P1 = StoModel(model1, parameters, maxrna)

    timevec = collect(1:1:200)
    timevec_on = collect(1:1:10)
    timevec_int = collect(1:1:20)

    (survival_burst, survival_interburst, survival_nextburst, prob_burst, mean_nascentrna, correlation_interburst, intensity_burst) =  StoThyLiveCell.ModelOutput_wosinglet(model1, parameters, maxrna,detectionLimitLC, detectionLimitNS, timevec_on[end],timevec[end],timevec[end],timevec_int[end])


    (nascentbin, P, ssp, stateTr, stateTr_on, stateAbs_on, totnbs, Pwos, stateAbs_on_wos, statePre_on_wos, weightsAbs_off_wos, sspTr_off_wos, weightsAbs_on, sspPreB, weightsTr_on, PabsOff, weightsTr_on_wos, weightsAbsorbed_off_wos, sspwos, weightsPre_on_and_on, Rn, NR, Nc, Qn, Nn, weightsTr_off_wos, Pabs_wos, weightsON_wos, rnanbvec_on, weightsPre_on_wos) = StoThyLiveCell.mo_basics_wosinglet(model1, parameters, maxrna, detectionLimitLC, detectionLimitNS) 

    StoThyLiveCell.mo_basics_wosinglet!(model1, parameters, maxrna, P,ssp,stateTr_on, stateAbs_on, totnbs, Pwos, stateAbs_on_wos, statePre_on_wos, weightsAbs_off_wos, sspTr_off_wos, weightsAbs_on, sspPreB, weightsTr_on, PabsOff, weightsTr_on_wos, weightsAbsorbed_off_wos, sspwos, weightsPre_on_and_on, Rn, NR, Nc, Qn, Nn, weightsTr_off_wos, Pabs_wos, weightsON_wos) 

    mean_nascentrna2 = StoThyLiveCell.mean_nascentrna_wosinglet(ssp, nascentbin, stateTr, maxrna, detectionLimitNS) 
    survival_burst2 = StoThyLiveCell.survival_burst_wosinglet( P, stateTr_on, weightsTr_on,timevec_on) 
    survival_interburst2 = StoThyLiveCell.survival_interburst_wosinglet(PabsOff, weightsTr_on_wos,timevec)  
    survival_nextburst2 = StoThyLiveCell.survival_nextburst_wosinglet(weightsAbsorbed_off_wos,PabsOff, timevec)  
    prob_burst2 = StoThyLiveCell.prob_burst_wosinglet(sspwos,weightsPre_on_and_on, stateTr_on) 
    correlation_interburst2 = StoThyLiveCell.correlation_interburst_wosinglet(Rn, NR, Nc, Qn, Nn, weightsTr_off_wos, 15000) 
    intensity_burst2 = StoThyLiveCell.intensity_burst_wosinglet(rnanbvec_on, Pwos, Pabs_wos, weightsPre_on_wos, statePre_on_wos, stateTr_on, weightsON_wos,timevec_int)  


    @test mean_nascentrna ≈ mean_nascentrna2
    @test prob_burst ≈ prob_burst2
    @test correlation_interburst ≈  correlation_interburst2
    @test intensity_burst[1] ≈ intensity_burst2[1] 
    @test intensity_burst[5] ≈ intensity_burst2[5]
    @test survival_nextburst[2] ≈ survival_nextburst2[2]
    @test survival_nextburst[200] ≈  survival_nextburst2[200]
    @test survival_interburst[2] ≈ survival_interburst2[2]
    @test survival_interburst[200] ≈  survival_interburst2[200] 
    @test survival_burst[2] ≈ survival_burst2[2]
    @test survival_burst[10] ≈  survival_burst2[10] 

    #@test correlation_interburst  ≈  correlation_interburst_test 

end




@testset "Test optim error function, without singlet" begin

    datafile= load("./data_test.jld2") ;
    data_test = datafile["data_all"];

    datatype = (StoThyLiveCell.Survival_InterBurst(),StoThyLiveCell.Survival_Burst(),StoThyLiveCell.Mean_Nascent(), StoThyLiveCell.Prob_Burst(), StoThyLiveCell.Correlation_InterBurst(),)
    datalist = data_test[[2,1,5,6,7]]
    datagroup = StoThyLiveCell.LiveCellData()
    dist = (StoThyLiveCell.LsqSurvival(1.0), StoThyLiveCell.LsqSurvival(1.0), StoThyLiveCell.LsqNumber(1.0), StoThyLiveCell.LsqProb(1.0), StoThyLiveCell.LsqNumber(1.0),)
    maxrnaLC = 10
    maxrnaFC = 40
    detectionLimitLC = 1
    detectionLimitNS = 2
    burstsinglet = :without

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
    optim_struct = StoThyLiveCell.OptimStruct{typeof(data), typeof(dist), typeof(model)}(data,dist,model,true)

    SRange = [(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),]

    FRange = [(0,200),(0,7),(0,0),(0,0),(0,0),]
    fixedparameters = [1.]
    #indices of the free parameters
    freeparametersidx = [1,2,3,4,5,6,7,8,9]


    err_func = StoThyLiveCell.ini_optim_withoutsinglet(optim_struct, optim_struct.data.datagroup)
    
    data_fit = StoThyLiveCell.ini_data(optim_struct, FRange)


    freeparameters = [.01,.01,.01,.01,.1,.1,.1,.1,10]

    (nascentbin, P, ssp, stateTr, stateTr_on, stateAbs_on, totnbs, Pwos, stateAbs_on_wos, statePre_on_wos, weightsAbs_off_wos, sspTr_off_wos, weightsAbs_on, sspPreB, weightsTr_on, PabsOff, weightsTr_on_wos, weightsAbsorbed_off_wos, sspwos, weightsPre_on_and_on, Rn, NR, Nc, Qn, Nn, weightsTr_off_wos, Pabs_wos, weightsON_wos, rnanbvec_on, weightsPre_on_wos) = StoThyLiveCell.mo_basics_wosinglet(model, ones(model.nbparameters+model.nbkini+1), maxrnaLC, data.detectionLimitLC, data.detectionLimitNS) 
    Qrna = zeros(model.nbstate*(maxrnaFC+1),model.nbstate*(maxrnaFC+1))
    utileMat = (nascentbin=nascentbin, P=P, ssp=ssp, stateTr=stateTr, stateTr_on=stateTr_on, stateAbs_on=stateAbs_on, totnbs=totnbs, Pwos=Pwos, stateAbs_on_wos=stateAbs_on_wos, statePre_on_wos=statePre_on_wos, weightsAbs_off_wos=weightsAbs_off_wos, sspTr_off_wos=sspTr_off_wos, weightsAbs_on=weightsAbs_on, sspPreB=sspPreB, weightsTr_on=weightsTr_on, PabsOff=PabsOff, weightsTr_on_wos=weightsTr_on_wos, weightsAbsorbed_off_wos=weightsAbsorbed_off_wos, sspwos=sspwos, weightsPre_on_and_on=weightsPre_on_and_on, Rn=Rn, NR=NR, Nc=Nc, Qn=Qn, Nn=Nn, weightsTr_off_wos=weightsTr_off_wos, Pabs_wos=Pabs_wos, weightsON_wos=weightsON_wos, rnanbvec_on=rnanbvec_on, weightsPre_on_wos=weightsPre_on_wos, Qrna=Qrna)
    optim_struct_wrapper = StoThyLiveCell.OptimStructWrapper{typeof(optim_struct.data),typeof(optim_struct.dist), typeof(optim_struct.model),typeof(err_func), typeof(utileMat)}(optim_struct.data, data_fit, optim_struct.dist, optim_struct.model, SRange, maxrnaLC, maxrnaFC, freeparametersidx,fixedparameters, utileMat, err_func)

   @btime $err_func($freeparameters,$optim_struct_wrapper )

    @test err_func(freeparameters,optim_struct_wrapper ) >= 0
end 


@testset "Test optim for live cells without singlets TWO MODELS" begin

    datafile= load("./data_test.jld2") ;
    data_test = datafile["data_all"];

    datatype_m1 = (StoThyLiveCell.Survival_InterBurst(),StoThyLiveCell.Survival_Burst(),StoThyLiveCell.Mean_Nascent(), StoThyLiveCell.Prob_Burst(), StoThyLiveCell.Correlation_InterBurst(),)
    datalist_m1 = data_test[[2,1,5,6,7]]
    datagroup_m1 = StoThyLiveCell.LiveCellData()
    dist_m1 = (StoThyLiveCell.LsqSurvival(1.0), StoThyLiveCell.LsqSurvival(1.0), StoThyLiveCell.LsqNumber(1.0), StoThyLiveCell.LsqProb(1.0), StoThyLiveCell.LsqNumber(1.0),)
    
    datatype_m2 = (StoThyLiveCell.Survival_Burst(),)
    datalist_m2 = data_test[[1]]
    datagroup_m2 = StoThyLiveCell.LiveCellData()
    dist_m2 = (StoThyLiveCell.LsqSurvival(1.0),)
    

    maxrnaLC = 10
    maxrnaFC = 40
    detectionLimitLC = 1
    detectionLimitNS = 2
    burstsinglet = :without

    data_m1 = StoThyLiveCell.DataFit{typeof(datatype_m1),typeof(datalist_m1)}(datatype_m1, datagroup_m1, datalist_m1,detectionLimitLC, detectionLimitNS, burstsinglet)
    data_m2= StoThyLiveCell.DataFit{typeof(datatype_m2),typeof(datalist_m2)}(datatype_m2, datagroup_m2, datalist_m2,detectionLimitLC, detectionLimitNS, burstsinglet)

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
    optimtest_m1 = StoThyLiveCell.OptimStruct{typeof(data_m1), typeof(dist_m1), typeof(model)}(data_m1,dist_m1,model,false)
    optimtest_m2 = StoThyLiveCell.OptimStruct{typeof(data_m2), typeof(dist_m2), typeof(model)}(data_m2,dist_m2,model,false)

    SRange = [(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),(0.0,50.0),]

    FRange_m1 = [(0,200),(0,7),(0,0),(0,0),(0,0),]
    FRange_m2 = [(0,7),]
   
    
    freeparametersidx_m1 =[1,2,3,4,5,6,7,8,9,10]
    freeparametersidx_m2 =[1,2,3,4,5,6,7,8,9,10]
 
    paramidx_m1=[1,2,3,4,5,6,7,8,9,10]
    paramidx_m2=[11,2,3,4,5,6,7,8,9,10]
    
    pm1_infreeparam = [1,2,3,4,5,6,7,8,9,10]
    pm2_infreeparam = [11,2,3,4,5,6,7,8,9,10]

    fixedparameters_m1 =[-1]
    fixedparameters_m2 =[-1]


    sol, bfparameters, minval, minidx, estimate_signal = StoThyLiveCell.optim_function_multipleModels(SRange, FRange_m1,FRange_m2, optimtest_m1, optimtest_m2; NbOptim=2, maxtime=1, maxiters=1,freeparametersidx_m1 =freeparametersidx_m1, freeparametersidx_m2 =freeparametersidx_m2, fixedparameters_m1 =fixedparameters_m1, fixedparameters_m2 =fixedparameters_m2, paramidx_m1=paramidx_m1, paramidx_m2=paramidx_m2,  pm1_infreeparam= pm1_infreeparam,  pm2_infreeparam= pm2_infreeparam,maxrnaLC=maxrnaLC, maxrnaFC=maxrnaFC,Method=SAMIN())

    println(sol[1].objective)
    @test typeof(sol[2].u) <: Vector
end