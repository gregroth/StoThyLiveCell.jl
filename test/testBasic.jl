using StoThyLiveCell
using Test

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
    model1 = StoThyLiveCell.StandardStoModel(6,8,paramToRate_idx,paramToRate_val,[1,3,5])

    #creating an instance 
    parameters = [0.0178504,  0.0436684,  0.0543096,  0.427785,  0.023986,  0.308174,  2.24418,  1.28387]
    kini = 3.80846
    delta = 1.
    maxrna = 25

    P1 = StoModel(model1, parameters,kini,delta, maxrna)

    timevec = 1:1:200
    timevec_on = 1:1:10
    timevec_int = 1:1:20

    (mnascentmrna_model, pburst_model, survivalspot_model,survivaldark_model, survivalnextburst_model, corr_interburst_model, intensity_model) =  ModelOutput(model1, parameters,kini,delta, maxrna,timevec_on[end],timevec[end],timevec[end],timevec_int[end])


    @test mnascentmrna_model ≈ 3.0337754651143656
    @test pburst_model ≈ 0.16568033517579464
    @test corr_interburst_model ≈ 0.10866194618122918
    @test intensity_model[1] ≈ 1.9075373067561898
    @test intensity_model[20] ≈ 0.0021212008496347893
    @test survivalnextburst_model[1] ≈ 0.9409352782206553
    @test survivalnextburst_model[200] ≈  0.0020411810555407066
    @test survivaldark_model[1] ≈ 0.8346193347912076
    @test survivaldark_model[200] ≈  0.0010301814686152551
    @test survivalspot_model[1] ≈ 0.742997360615419
    @test survivalspot_model[10] ≈  0.024598091184660095
end

@testset "Test single output fcts" begin
    (P,ssp, stateTr, stateTr_on, stateAbs_on, weightsTr_off) = mo_basics(model1, parameters,kini,delta, maxrna) 

    mnascent_s = mo_mnascent(ssp, maxrna, stateTr, model1.nbstate) 
    survivalon_s = mo_ontime(P, ssp,stateTr_on, stateAbs_on,timevec_on)
    survivaloff_s = mo_offtime(P, ssp,stateAbs_on, weightsTr_off,timevec)
    survivalnb_s = mo_nextbursttime(ssp,stateAbs_on,timevec)
    pburst_s = mo_pon(ssp,stateTr_on)
    corr_s = mo_interburstcorr(P, weightsTr_off,stateAbs_on, stateTr_on, 15000) 
    avgint_s =mo_avgintensity(P, sspTr_off,stateAbs_on, stateTr_on,timevec_intensity, nbstate, maxrna)
    @test mnascentmrna_model ≈ mnascent_s
    @test pburst_model ≈ pburst_s
    @test corr_interburst_model ≈ corr_s
    @test intensity_model[1] ≈ avgint_s[1]
    @test intensity_model[20] ≈ avgint_s[20]
    @test survivalnextburst_model[1] ≈ survivalnb_s[1]
    @test survivalnextburst_model[200] ≈  survivalnb_s[200]
    @test survivaldark_model[1] ≈  survivaloff_s[1]
    @test survivaldark_model[200] ≈   survivaloff_s[200]
    @test survivalspot_model[1] ≈ survivalon_s[1]
    @test survivalspot_model[10] ≈  survivalon_s[10]
end
#= @testset "StoThyLiveCell.jl" begin
    # Write your tests here.
end
 =#
