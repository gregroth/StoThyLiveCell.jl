using StoThyLiveCell
using Test


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


#= @testset "StoThyLiveCell.jl" begin
    # Write your tests here.
end
 =#