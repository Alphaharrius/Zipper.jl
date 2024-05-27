using Zipper
using Plots
plotlyjs()


function focktraceL1norm(fockmap,volume)
    @info("Calculating L1norm...")
    modeevalpairs = fockmap|>eigspech|>geteigenvalues
    return(sum([abs(modeevalpair[2]) for modeevalpair in modeevalpairs])/volume)
end

function startingcorrelations(correlations::FockMap, unitcellscale::Number)
    crystalfock = correlations|>getoutspace
    scale = Scale([unitcellscale 0; 0 unitcellscale], crystalfock|>getcrystal|>getspace)
    @info("Performing unitcellblocking...")
    @info("Generating unitcellblocking transformation...")
    blocker = @time scale * crystalfock
    @info("Performing unitcellblocking on correlations...")
    blockedcorrelations = @time blocker * correlations * blocker'
    return blockedcorrelations
end

# Trivial fixed RG level 

# -0.1, -0.9
# 0.1 0.9 
# Local size 2
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size48/localsize2/RGstep1(8,4,2,2)")
traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc2_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 3
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size48/localsize3/RGstep1(24,6,4,2)")
traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc3_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 4
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size48/localsize4/RGstep1(48,8,6,2)")
traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc4_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 6
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size48/localsize6/RGstep1(120,12,10,2)")
traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc6_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 8
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size48/localsize8/RGstep1(224,16,14,2)")
traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc8_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 12
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size48/localsize12/RGstep1(528,24,22,2)")
traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc12_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)


errorlist_trivial_t_a_neg_01_t_b_neg_09_size48 = log.([traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc2_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc3_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc4_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc6_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc8_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc12_RGstep1])

# radius48 = ([8/(2^2*2^2*2),24/(2^2*3^2*2),48/(2^2*4^2*2),120/(2^2*6^2*2),224/(2^2*8^2*2),528/(2^2*12^2*2)])
radius48 = log.([2,3,4,6,8,12])
scatter(radius48,errorlist_trivial_t_a_neg_01_t_b_neg_09_size48 ,mode="markers")


#L2 norm
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size48/localsize2/RGstep1(8,4,2,2)")
traceL2norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc2_RGstep1 = fioload("traceL2norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 3
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size48/localsize3/RGstep1(24,6,4,2)")
traceL2norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc3_RGstep1 = fioload("traceL2norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 4
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size48/localsize4/RGstep1(48,8,6,2)")
traceL2norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc4_RGstep1 = fioload("traceL2norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 6
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size48/localsize6/RGstep1(120,12,10,2)")
traceL2norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc6_RGstep1 = fioload("traceL2norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 8
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size48/localsize8/RGstep1(224,16,14,2)")
traceL2norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc8_RGstep1 = fioload("traceL2norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 12
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size48/localsize12/RGstep1(528,24,22,2)")
traceL2norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc12_RGstep1 = fioload("traceL2norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)


errorlist_trivial_t_a_neg_01_t_b_neg_09_size48_L2norm = log.([traceL2norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc2_RGstep1,
                                                    traceL2norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc3_RGstep1,
                                                    traceL2norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc4_RGstep1,
                                                    traceL2norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc6_RGstep1,
                                                    traceL2norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc8_RGstep1,
                                                    traceL2norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc12_RGstep1])

# radius48 = ([8/(2^2*2^2*2),24/(2^2*3^2*2),48/(2^2*4^2*2),120/(2^2*6^2*2),224/(2^2*8^2*2),528/(2^2*12^2*2)])
radius48 = log.([2,3,4,6,8,12])
scatter(radius48,errorlist_trivial_t_a_neg_01_t_b_neg_09_size48_L2norm ,mode="markers")

#across RG L1 norm

fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size48/localsize2/RGstep1(4,4,4,4)")
traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc2_RGstep1 = fioload("traceL1norm")

fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size48/localsize2/RGstep2(4,4,4,4)")
traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc2_RGstep2 = fioload("traceL1norm")

fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size48/localsize2/RGstep3(4,4,4,4)")
traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc2_RGstep3 = fioload("traceL1norm")

fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size48/localsize2/RGstep4(4,4,4,4)")
traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc2_RGstep4 = fioload("traceL1norm")

fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_02_t_b_neg_08/size48/localsize2/RGstep1(4,4,4,4)")
traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size48_loc2_RGstep1 = fioload("traceL1norm")

fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_02_t_b_neg_08/size48/localsize2/RGstep2(4,4,4,4)")
traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size48_loc2_RGstep2 = fioload("traceL1norm")

fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_02_t_b_neg_08/size48/localsize2/RGstep3(4,4,4,4)")
traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size48_loc2_RGstep3 = fioload("traceL1norm")

fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_02_t_b_neg_08/size48/localsize2/RGstep4(4,4,4,4)")
traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size48_loc2_RGstep4 = fioload("traceL1norm")

fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_03_t_b_neg_07/size48/localsize2/RGstep1(4,4,4,4)")
traceL1norm_trivial_t_a_neg_03_t_b_neg_07_size48_loc2_RGstep1 = fioload("traceL1norm")

fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_03_t_b_neg_07/size48/localsize2/RGstep2(4,4,4,4)")
traceL1norm_trivial_t_a_neg_03_t_b_neg_07_size48_loc2_RGstep2 = fioload("traceL1norm")

fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_03_t_b_neg_07/size48/localsize2/RGstep3(4,4,4,4)")
traceL1norm_trivial_t_a_neg_03_t_b_neg_07_size48_loc2_RGstep3 = fioload("traceL1norm")

fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_03_t_b_neg_07/size48/localsize2/RGstep4(4,4,4,4)")
traceL1norm_trivial_t_a_neg_03_t_b_neg_07_size48_loc2_RGstep4 = fioload("traceL1norm")

fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_04_t_b_neg_06/size48/localsize2/RGstep1(4,4,4,4)")
traceL1norm_trivial_t_a_neg_04_t_b_neg_06_size48_loc2_RGstep1 = fioload("traceL1norm")

fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_04_t_b_neg_06/size48/localsize2/RGstep2(4,4,4,4)")
traceL1norm_trivial_t_a_neg_04_t_b_neg_06_size48_loc2_RGstep2 = fioload("traceL1norm")

fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_04_t_b_neg_06/size48/localsize2/RGstep3(4,4,4,4)")
traceL1norm_trivial_t_a_neg_04_t_b_neg_06_size48_loc2_RGstep3 = fioload("traceL1norm")

fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_04_t_b_neg_06/size48/localsize2/RGstep4(4,4,4,4)")
traceL1norm_trivial_t_a_neg_04_t_b_neg_06_size48_loc2_RGstep4 = fioload("traceL1norm")

fiodir("/Users/slwongag/Desktop/ZERData/GMERA/diracsemimetal/size48/localsize2/RGstep1(4,4,4,4)")
traceL1norm_diracsemimetal_size48_loc2_RGstep1 = fioload("traceL1norm")

fiodir("/Users/slwongag/Desktop/ZERData/GMERA/diracsemimetal/size48/localsize2/RGstep2(4,4,4,4)")
traceL1norm_diracsemimetal_size48_loc2_RGstep2 = fioload("traceL1norm")

fiodir("/Users/slwongag/Desktop/ZERData/GMERA/diracsemimetal/size48/localsize2/RGstep3(4,4,4,4)")
traceL1norm_diracsemimetal_size48_loc2_RGstep3 = fioload("traceL1norm")

fiodir("/Users/slwongag/Desktop/ZERData/GMERA/diracsemimetal/size48/localsize2/RGstep4(4,4,4,4)")
traceL1norm_diracsemimetal_size48_loc2_RGstep4 = fioload("traceL1norm")

fiodir("/Users/slwongag/Desktop/ZERData/GMERA/cherninsulator/tnn_01im/size48/RGstep1(4,4,4,4)")
traceL1norm_cherninsulator_tnn_01im_size48_loc2_RGstep1 = fioload("traceL1norm")

fiodir("/Users/slwongag/Desktop/ZERData/GMERA/cherninsulator/tnn_01im/size48/RGstep2(4,4,4,4)")
traceL1norm_cherninsulator_tnn_01im_size48_loc2_RGstep2 = fioload("traceL1norm")

fiodir("/Users/slwongag/Desktop/ZERData/GMERA/cherninsulator/tnn_01im/size48/RGstep3(4,4,4,4)")
traceL1norm_cherninsulator_tnn_01im_size48_loc2_RGstep3 = fioload("traceL1norm")

fiodir("/Users/slwongag/Desktop/ZERData/GMERA/cherninsulator/tnn_01im/size48/RGstep4(4,4,4,4)")
traceL1norm_cherninsulator_tnn_01im_size48_loc2_RGstep4 = fioload("traceL1norm")

fiodir("/Users/slwongag/Desktop/ZERData/GMERA/cherninsulator/tnn_02im/size48/RGstep1(4,4,4,4)")
traceL1norm_cherninsulator_tnn_02im_size48_loc2_RGstep1 = fioload("traceL1norm")

fiodir("/Users/slwongag/Desktop/ZERData/GMERA/cherninsulator/tnn_02im/size48/RGstep2(4,4,4,4)")
traceL1norm_cherninsulator_tnn_02im_size48_loc2_RGstep2 = fioload("traceL1norm")

fiodir("/Users/slwongag/Desktop/ZERData/GMERA/cherninsulator/tnn_02im/size48/RGstep3(4,4,4,4)")
traceL1norm_cherninsulator_tnn_02im_size48_loc2_RGstep3 = fioload("traceL1norm")

fiodir("/Users/slwongag/Desktop/ZERData/GMERA/cherninsulator/tnn_02im/size48/RGstep4(4,4,4,4)")
traceL1norm_cherninsulator_tnn_02im_size48_loc2_RGstep4 = fioload("traceL1norm")

errorlist_trivial_t_a_neg_01_t_b_neg_09_size48_across_RG = log.([traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc2_RGstep1,
                                                                traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc2_RGstep2,
                                                                traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc2_RGstep3,
                                                                traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size48_loc2_RGstep4])

errorlist_trivial_t_a_neg_02_t_b_neg_08_size48_across_RG = log.([traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size48_loc2_RGstep1,
                                                                traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size48_loc2_RGstep2,
                                                                traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size48_loc2_RGstep3,
                                                                traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size48_loc2_RGstep4])

errorlist_trivial_t_a_neg_03_t_b_neg_07_size48_across_RG = log.([traceL1norm_trivial_t_a_neg_03_t_b_neg_07_size48_loc2_RGstep1,
                                                                traceL1norm_trivial_t_a_neg_03_t_b_neg_07_size48_loc2_RGstep2,
                                                                traceL1norm_trivial_t_a_neg_03_t_b_neg_07_size48_loc2_RGstep3,
                                                                traceL1norm_trivial_t_a_neg_03_t_b_neg_07_size48_loc2_RGstep4])

errorlist_trivial_t_a_neg_04_t_b_neg_06_size48_across_RG = log.([traceL1norm_trivial_t_a_neg_04_t_b_neg_06_size48_loc2_RGstep1,
                                                                traceL1norm_trivial_t_a_neg_04_t_b_neg_06_size48_loc2_RGstep2,
                                                                traceL1norm_trivial_t_a_neg_04_t_b_neg_06_size48_loc2_RGstep3,
                                                                traceL1norm_trivial_t_a_neg_04_t_b_neg_06_size48_loc2_RGstep4])

errorlist_diracsemimetal_size48_across_RG = log.([traceL1norm_diracsemimetal_size48_loc2_RGstep1,
                                                traceL1norm_diracsemimetal_size48_loc2_RGstep2,
                                                traceL1norm_diracsemimetal_size48_loc2_RGstep3,
                                                traceL1norm_diracsemimetal_size48_loc2_RGstep4])

errorlist_cherninsulator_tnn_01im_size48_acrossRG = log.([traceL1norm_cherninsulator_tnn_01im_size48_loc2_RGstep1,
                                                        traceL1norm_cherninsulator_tnn_01im_size48_loc2_RGstep2,
                                                        traceL1norm_cherninsulator_tnn_01im_size48_loc2_RGstep3,
                                                        traceL1norm_cherninsulator_tnn_01im_size48_loc2_RGstep4])

errorlist_cherninsulator_tnn_02im_size48_acrossRG = log.([traceL1norm_cherninsulator_tnn_02im_size48_loc2_RGstep1,
                                                        traceL1norm_cherninsulator_tnn_02im_size48_loc2_RGstep2,
                                                        traceL1norm_cherninsulator_tnn_02im_size48_loc2_RGstep3,
                                                        traceL1norm_cherninsulator_tnn_02im_size48_loc2_RGstep4])

RGstep = log.([1,2,3,4])

scatter!(RGstep,errorlist_trivial_t_a_neg_01_t_b_neg_09_size48_across_RG ,mode="markers")
scatter!(RGstep,errorlist_trivial_t_a_neg_02_t_b_neg_08_size48_across_RG ,mode="markers")
scatter!(RGstep,errorlist_trivial_t_a_neg_03_t_b_neg_07_size48_across_RG ,mode="markers")
scatter!(RGstep,errorlist_trivial_t_a_neg_04_t_b_neg_06_size48_across_RG ,mode="markers")
scatter(RGstep,errorlist_diracsemimetal_size48_across_RG ,mode="markers")
scatter!(RGstep,errorlist_cherninsulator_tnn_01im_size48_acrossRG ,mode="markers")
scatter!(RGstep,errorlist_cherninsulator_tnn_02im_size48_acrossRG ,mode="markers")