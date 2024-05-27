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

# 0.1 0.9 
# Local size 2
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size96/localsize2/RGstep1(8,4,2,2)")
traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size96_loc2_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 3
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size96/localsize3/RGstep1(24,6,4,2)")
traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size96_loc3_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 4
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size96/localsize4/RGstep1(48,8,6,2)")
traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size96_loc4_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 6
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size96/localsize6/RGstep1(120,12,10,2)")
traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size96_loc6_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 8
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size96/localsize8/RGstep1(224,16,14,2)")
traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size96_loc8_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 12
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size96/localsize12/RGstep1(528,24,22,2)")
traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size96_loc12_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 16
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size96/localsize16/RGstep1(960,32,30,2)")
traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size96_loc16_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 24
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size96/localsize24/RGstep1(2208,48,46,2)")
traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size96_loc24_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

errorlist_trivial_t_a_neg_01_t_b_neg_09_size96 = log.([traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size96_loc2_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size96_loc3_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size96_loc4_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size96_loc6_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size96_loc8_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size96_loc12_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size96_loc16_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size96_loc24_RGstep1])
# radius96 = ([8/(2^2*2^2*2),24/(2^2*3^2*2),48/(2^2*4^2*2),120/(2^2*6^2*2),224/(2^2*8^2*2),528/(2^2*12^2*2),960/(2^2*16^2*2),2208/(2^2*24^2*2)])
radius96 = log.([2,3,4,6,8,12,16,24])
scatter(radius96,errorlist_trivial_t_a_neg_01_t_b_neg_09_size96 ,mode="markers")
A3 = hcat(radius,ones(length(radius)))
slope3, intercept3 = A3 \ errorlist_trivial_t_a_neg_01_t_b_neg_09_size96 

fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size96/localsize2/RGstep1")
traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size96_loc2_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

#(1/2)
# Local size 3
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size96/localsize3/RGstep1(18,9,5,4)")
traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size96_loc3_RGstep1_ratio = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 4
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size96/localsize4/RGstep1(32,16,8,8)")
traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size96_loc4_RGstep1_ratio = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 6
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size96/localsize6/RGstep1(72,36,18,18)")
traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size96_loc6_RGstep1_ratio = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 8
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size96/localsize8/RGstep1(128,64,32,32)")
traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size96_loc8_RGstep1_ratio = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 12
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size96/localsize12/RGstep1(288,144,72,72)")
traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size96_loc12_RGstep1_ratio = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

errorlist_trivial_t_a_neg_01_t_b_neg_09_size96 = log.([
                                                    traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size96_loc3_RGstep1_ratio,
                                                    traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size96_loc4_RGstep1_ratio,
                                                    traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size96_loc6_RGstep1_ratio,
                                                    traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size96_loc8_RGstep1_ratio,
                                                    traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size96_loc12_RGstep1_ratio
                                                    ])
radius = ([3,4,6,8,12])
scatter(radius,errorlist_trivial_t_a_neg_01_t_b_neg_09_size96 ,mode="markers")

# 0.2 0.8 
# Local size 2
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_02_t_b_neg_08/size96/localsize2/RGstep1")
traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc2_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 3
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_02_t_b_neg_08/size96/localsize3/RGstep1")
traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc3_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 4
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_02_t_b_neg_08/size96/localsize4/RGstep1")
traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc4_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 6
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_02_t_b_neg_08/size96/localsize6/RGstep1")
traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc6_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 8
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_02_t_b_neg_08/size96/localsize8/RGstep1")
traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc8_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 12
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_02_t_b_neg_08/size96/localsize12/RGstep1")
traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc12_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 16
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_02_t_b_neg_08/size96/localsize16/RGstep1")
traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc16_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 24
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_02_t_b_neg_08/size96/localsize24/RGstep1")
traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc24_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

errorlist_trivial_t_a_neg_02_t_b_neg_08_size96 = log.([traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc2_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc3_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc4_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc6_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc8_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc12_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc16_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc24_RGstep1])
radius = ([2,3,4,6,8,12,16,24])
scatter!(radius,errorlist_trivial_t_a_neg_02_t_b_neg_08_size96 ,mode="markers")


# 0.3 0.7 
# Local size 2
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_03_t_b_neg_07/size96/localsize2/RGstep1")
traceL1norm_trivial_t_a_neg_03_t_b_neg_07_size96_loc2_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 3
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_03_t_b_neg_07/size96/localsize3/RGstep1")
traceL1norm_trivial_t_a_neg_03_t_b_neg_07_size96_loc3_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 4
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_03_t_b_neg_07/size96/localsize4/RGstep1")
traceL1norm_trivial_t_a_neg_03_t_b_neg_07_size96_loc4_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 6
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_02_t_b_neg_08/size96/localsize6/RGstep1")
traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc6_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 8
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_02_t_b_neg_08/size96/localsize8/RGstep1")
traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc8_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 12
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_02_t_b_neg_08/size96/localsize12/RGstep1")
traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc12_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 16
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_02_t_b_neg_08/size96/localsize16/RGstep1")
traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc16_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 24
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_02_t_b_neg_08/size96/localsize24/RGstep1")
traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc24_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

errorlist_trivial_t_a_neg_02_t_b_neg_08_size96 = log.([traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc2_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc3_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc4_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc6_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc8_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc12_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc16_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc24_RGstep1])
radius = ([2,3,4,6,8,12,16,24])
scatter!(radius,errorlist_trivial_t_a_neg_02_t_b_neg_08_size96 ,mode="markers")


# 0.4 0.6 
# Local size 2
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_02_t_b_neg_08/size96/localsize2/RGstep1")
traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc2_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 3
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_02_t_b_neg_08/size96/localsize3/RGstep1")
traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc3_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 4
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_02_t_b_neg_08/size96/localsize4/RGstep1")
traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc4_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 6
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_02_t_b_neg_08/size96/localsize6/RGstep1")
traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc6_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 8
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_02_t_b_neg_08/size96/localsize8/RGstep1")
traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc8_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 12
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_02_t_b_neg_08/size96/localsize12/RGstep1")
traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc12_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 16
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_02_t_b_neg_08/size96/localsize16/RGstep1")
traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc16_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 24
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_02_t_b_neg_08/size96/localsize24/RGstep1")
traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc24_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

errorlist_trivial_t_a_neg_02_t_b_neg_08_size96 = log.([traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc2_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc3_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc4_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc6_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc8_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc12_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc16_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_02_t_b_neg_08_size96_loc24_RGstep1])
radius = ([2,3,4,6,8,12,16,24])
scatter!(radius,errorlist_trivial_t_a_neg_02_t_b_neg_08_size96 ,mode="markers")

