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
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size24/localsize2/RGstep1(8,4,2,2)")
traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size24_loc2_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 3
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size24/localsize3/RGstep1(24,6,4,2)")
traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size24_loc3_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 4
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size24/localsize4/RGstep1(48,8,6,2)")
traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size24_loc4_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

# Local size 6
fiodir("/Users/slwongag/Desktop/ZERData/GMERA/trivial/t_a_neg_01_t_b_neg_09/size24/localsize6/RGstep1(120,12,10,2)")
traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size24_loc6_RGstep1 = fioload("traceL1norm")
rg1firstcorrelationsresults=fioload("rg1firstcorrelationsresults")
visualize(rg1firstcorrelationsresults|>getinspace|>getcrystal|>getunitcell)

errorlist_trivial_t_a_neg_01_t_b_neg_09_size24 = ([traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size24_loc2_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size24_loc3_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size24_loc4_RGstep1,
                                                    traceL1norm_trivial_t_a_neg_01_t_b_neg_09_size24_loc6_RGstep1])

radius24 = ([8/(2^2*2^2*2),24/(2^2*3^2*2),48/(2^2*4^2*2),120/(2^2*6^2*2)])
scatter(radius24,errorlist_trivial_t_a_neg_01_t_b_neg_09_size24,mode="markers")