using LinearAlgebra
using Zipper, Plots
plotlyjs()

fiodir("/Users/slwongag/Zipper.jl/ZERData/GMERA/threshold_1_001/trivial/t_a_neg_02_t_b_neg_08")

fioload("noofrgsteploc2size96")

ref = FockMap(fioload("correlations_size12"))
ref |>eigspech|>geteigenvalues
fioload("correlations_size12")|>eigspech|>geteigenvalues
fioload("errorloc2size12L1norm")
errorloc2size12L1norm = fioload("errorloc2size12L1norm")
errorloc2size24L1norm = fioload("errorloc2size24L1norm")
errorloc2size48L1norm = fioload("errorloc2size48L1norm")
errorloc2size96L1norm = fioload("errorloc2size96L1norm")

errorloc3size12L1norm = fioload("errorloc3size12L1norm")
errorloc3size24L1norm = fioload("errorloc3size24L1norm")
errorloc3size48L1norm = fioload("errorloc3size48L1norm")
# errorloc3size96L1norm = fioload("errorloc3size96L1norm")


errorloc4size24L1norm = fioload("errorloc4size24L1norm")
errorloc4size48L1norm = fioload("errorloc4size48L1norm")
# errorloc4size96L1norm = fioload("errorloc4size96L1norm")

errorloc6size24L1norm = fioload("errorloc6size24L1norm")
errorloc6size48L1norm = fioload("errorloc6size48L1norm")
# errorloc4size96L1norm = fioload("errorloc4size96L1norm")

fioload("noofcoremodesloc6size48")
errorloc8size48L1norm = fioload("errorloc8size48L1norm")
# errorloc8size96L1norm = fioload("errorloc8size96L1norm")

errorloc12size48L1norm = fioload("errorloc12size48L1norm")
# errorloc4size96L1norm = fioload("errorloc4size96L1norm")

radiuslist12 = [2,3]
errorlist12 = log.([errorloc2size12L1norm,errorloc3size12L1norm])
radiuslist24 = [2,3,4,6]
errorlist24 = log.([errorloc2size24L1norm,errorloc3size24L1norm,errorloc4size24L1norm,errorloc6size24L1norm])
radiuslist48 = [2,3,4,6,8,12]
errorlist48 = log.([errorloc2size48L1norm,errorloc3size48L1norm,errorloc4size48L1norm,errorloc6size48L1norm,errorloc8size48L1norm,errorloc12size48L1norm])


scatter(radiuslist12, errorlist12, mode="markers", name="size12")
scatter!(radiuslist24, errorlist24, mode="markers", name="size24")
scatter!(radiuslist48, errorlist48, mode="markers", name="size48")

plot([scatter(x=radiuslist12, y=errorlist12, mode="markers", name="size12"),scatter(x=radiuslist24, y=errorlist24, mode="markers", name="size24"),scatter(x=radiuslist48, y=errorlist48, mode="markers", name="size48")],
Layout(
        title="log(Tr(sqrt(C_orig-C_approx))) vs radius for trivial insulator ",
        template="simple_white"
    ))
