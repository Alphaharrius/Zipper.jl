using Pkg
Pkg.resolve()
Pkg.Registry.update()
Pkg.instantiate()
Pkg.activate("../../../../")
Pkg.resolve()

using Zipper
using LinearAlgebra,Plots

power = 3
onsitepotential = 0
systemsize=2^power*3

power = power+1
scalingswifnoofrgstepsfor2 = [(2^i,power-i+1) for i in 1:power-1]
scalingswifnoofrgstepsfor3 = [(2^(i-1)*3,power-i+1) for i in 2:power-1]
scalingswifnoofrgsteps = vcat(scalingswifnoofrgstepsfor2,scalingswifnoofrgstepsfor3)

for (scaling,noofrgsteps) in scalingswifnoofrgsteps
    fiodir("../../../../../../../storage/data/slwongag/GMERA/dirac/systemsize$systemsize/localsize$scaling")
    for rgstep in range(1,noofrgsteps-1)
        blockedcorrelations = fioload("blockedcorrelations")
        distilled = distillation(fioload("gmera$rgstep"*"thirdcouriercorrelations")|>crystalspectrum, :filled=>(v->v < 0.5), :empty=>(v->v > 0.5))
        emptybands = distilled[:empty]
        emptyisometry = crystalisometry(emptybands)|>CrystalDenseMap
        couriercomposemap = fioload("gmera$rgstep"*"couriercomposemapgmera")
        terminatedcouriercorrelation = (couriercomposemap*emptyisometry)*(couriercomposemap*emptyisometry)'
        if rgstep == 1
            approximatecorrelation = fioload("gmera1approximatecorrelation")
        else
            approximatecorrelation = fioload("gmera$rgstep"*"approximatecorrelationsofar")
        end
        terminatedapproximatecorrelation = terminatedcouriercorrelation + approximatecorrelation
        fiosave(terminatedapproximatecorrelation, name="gmera$rgstep"*"terminatedapproximatecorrelation")

        terminatedalldiff = blockedcorrelations - terminatedapproximatecorrelation
        fiosave(terminatedalldiff, name="gmera$rgstep"*"terminatedalldiff")

        terminatedalldiffL1norm = focktraceL1norm(terminatedalldiff,systemsize^2*6)
        fiosave(terminatedalldiffL1norm, name="gmera$rgstep"*"terminatedalldiffL1norm")
    end
end
