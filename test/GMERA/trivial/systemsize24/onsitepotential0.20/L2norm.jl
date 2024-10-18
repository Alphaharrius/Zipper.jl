using Pkg
Pkg.resolve()
Pkg.Registry.update()
Pkg.instantiate()
Pkg.activate("../../../../../")
Pkg.resolve()

using Zipper
using LinearAlgebra,Plots

power = 3
onsitepotential = 0.2
nnhopping = 0
systemsize=2^power*3

power = power+1
scalingswifnoofrgstepsfor2 = [(2^i,power-i+1) for i in 1:power-1]
scalingswifnoofrgstepsfor3 = [(2^(i-1)*3,power-i+1) for i in 2:power-2]
scalingswifnoofrgsteps = vcat(scalingswifnoofrgstepsfor2,scalingswifnoofrgstepsfor3)

for (scaling,noofrgsteps) in scalingswifnoofrgsteps
    fiodir("../../../../../../../../storage/data/slwongag/GMERA/trivial/systemsize$systemsize/onsitepotential$onsitepotential-localsize$scaling")
    for rgstep in range(1,noofrgsteps-1)
        gmerafirstalldiff = fioload("gmera$rgstep"*"firstalldiff")
        gmerasecondalldiff = fioload("gmera$rgstep"*"secondalldiff")
        gmerathirdalldiff = fioload("gmera$rgstep"*"thirdalldiff")

        gmerafirstalldiffL2norm = matrixL2norm(gmerafirstalldiff,systemsize)
        gmerasecondalldiffL2norm = matrixL2norm(gmerasecondalldiff,systemsize)
        gmerathirdalldiffL2norm = matrixL2norm(gmerathirdalldiff,systemsize)

        fiosave(gmerafirstalldiffL2norm, name="gmera$rgstep"*"firstalldiffL2norm")
        fiosave(gmerasecondalldiffL2norm, name="gmera$rgstep"*"secondalldiffL2norm")
        fiosave(gmerathirdalldiffL2norm, name="gmera$rgstep"*"thirdalldiffL2norm")
    end
end