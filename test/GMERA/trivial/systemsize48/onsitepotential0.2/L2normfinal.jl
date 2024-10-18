using Pkg
Pkg.resolve()
Pkg.Registry.update()
Pkg.instantiate()
Pkg.activate("../../../../../")
Pkg.resolve()

using Zipper
using LinearAlgebra,Plots

power = 4
onsitepotential = 0.2
nnhopping = 0
systemsize=2^power*3

power = power+1
scalingswifnoofrgstepsfor2 = [(2^i,power-i+1) for i in 1:power-2]
scalingswifnoofrgstepsfor3 = [(2^(i-1)*3,power-i+1) for i in 2:power-3]
scalingswifnoofrgsteps = vcat(scalingswifnoofrgstepsfor2,scalingswifnoofrgstepsfor3)

for (scaling,noofrgsteps) in scalingswifnoofrgsteps[4:4]
    fiodir("../../../../../../../../storage/data/slwongag/GMERA/trivial/systemsize$systemsize/onsitepotential$onsitepotential-localsize$scaling")
    gmeraalldiff = fioload("gmeraalldiff")
    gmeraalldiffL2norm = matrixL2normmod(gmeraalldiff,systemsize)
    fiosave(gmeraalldiffL2norm, name="gmeraalldiffL2norm")
end