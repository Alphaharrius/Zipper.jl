using Pkg
Pkg.resolve()
Pkg.Registry.update()
Pkg.instantiate()
Pkg.activate("../../../../../")
Pkg.resolve()

using Zipper
using LinearAlgebra,Plots

power = 3
onsitepotential = 0
nnhopping = 0.3im
systemsize=2^power*3

power = power+1
scalingswifnoofrgstepsfor2 = [(2^i,power-i+1) for i in 1:power-1]
scalingswifnoofrgstepsfor3 = [(2^(i-1)*3,power-i+1) for i in 2:power-2]
scalingswifnoofrgsteps = vcat(scalingswifnoofrgstepsfor2,scalingswifnoofrgstepsfor3)

for (scaling,noofrgsteps) in scalingswifnoofrgsteps
    fiodir("../../../../../../../../storage/data/slwongag/GMERA/chern/systemsize$systemsize/nnhopping$nnhopping-localsize$scaling")
    gmeraalldiff = fioload("gmeraalldiff")
    gmeraalldiffL2norm = matrixL2norm(gmeraalldiff,systemsize)
    fiosave(gmeraalldiffL2norm, name="gmeraalldiffL2norm")
end