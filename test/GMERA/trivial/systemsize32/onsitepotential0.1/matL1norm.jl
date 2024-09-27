using Pkg
Pkg.resolve()
Pkg.Registry.update()
Pkg.instantiate()
Pkg.activate("../../../../../")
Pkg.resolve()

using Zipper
using LinearAlgebra,Plots

setmaxthreads(setmaxthreads(Threads.nthreads()))
usecrystaldensemap()

power = 5
onsitepotential = 0.1
systemsize=2^power

scalingswifnoofrgsteps = [(2^i,power-i+1) for i in 1:power-1]

for (scaling,noofrgsteps) in scalingswifnoofrgsteps
    fiodir("../../../../../../../../storage/data/slwongag/GMERA/trivial/systemsize$systemsize/onsitepotential$onsitepotential-localsize$scaling")
    gmeraalldiff = fioload("gmeraalldiff")
    gmeraalldiffmatL1norm =  focktraceL1norm(gmeraalldiff,systemsize^2*6)
    fiosave(gmeraalldiffmatL1norm, name="gmeraalldiffmatL1norm")
end