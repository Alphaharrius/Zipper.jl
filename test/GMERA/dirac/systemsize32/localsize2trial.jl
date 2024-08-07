using Pkg
Pkg.resolve()
Pkg.Registry.update()
Pkg.instantiate()
Pkg.activate("../../../../")
Pkg.resolve()

using Zipper
using LinearAlgebra,Plots

setmaxthreads(setmaxthreads(Threads.nthreads()))
usecrystaldensemap()
# plotlyjs()

power = 5
onsitepotential = 0
nnhopping = 0
systemsize=2^power
correlations,H = generatesystem(-onsitepotential+0im, onsitepotential+0im,-1 + 0im,-nnhopping,systemsize)
crystalfock = correlations|>getoutspace

scaling = 2
fiodir("../../../../../../../storage/data/slwongag/GMERA/dirac/systemsize$systemsize/localsize$scaling/testingdifferentnoofdistillmodes")

blockedcorrelations,blockedH,blocker = blocking(correlations,H,scaling)
fiosave(blockedcorrelations, name="blockedcorrelations")
fiosave(blockedH, name="blockedH")
fiosave(blocker, name="blocker")

gmera1firststepdata = gmerafirststepbycount(blockedcorrelations,blockedH,384,1)
gmera1secondstepdata = gmerasecondstepbycount(gmera1firststepdata[:couriercorrelations],gmera1firststepdata[:courierH],96,(384/6)|>Int)


fiosave(gmera1firststepdata[:localcorrelations], name="gmera1firstlocalcorrelations")
fiosave(gmera1secondstepdata[:localcorrelations], name="gmera1secondlocalcorrelations")
