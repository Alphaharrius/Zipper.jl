using Pkg
Pkg.resolve()
Pkg.Registry.update()
Pkg.instantiate()
Pkg.activate("../../../../../../../")
Pkg.resolve()

using Zipper
using LinearAlgebra,Plots

setmaxthreads(8)
usecrystaldensemap()
# plotlyjs()

power = 5
onsitepotential = 0.1
nnhopping = 0
systemsize=2^power
scaling = 4
correlations,H = generatesystem(-onsitepotential+0im, onsitepotential+0im,-1 + 0im,-nnhopping,systemsize)
crystalfock = correlations|>getoutspace
# fiodir("/Users/slwongag/Desktop/data/trivial/systemsize32/localsize4/onsitepotential_0.1")
fiodir("../../../../data/trivial/systemsize$systemsize/localsize4/onsitepotential_$onsitepotential")

rg1H,rg1correlations,courier1composemap,gmera1approximatecorrelation,blockedcorrelations,noofflavourpermodeforintermediate = firstgmerastep(correlations,H,scaling)
rg2H,rg2correlations,courier12composemap,gmera12approximatecorrelation = intermediategmerastep(rg1correlations,rg1H,courier1composemap,gmera1approximatecorrelation,2,noofflavourpermodeforintermediate)
rg3H,rg3correlations,courier123composemap,gmera123approximatecorrelation = intermediategmerastep(rg2correlations,rg2H,courier12composemap,gmera12approximatecorrelation,3,noofflavourpermodeforintermediate)
finalgmerastep(rg3correlations,rg3H,courier123composemap,gmera123approximatecorrelation,blockedcorrelations,4,noofflavourpermodeforintermediate,systemsize)

