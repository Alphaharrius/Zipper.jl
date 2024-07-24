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
onsitepotential = 0
nnhopping = 0.4im
systemsize=2^power
scaling = 2
correlations,H = generatesystem(-onsitepotential+0im, onsitepotential+0im,-1 + 0im,-nnhopping,systemsize)
crystalfock = correlations|>getoutspace

fiodir("../../../../data/chern/systemsize$systemsize/localsize2/nnhopping_$nnhopping")

rg1H,rg1correlations,courier1composemap,gmera1approximatecorrelation,blockedcorrelations,noofflavourpermodeforintermediate = firstgmerastep(correlations,H,scaling)
rg2H,rg2correlations,courier12composemap,gmera12approximatecorrelation = intermediategmerastep(rg1correlations,rg1H,courier1composemap,gmera1approximatecorrelation,2,noofflavourpermodeforintermediate)
rg3H,rg3correlations,courier123composemap,gmera123approximatecorrelation = intermediategmerastep(rg2correlations,rg2H,courier12composemap,gmera12approximatecorrelation,3,noofflavourpermodeforintermediate)
rg4H,rg4correlations,courier1234composemap,gmera1234approximatecorrelation = intermediategmerastep(rg3correlations,rg3H,courier123composemap,gmera123approximatecorrelation,4,noofflavourpermodeforintermediate)
finalgmerastep(rg4correlations,rg4H,courier1234composemap,gmera1234approximatecorrelation,blockedcorrelations,5,noofflavourpermodeforintermediate,systemsize)

