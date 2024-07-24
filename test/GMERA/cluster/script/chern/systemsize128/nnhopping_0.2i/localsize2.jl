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

power = 7
onsitepotential = 0
nnhopping = 0.2im
systemsize=2^power
scaling = 2
correlations,H = generatesystem(-onsitepotential+0im, onsitepotential+0im,-1 + 0im,-nnhopping,systemsize)
crystalfock = correlations|>getoutspace

fiodir("../../../../data/chern/systemsize$systemsize/localsize2/nnhopping_$nnhopping")

rg1H,rg1correlations,courier1composemap,gmera1approximatecorrelation,blockedcorrelations,noofflavourpermodeforintermediate = firstgmerastep(correlations,H,scaling)
rg2H,rg2correlations,courier12composemap,gmera12approximatecorrelation = intermediategmerastep(rg1correlations,rg1H,courier1composemap,gmera1approximatecorrelation,2,noofflavourpermodeforintermediate)
rg3H,rg3correlations,courier123composemap,gmera123approximatecorrelation = intermediategmerastep(rg2correlations,rg2H,courier12composemap,gmera12approximatecorrelation,3,noofflavourpermodeforintermediate)
rg4H,rg4correlations,courier1234composemap,gmera1234approximatecorrelation = intermediategmerastep(rg3correlations,rg3H,courier123composemap,gmera123approximatecorrelation,4,noofflavourpermodeforintermediate)
rg5H,rg5correlations,courier12345composemap,gmera12345approximatecorrelation = intermediategmerastep(rg4correlations,rg4H,courier1234composemap,gmera1234approximatecorrelation,5,noofflavourpermodeforintermediate)
rg6H,rg6correlations,courier123456composemap,gmera123456approximatecorrelation = intermediategmerastep(rg5correlations,rg5H,courier12345composemap,gmera12345approximatecorrelation,6,noofflavourpermodeforintermediate)
finalgmerastep(rg6correlations,rg6H,courier123456composemap,gmera123456approximatecorrelation,blockedcorrelations,7,noofflavourpermodeforintermediate,systemsize)

