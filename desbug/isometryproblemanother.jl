using Zipper
using LinearAlgebra,Plots

setmaxthreads(setmaxthreads(Threads.nthreads()))
usecrystaldensemap()
plotlyjs()

power = 3
onsitepotential = 0.05
nnhopping = 0
systemsize=2^power*3
correlations,H = generatesystem(-onsitepotential+0im, onsitepotential+0im,-1 + 0im,-nnhopping,systemsize)
crystalfock = correlations|>getoutspace
crystal = crystalfock|>getcrystal

scaling = 6

rg1H,rg1correlations,couriercomposemapgmera1,gmera1approximatecorrelation,blockedcorrelations,noofflavourpermodeforlaterrg = firstgmerastep(correlations,H,scaling)
rg2H,rg2correlations,couriercomposemap,gmeraapproximatecorrelationsofar = intermediategmerastep(rg1correlations,rg1H,couriercomposemapgmera1,gmera1approximatecorrelation,2,noofflavourpermodeforlaterrg)
rg3H,rg3correlations,couriercomposemap,gmeraapproximatecorrelationsofar = intermediategmerastep(rg2correlations,rg2H,couriercomposemap,gmeraapproximatecorrelationsofar,3,noofflavourpermodeforlaterrg)
finalgmerastep(rg1correlations,rg1H,couriercomposemapgmera1,gmera1approximatecorrelation,blockedcorrelations,4,24,noofflavourpermodeforlaterrg)

2: 12.313461843669518
4: 8.183303176393963
6: 7.782090380593459
8: 2.523282106354092
12: 3.362318674751047

rg2correlations|>getinspace|>getcrystal|>size