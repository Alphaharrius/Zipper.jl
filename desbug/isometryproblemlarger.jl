using Zipper
using LinearAlgebra,Plots


setmaxthreads(setmaxthreads(Threads.nthreads()))
usecrystaldensemap()
plotlyjs()

power = 5
onsitepotential = 0
nnhopping = 0
systemsize=2^power
correlations,H = generatesystem(-onsitepotential+0im, onsitepotential+0im,-1 + 0im,-nnhopping,systemsize)
crystalfock = correlations|>getoutspace
crystal = crystalfock|>getcrystal

scaling = 8

blockedcorrelations,blockedH,blocker = blocking(correlations,H,scaling)

firstgmeraresult = gmerafirststep(blockedcorrelations,blockedH)
secondgmeraresult = gmerasecondstep(firstgmeraresult[:firstgmeracorrelations],firstgmeraresult[:firstgmeraH])
thirdgmeraresult = gmerathirdstep(secondgmeraresult[:secondgmeracorrelations],secondgmeraresult[:secondgmeraH])

thirdgmeraresult[:thirdgmeracorrelations]|>getinspace|>getcrystal|>getunitcell|>visualize
rg2thirdgmeraresult[:thirdgmeracorrelations]|>getinspace|>getcrystal|>getunitcell|>visualize

subscaling = 2

blockedcorrelations2,blockedH2,blocker2 = blocking(thirdgmeraresult[:thirdgmeracorrelations],thirdgmeraresult[:thirdgmeraH],subscaling)
rg2firstgmeraresult = gmerafirststep(blockedcorrelations2,blockedH2)
rg2secondgmeraresult = gmerasecondstep(rg2firstgmeraresult[:firstgmeracorrelations],rg2firstgmeraresult[:firstgmeraH])
rg2thirdgmeraresult = gmerathirdstep(rg2secondgmeraresult[:secondgmeracorrelations],rg2secondgmeraresult[:secondgmeraH])

blockedcorrelations-blockedcorrelations
sum(abs.((blockedcorrelations|>CrystalFockMap)|>rep))

firstgmeraresult
secondgmeraresult[:secondgmeracorrelations]
thirdgmeraresult[:rawcouriercorrelations]

blockedcorrelations3,blockedH3,blocker3 = blocking(rg2thirdgmeraresult[:thirdgmeracorrelations],rg2thirdgmeraresult[:thirdgmeraH],subscaling)
rg3firstgmeraresult = gmerafirststep(blockedcorrelations3,blockedH3)
rg3secondgmeraresult = gmerasecondstep(rg3firstgmeraresult[:firstgmeracorrelations],rg3firstgmeraresult[:firstgmeraH])
rg3thirdgmeraresult = gmerathirdstep(rg3secondgmeraresult[:secondgmeracorrelations],rg3secondgmeraresult[:secondgmeraH])

blockedcorrelations3,blockedH3,blocker3 = blocking(rg2thirdgmeraresult[:thirdgmeracorrelations],rg2thirdgmeraresult[:thirdgmeraH],subscaling)
gmerafinalstep(blockedcorrelations3,blockedH3)

