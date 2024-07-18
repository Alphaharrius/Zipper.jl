using Pkg
Pkg.resolve()
Pkg.Registry.update()
Pkg.instantiate()
Pkg.activate("../../../../../../../../")
Pkg.resolve()

using Zipper
using LinearAlgebra,Plots

setmaxthreads(8)

# plotlyjs()

function focktraceL1norm(fockmap,systemsize)
    modeevalpairs = fockmap|>eigspech|>geteigenvalues
    return(sum([abs(modeevalpair[2]) for modeevalpair in modeevalpairs])/systemsize)
end

function blocking(correlations,H,scaling)
    crystalfock = correlations|>getoutspace
    scale = Scale([scaling 0; 0 scaling], crystalfock|>getcrystal|>getspace)
    @info("Performing rgblocking...",scale)
    @info("Generating rgblocking transformation...")
    blocker = @time scale * crystalfock
    @info("Performing rgblocking on correlations...")
    blockedcorrelations = @time blocker * correlations * blocker'
    blockedH = @time blocker * H * blocker'
    return blockedcorrelations,blockedH,blocker
end

power = 5
onsitepotential = 0.2
systemsize=2^power
scaling = 16
correlations,H = generatesystem(-onsitepotential+0im, onsitepotential+0im,-1 + 0im,0im,systemsize)
crystalfock = correlations|>getoutspace

fiodir("../../../../../data/trivial/systemsize32/localsize16/onsitepotential_$onsitepotential")

blockedcorrelations,blockedH,blocer = blocking(correlations,H,scaling)
fiosave(blockedcorrelations, name="blockedcorrelations")
fiosave(blockedH, name="blockedH")

gmera1firststepdata = gmerafirststepbycount(blockedcorrelations,blockedH,1152,1)
gmera1secondstepdata = gmerasecondstepbycount(gmera1firststepdata[:couriercorrelations],gmera1firststepdata[:courierH],768,192)
gmera1thirdstepdata = gmerathirdstepbycount(gmera1secondstepdata[:couriercorrelations],gmera1secondstepdata[:courierH],384,128)
gmera1thirdstepterminatedata = gmerathirdstepbycount(gmera1secondstepdata[:couriercorrelations],gmera1secondstepdata[:courierH],0,128)

gmera1firstapproximation = gmera1firststepdata[:globalemptyisometry]*gmera1firststepdata[:globalemptyisometry]'
gmera1secondapproximation = (gmera1firststepdata[:wanniercourierisometry]*gmera1secondstepdata[:globalemptyisometry])*(gmera1firststepdata[:wanniercourierisometry]*gmera1secondstepdata[:globalemptyisometry])'
gmera1thirdapproximation = (gmera1firststepdata[:wanniercourierisometry]*gmera1secondstepdata[:wanniercourierisometry]*gmera1thirdstepdata[:globalemptyisometry])*(gmera1firststepdata[:wanniercourierisometry]*gmera1secondstepdata[:wanniercourierisometry]*gmera1thirdstepdata[:globalemptyisometry])'
gmera1thirdterminateapproximation = (gmera1firststepdata[:wanniercourierisometry]*gmera1secondstepdata[:wanniercourierisometry]*gmera1thirdstepterminatedata[:globalemptyisometry])*(gmera1firststepdata[:wanniercourierisometry]*gmera1secondstepdata[:wanniercourierisometry]*gmera1thirdstepterminatedata[:globalemptyisometry])'

couriercomposemapgmera1 = (gmera1firststepdata[:wanniercourierisometry]*gmera1secondstepdata[:wanniercourierisometry]*gmera1thirdstepdata[:wanniercourierisometry])

gmera1approximatecorrelation = gmera1firstapproximation + gmera1secondapproximation  + gmera1thirdapproximation
gmera1terminateapproximatecorrelation = gmera1firstapproximation + gmera1secondapproximation + gmera1thirdterminateapproximation
fiosave(gmera1approximatecorrelation, name="gmera1approximatecorrelation")
fiosave(gmera1terminateapproximatecorrelation, name="gmera1terminateapproximatecorrelation")

gmera1terminatediff = blockedcorrelations - gmera1terminateapproximatecorrelation
gmera1terminatediffnorm = focktraceL1norm(gmera1terminatediff,systemsize^2*6)
fiosave(gmera1terminatediffnorm, name="gmera1terminatetraceL1normofdiff")

fiosave(gmera1firststepdata[:localcorrelations], name="gmera1firstlocalcorrelations")
fiosave(gmera1firststepdata[:couriercorrelations], name="gmera1firstcouriercorrelations")
fiosave(gmera1firststepdata[:filledcorrelations], name="gmera1firstfilledcorrelations")
fiosave(gmera1firststepdata[:emptycorrelations], name="gmera1firstemptycorrelations")
fiosave(gmera1firststepdata[:courierH], name="gmera1firstcourierH")
fiosave(gmera1firststepdata[:filledH], name="gmera1firstfilledH")
fiosave(gmera1firststepdata[:emptyH], name="gmera1firstemptyH")
fiosave(gmera1firststepdata[:rawcouriercorrelations], name="gmera1firstrawcouriercorrelations")
fiosave(gmera1firststepdata[:wanniercourierisometry], name="gmera1firstwanniercourierisometry")
fiosave(gmera1firststepdata[:globalemptyisometry], name="gmera1firstglobalemptyisometry")
fiosave(gmera1firststepdata[:globalfilledisometry], name="gmera1firstglobalfilledisometry")
fiosave(gmera1firststepdata[:wanniercourierstates], name="gmera1firstwanniercourierstates")
fiosave(gmera1firststepdata[:minsvdforcourier], name="gmera1firstminsvdforcourier")

fiosave(gmera1secondstepdata[:localcorrelations], name="gmera1secondlocalcorrelations")
fiosave(gmera1secondstepdata[:couriercorrelations], name="gmera1secondcouriercorrelations")
fiosave(gmera1secondstepdata[:filledcorrelations], name="gmera1secondfilledcorrelations")
fiosave(gmera1secondstepdata[:emptycorrelations], name="gmera1secondemptycorrelations")
fiosave(gmera1secondstepdata[:courierH], name="gmera1secondcourierH")
fiosave(gmera1secondstepdata[:filledH], name="gmera1secondfilledH")
fiosave(gmera1secondstepdata[:emptyH], name="gmera1secondemptyH")
fiosave(gmera1secondstepdata[:rawcouriercorrelations], name="gmera1secondrawcouriercorrelations")
fiosave(gmera1secondstepdata[:wanniercourierisometry], name="gmera1secondwanniercourierisometry")
fiosave(gmera1secondstepdata[:globalemptyisometry], name="gmera1secondglobalemptyisometry")
fiosave(gmera1secondstepdata[:globalfilledisometry], name="gmera1secondglobalfilledisometry")
fiosave(gmera1secondstepdata[:wanniercourierstates], name="gmera1secondwanniercourierstates")
fiosave(gmera1secondstepdata[:minsvdforcourier], name="gmera1secondminsvdforcourier")

fiosave(gmera1thirdstepdata[:localcorrelations], name="gmera1thirdlocalcorrelations")
fiosave(gmera1thirdstepdata[:couriercorrelations], name="gmera1thirdcouriercorrelations")
fiosave(gmera1thirdstepdata[:filledcorrelations], name="gmera1thirdfilledcorrelations")
fiosave(gmera1thirdstepdata[:emptycorrelations], name="gmera1thirdemptycorrelations")
fiosave(gmera1thirdstepdata[:courierH], name="gmera1thirdcourierH")
fiosave(gmera1thirdstepdata[:filledH], name="gmera1thirdfilledH")
fiosave(gmera1thirdstepdata[:emptyH], name="gmera1thirdemptyH")
fiosave(gmera1thirdstepdata[:rawcouriercorrelations], name="gmera1thirdrawcouriercorrelations")
fiosave(gmera1thirdstepdata[:wanniercourierisometry], name="gmera1thirdwanniercourierisometry")
fiosave(gmera1thirdstepdata[:globalemptyisometry], name="gmera1thirdglobalemptyisometry")
fiosave(gmera1thirdstepdata[:globalfilledisometry], name="gmera1thirdglobalfilledisometry")
fiosave(gmera1thirdstepdata[:wanniercourierstates], name="gmera1thirdwanniercourierstates")
fiosave(gmera1thirdstepdata[:minsvdforcourier], name="gmera1thirdminsvdforcourier")

rg1correlations = gmera1thirdstepdata[:couriercorrelations]
rg1H = gmera1thirdstepdata[:courierH]
rg1blockedcorrelations,rg1blockedH,rg1blocker = blocking(rg1correlations,rg1H,2)
fiosave(rg1blockedcorrelations, name="rg1blockedcorrelations")
fiosave(rg1blockedH, name="rg1blockedH")
couriercomposemapgmera1 = couriercomposemapgmera1*rg1blocker'

gmera2finalstepdata = gmerafinalstep(rg1blockedcorrelations,rg1blockedH,64)
fiosave(gmera2finalstepdata[:localcorrelations], name="gmera2finallocalcorrelations")
fiosave(gmera2finalstepdata[:filledcorrelations], name="gmera2finalfilledcorrelations")
fiosave(gmera2finalstepdata[:emptycorrelations], name="gmera2finalemptycorrelations")
fiosave(gmera2finalstepdata[:filledH], name="gmera2finalfilledH")
fiosave(gmera2finalstepdata[:emptyH], name="gmera2finalemptyH")
fiosave(gmera2finalstepdata[:globalemptyisometry], name="gmera2finalglobalemptyisometry")
fiosave(gmera2finalstepdata[:globalfilledisometry], name="gmera2finalglobalfilledisometry")

gmera2finalapproximation = (couriercomposemapgmera1*gmera2finalstepdata[:globalemptyisometry])*(couriercomposemapgmera1*gmera2finalstepdata[:globalemptyisometry])'
fiosave(gmera2finalapproximation, name="gmera2finalapproximation")
gmera12allapproximatecorrelation = gmera1approximatecorrelation+gmera2finalapproximation
fiosave(gmera12allapproximatecorrelation, name="gmera12allapproximatecorrelation")

gmera12alldiff = blockedcorrelations - gmera12allapproximatecorrelation
gmera12alldiffnorm = focktraceL1norm(gmera12alldiff,systemsize^2*6)
fiosave(gmera12alldiffnorm, name="gmera12alldiffnorm")

