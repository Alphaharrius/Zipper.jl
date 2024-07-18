using Pkg
Pkg.resolve()
Pkg.Registry.update()
Pkg.instantiate()
Pkg.activate("../../../../../../../../")
Pkg.resolve()

using Zipper
using LinearAlgebra,Plots

setmaxthreads(8)
usecrystaldensemap()
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

function firstgmerastep(correlations,H,scaling)
    blockedcorrelations,blockedH,blocer = blocking(correlations,H,scaling)
    fiosave(blockedcorrelations, name="blockedcorrelations")
    fiosave(blockedH, name="blockedH")

    gmera1firststepdata = gmerafirststepbycount(blockedcorrelations,blockedH,18,1)
    gmera1secondstepdata = gmerasecondstepbycount(gmera1firststepdata[:couriercorrelations],gmera1firststepdata[:courierH],12,3)
    gmera1thirdstepdata = gmerathirdstepbycount(gmera1secondstepdata[:couriercorrelations],gmera1secondstepdata[:courierH],6,2)

    gmera1firstapproximation = gmera1firststepdata[:globalemptyisometry]*gmera1firststepdata[:globalemptyisometry]'
    gmera1secondapproximation = (gmera1firststepdata[:wanniercourierisometry]*gmera1secondstepdata[:globalemptyisometry])*(gmera1firststepdata[:wanniercourierisometry]*gmera1secondstepdata[:globalemptyisometry])'
    gmera1thirdapproximation = (gmera1firststepdata[:wanniercourierisometry]*gmera1secondstepdata[:wanniercourierisometry]*gmera1thirdstepdata[:globalemptyisometry])*(gmera1firststepdata[:wanniercourierisometry]*gmera1secondstepdata[:wanniercourierisometry]*gmera1thirdstepdata[:globalemptyisometry])'

    couriercomposemapgmera1 = (gmera1firststepdata[:wanniercourierisometry]*gmera1secondstepdata[:wanniercourierisometry]*gmera1thirdstepdata[:wanniercourierisometry])

    gmera1approximatecorrelation = gmera1firstapproximation + gmera1secondapproximation + gmera1thirdapproximation
    fiosave(gmera1approximatecorrelation, name="gmera1approximatecorrelation")

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
    return  rg1H,rg1correlations,couriercomposemapgmera1,gmera1approximatecorrelation,blockedcorrelations
end

function intermediategmerastep(rgcorrelations,rgH,couriercomposemap,gmerprevsumapproximatecorrelation,rgstep)
    rgblockedcorrelations,rgblockedH,rgblocker = blocking(rgcorrelations,rgH,2)
    prevrgstep = rgstep-1
    fiosave(rgblockedcorrelations, name="rg$prevrgstep"*"blockedcorrelations")
    fiosave(rgblockedH, name="rg$prevrgstep"*"blockedH")
    couriercomposemap = couriercomposemap*rgblocker'

    gmerafirststepdata = gmerafirststepbycount(rgblockedcorrelations,rgblockedH,18,1)
    gmerasecondstepdata = gmerasecondstepbycount(gmerafirststepdata[:couriercorrelations],gmerafirststepdata[:courierH],12,3)
    gmerathirdstepdata = gmerathirdstepbycount(gmerasecondstepdata[:couriercorrelations],gmerasecondstepdata[:courierH],6,2)

    
    gmerafirstapproximation = (couriercomposemap*gmerafirststepdata[:globalemptyisometry])*(couriercomposemap*gmerafirststepdata[:globalemptyisometry])'
    gmerasecondapproximation = (couriercomposemap*gmerafirststepdata[:wanniercourierisometry]*gmerasecondstepdata[:globalemptyisometry])*(couriercomposemap*gmerafirststepdata[:wanniercourierisometry]*gmerasecondstepdata[:globalemptyisometry])'
    gmerathirdapproximation = (couriercomposemap*gmerafirststepdata[:wanniercourierisometry]*gmerasecondstepdata[:wanniercourierisometry]*gmerathirdstepdata[:globalemptyisometry])*(couriercomposemap*gmerafirststepdata[:wanniercourierisometry]*gmerasecondstepdata[:wanniercourierisometry]*gmerathirdstepdata[:globalemptyisometry])'
    
    couriercomposemapgmera = couriercomposemap*(gmerafirststepdata[:wanniercourierisometry]*gmerasecondstepdata[:wanniercourierisometry]*gmerathirdstepdata[:wanniercourierisometry])

    gmeraapproximatecorrelation = gmerafirstapproximation + gmerasecondapproximation  + gmerathirdapproximation
    fiosave(gmeraapproximatecorrelation, name="gmera$rgstep"*"approximatecorrelation")

    fiosave(gmerafirststepdata[:localcorrelations], name="gmera$rgstep"*"firstlocalcorrelations")
    fiosave(gmerafirststepdata[:couriercorrelations], name="gmera$rgstep"*"firstcouriercorrelations")
    fiosave(gmerafirststepdata[:filledcorrelations], name="gmera$rgstep"*"firstfilledcorrelations")
    fiosave(gmerafirststepdata[:emptycorrelations], name="gmera$rgstep"*"firstemptycorrelations")
    fiosave(gmerafirststepdata[:courierH], name="gmera$rgstep"*"firstcourierH")
    fiosave(gmerafirststepdata[:filledH], name="gmera$rgstep"*"firstfilledH")
    fiosave(gmerafirststepdata[:emptyH], name="gmera$rgstep"*"firstemptyH")
    fiosave(gmerafirststepdata[:rawcouriercorrelations], name="gmera$rgstep"*"firstrawcouriercorrelations")
    fiosave(gmerafirststepdata[:wanniercourierisometry], name="gmera$rgstep"*"firstwanniercourierisometry")
    fiosave(gmerafirststepdata[:globalemptyisometry], name="gmera$rgstep"*"firstglobalemptyisometry")
    fiosave(gmerafirststepdata[:globalfilledisometry], name="gmera$rgstep"*"firstglobalfilledisometry")
    fiosave(gmerafirststepdata[:wanniercourierstates], name="gmera$rgstep"*"firstwanniercourierstates")
    fiosave(gmerafirststepdata[:minsvdforcourier], name="gmera$rgstep"*"firstminsvdforcourier")

    fiosave(gmerasecondstepdata[:localcorrelations], name="gmera$rgstep"*"secondlocalcorrelations")
    fiosave(gmerasecondstepdata[:couriercorrelations], name="gmera$rgstep"*"secondcouriercorrelations")
    fiosave(gmerasecondstepdata[:filledcorrelations], name="gmera$rgstep"*"secondfilledcorrelations")
    fiosave(gmerasecondstepdata[:emptycorrelations], name="gmera$rgstep"*"secondemptycorrelations")
    fiosave(gmerasecondstepdata[:courierH], name="gmera$rgstep"*"secondcourierH")
    fiosave(gmerasecondstepdata[:filledH], name="gmera$rgstep"*"secondfilledH")
    fiosave(gmerasecondstepdata[:emptyH], name="gmera$rgstep"*"secondemptyH")
    fiosave(gmerasecondstepdata[:rawcouriercorrelations], name="gmera$rgstep"*"secondrawcouriercorrelations")
    fiosave(gmerasecondstepdata[:wanniercourierisometry], name="gmera$rgstep"*"secondwanniercourierisometry")
    fiosave(gmerasecondstepdata[:globalemptyisometry], name="gmera$rgstep"*"secondglobalemptyisometry")
    fiosave(gmerasecondstepdata[:globalfilledisometry], name="gmera$rgstep"*"secondglobalfilledisometry")
    fiosave(gmerasecondstepdata[:wanniercourierstates], name="gmera$rgstep"*"secondwanniercourierstates")
    fiosave(gmerasecondstepdata[:minsvdforcourier], name="gmera$rgstep"*"secondminsvdforcourier")

    fiosave(gmerathirdstepdata[:localcorrelations], name="gmera$rgstep"*"thirdlocalcorrelations")
    fiosave(gmerathirdstepdata[:couriercorrelations], name="gmera$rgstep"*"thirdcouriercorrelations")
    fiosave(gmerathirdstepdata[:filledcorrelations], name="gmera$rgstep"*"thirdfilledcorrelations")
    fiosave(gmerathirdstepdata[:emptycorrelations], name="gmera$rgstep"*"thirdemptycorrelations")
    fiosave(gmerathirdstepdata[:courierH], name="gmera$rgstep"*"thirdcourierH")
    fiosave(gmerathirdstepdata[:filledH], name="gmera$rgstep"*"thirdfilledH")
    fiosave(gmerathirdstepdata[:emptyH], name="gmera$rgstep"*"thirdemptyH")
    fiosave(gmerathirdstepdata[:rawcouriercorrelations], name="gmera$rgstep"*"thirdrawcouriercorrelations")
    fiosave(gmerathirdstepdata[:wanniercourierisometry], name="gmera$rgstep"*"thirdwanniercourierisometry")
    fiosave(gmerathirdstepdata[:globalemptyisometry], name="gmera$rgstep"*"thirdglobalemptyisometry")
    fiosave(gmerathirdstepdata[:globalfilledisometry], name="gmera$rgstep"*"thirdglobalfilledisometry")
    fiosave(gmerathirdstepdata[:wanniercourierstates], name="gmera$rgstep"*"thirdwanniercourierstates")
    fiosave(gmerathirdstepdata[:minsvdforcourier], name="gmera$rgstep"*"thirdminsvdforcourier")

    rgcorrelations = gmerathirdstepdata[:couriercorrelations]
    rgH = gmerathirdstepdata[:courierH]   
    gmerasumapproximatecorrelation = gmerprevsumapproximatecorrelation+gmeraapproximatecorrelation
    return  rgH,rgcorrelations,couriercomposemapgmera,gmerasumapproximatecorrelation
end

function finalgmerastep(rgcorrelations,rgH,couriercomposemap,gmerprevsumapproximatecorrelation,blockedcorrelations,rgstep)
    rgblockedcorrelations,rgblockedH,rgblocker = blocking(rgcorrelations,rgH,2)
    prevrgstep = rgstep-1
    fiosave(rgblockedcorrelations, name="rg$prevrgstep"*"blockedcorrelations")
    fiosave(rgblockedH, name="rg$prevrgstep"*"blockedH")
    couriercomposemap = couriercomposemap*rgblocker'

    gmerafinalstepdata = gmerafinalstep(rgblockedcorrelations,rgblockedH,1)
    fiosave(gmerafinalstepdata[:localcorrelations], name="gmera$rgstep"*"finallocalcorrelations")
    fiosave(gmerafinalstepdata[:filledcorrelations], name="gmera$rgstep"*"finalfilledcorrelations")
    fiosave(gmerafinalstepdata[:emptycorrelations], name="gmera$rgstep"*"finalemptycorrelations")
    fiosave(gmerafinalstepdata[:filledH], name="gmera$rgstep"*"finalfilledH")
    fiosave(gmerafinalstepdata[:emptyH], name="gmera$rgstep"*"finalemptyH")
    fiosave(gmerafinalstepdata[:globalemptyisometry], name="gmera$rgstep"*"finalglobalemptyisometry")
    fiosave(gmerafinalstepdata[:globalfilledisometry], name="gmera$rgstep"*"finalglobalfilledisometry")

    gmerafinalapproximation = (couriercomposemap*gmerafinalstepdata[:globalemptyisometry])*(couriercomposemap*gmerafinalstepdata[:globalemptyisometry])'
    fiosave(gmerafinalapproximation, name="gmera$rgstep"*"finalapproximation")
    gmeraallapproximatecorrelation = gmerprevsumapproximatecorrelation+gmerafinalapproximation
    fiosave(gmeraallapproximatecorrelation, name="gmeraallapproximatecorrelation")

    gmeraalldiff = blockedcorrelations - gmeraallapproximatecorrelation
    gmeraalldiff = gmeraalldiff|>CrystalFockMap
    gmeraalldiffnorm = focktraceL1norm(gmeraalldiff,systemsize^2*6)
    fiosave(gmeraalldiffnorm, name="gmeraalldiffnorm")
end

power = 5
onsitepotential = 0
nnhopping = 0.4im
systemsize=2^power
scaling = 2
correlations,H = generatesystem(-onsitepotential+0im, onsitepotential+0im,-1 + 0im,-nnhopping,systemsize)
crystalfock = correlations|>getoutspace

fiodir("../../../../../data/chern/systemsize$systemsize/localsize2/nnhopping_$nnhopping")

rg1H,rg1correlations,courier1composemap,gmera1approximatecorrelation,blockedcorrelations = firstgmerastep(correlations,H,scaling)
rg2H,rg2correlations,courier12composemap,gmera12approximatecorrelation = intermediategmerastep(rg1correlations,rg1H,courier1composemap,gmera1approximatecorrelation,2)
rg3H,rg3correlations,courier123composemap,gmera123approximatecorrelation = intermediategmerastep(rg2correlations,rg2H,courier12composemap,gmera12approximatecorrelation,3)
rg4H,rg4correlations,courier1234composemap,gmera1234approximatecorrelation = intermediategmerastep(rg3correlations,rg3H,courier123composemap,gmera123approximatecorrelation,4)
finalgmerastep(rg4correlations,rg4H,courier1234composemap,gmera1234approximatecorrelation,blockedcorrelations,5)

