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
scaling = 2
correlations,H = generatesystem(-onsitepotential+0im, onsitepotential+0im,-1 + 0im,0im,systemsize)
crystalfock = correlations|>getoutspace

fiodir("../../../../../data/trivial/systemsize32/localsize2/onsitepotential_$onsitepotential")

blockedcorrelations,blockedH,blocer = blocking(correlations,H,scaling)
fiosave(blockedcorrelations, name="blockedcorrelations")
fiosave(blockedH, name="blockedH")

gmera1firststepdata = gmerafirststepbycount(blockedcorrelations,blockedH,18,1)
gmera1secondstepdata = gmerasecondstepbycount(gmera1firststepdata[:couriercorrelations],gmera1firststepdata[:courierH],12,3)
gmera1thirdstepdata = gmerathirdstepbycount(gmera1secondstepdata[:couriercorrelations],gmera1secondstepdata[:courierH],6,2)
gmera1thirdstepterminatedata = gmerathirdstepbycount(gmera1secondstepdata[:couriercorrelations],gmera1secondstepdata[:courierH],0,2)

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

gmera2firststepdata = gmerafirststepbycount(rg1blockedcorrelations,rg1blockedH,18,1)
gmera2secondstepdata = gmerasecondstepbycount(gmera2firststepdata[:couriercorrelations],gmera2firststepdata[:courierH],12,3)
gmera2thirdstepdata = gmerathirdstepbycount(gmera2secondstepdata[:couriercorrelations],gmera2secondstepdata[:courierH],6,2)

gmera2thirdstepterminatedata = gmerathirdstepbycount(gmera2secondstepdata[:couriercorrelations],gmera2secondstepdata[:courierH],0,2)

gmera2firstapproximation = (couriercomposemapgmera1*gmera2firststepdata[:globalemptyisometry])*(couriercomposemapgmera1*gmera2firststepdata[:globalemptyisometry])'
gmera2secondapproximation = (couriercomposemapgmera1*gmera2firststepdata[:wanniercourierisometry]*gmera2secondstepdata[:globalemptyisometry])*(couriercomposemapgmera1*gmera2firststepdata[:wanniercourierisometry]*gmera2secondstepdata[:globalemptyisometry])'
gmera2thirdapproximation = (couriercomposemapgmera1*gmera2firststepdata[:wanniercourierisometry]*gmera2secondstepdata[:wanniercourierisometry]*gmera2thirdstepdata[:globalemptyisometry])*(couriercomposemapgmera1*gmera2firststepdata[:wanniercourierisometry]*gmera2secondstepdata[:wanniercourierisometry]*gmera2thirdstepdata[:globalemptyisometry])'
gmera2thirdterminateapproximation = (couriercomposemapgmera1*gmera2firststepdata[:wanniercourierisometry]*gmera2secondstepdata[:wanniercourierisometry]*gmera2thirdstepterminatedata[:globalemptyisometry])*(couriercomposemapgmera1*gmera2firststepdata[:wanniercourierisometry]*gmera2secondstepdata[:wanniercourierisometry]*gmera2thirdstepterminatedata[:globalemptyisometry])'

couriercomposemapgmera12 = couriercomposemapgmera1*(gmera2firststepdata[:wanniercourierisometry]*gmera2secondstepdata[:wanniercourierisometry]*gmera2thirdstepdata[:wanniercourierisometry])

gmera2approximatecorrelation = gmera2firstapproximation + gmera2secondapproximation  + gmera2thirdapproximation
gmera12approximatecorrelation = gmera1approximatecorrelation+gmera2approximatecorrelation
gmera12terminateapproximatecorrelation = gmera1approximatecorrelation+ gmera2firstapproximation + gmera2secondapproximation + gmera2thirdterminateapproximation
fiosave(gmera2approximatecorrelation, name="gmera2approximatecorrelation")
fiosave(gmera12terminateapproximatecorrelation, name="gmera12terminateapproximatecorrelation")

gmera12terminatediff = blockedcorrelations - gmera12terminateapproximatecorrelation
gmera12terminatediffnorm = focktraceL1norm(gmera12terminatediff,systemsize^2*6)
fiosave(gmera12terminatediffnorm, name="gmera12terminatetraceL1normofdiff")

fiosave(gmera2firststepdata[:localcorrelations], name="gmera2firstlocalcorrelations")
fiosave(gmera2firststepdata[:couriercorrelations], name="gmera2firstcouriercorrelations")
fiosave(gmera2firststepdata[:filledcorrelations], name="gmera2firstfilledcorrelations")
fiosave(gmera2firststepdata[:emptycorrelations], name="gmera2firstemptycorrelations")
fiosave(gmera2firststepdata[:courierH], name="gmera2firstcourierH")
fiosave(gmera2firststepdata[:filledH], name="gmera2firstfilledH")
fiosave(gmera2firststepdata[:emptyH], name="gmera2firstemptyH")
fiosave(gmera2firststepdata[:rawcouriercorrelations], name="gmera2firstrawcouriercorrelations")
fiosave(gmera2firststepdata[:wanniercourierisometry], name="gmera2firstwanniercourierisometry")
fiosave(gmera2firststepdata[:globalemptyisometry], name="gmera2firstglobalemptyisometry")
fiosave(gmera2firststepdata[:globalfilledisometry], name="gmera2firstglobalfilledisometry")
fiosave(gmera2firststepdata[:wanniercourierstates], name="gmera2firstwanniercourierstates")
fiosave(gmera2firststepdata[:minsvdforcourier], name="gmera2firstminsvdforcourier")

fiosave(gmera2secondstepdata[:localcorrelations], name="gmera2secondlocalcorrelations")
fiosave(gmera2secondstepdata[:couriercorrelations], name="gmera2secondcouriercorrelations")
fiosave(gmera2secondstepdata[:filledcorrelations], name="gmera2secondfilledcorrelations")
fiosave(gmera2secondstepdata[:emptycorrelations], name="gmera2secondemptycorrelations")
fiosave(gmera2secondstepdata[:courierH], name="gmera2secondcourierH")
fiosave(gmera2secondstepdata[:filledH], name="gmera2secondfilledH")
fiosave(gmera2secondstepdata[:emptyH], name="gmera2secondemptyH")
fiosave(gmera2secondstepdata[:rawcouriercorrelations], name="gmera2secondrawcouriercorrelations")
fiosave(gmera2secondstepdata[:wanniercourierisometry], name="gmera2secondwanniercourierisometry")
fiosave(gmera2secondstepdata[:globalemptyisometry], name="gmera2secondglobalemptyisometry")
fiosave(gmera2secondstepdata[:globalfilledisometry], name="gmera2secondglobalfilledisometry")
fiosave(gmera2secondstepdata[:wanniercourierstates], name="gmera2secondwanniercourierstates")
fiosave(gmera2secondstepdata[:minsvdforcourier], name="gmera2secondminsvdforcourier")

fiosave(gmera2thirdstepdata[:localcorrelations], name="gmera2thirdlocalcorrelations")
fiosave(gmera2thirdstepdata[:couriercorrelations], name="gmera2thirdcouriercorrelations")
fiosave(gmera2thirdstepdata[:filledcorrelations], name="gmera2thirdfilledcorrelations")
fiosave(gmera2thirdstepdata[:emptycorrelations], name="gmera2thirdemptycorrelations")
fiosave(gmera2thirdstepdata[:courierH], name="gmera2thirdcourierH")
fiosave(gmera2thirdstepdata[:filledH], name="gmera2thirdfilledH")
fiosave(gmera2thirdstepdata[:emptyH], name="gmera2thirdemptyH")
fiosave(gmera2thirdstepdata[:rawcouriercorrelations], name="gmera2thirdrawcouriercorrelations")
fiosave(gmera2thirdstepdata[:wanniercourierisometry], name="gmera2thirdwanniercourierisometry")
fiosave(gmera2thirdstepdata[:globalemptyisometry], name="gmera2thirdglobalemptyisometry")
fiosave(gmera2thirdstepdata[:globalfilledisometry], name="gmera2thirdglobalfilledisometry")
fiosave(gmera2thirdstepdata[:wanniercourierstates], name="gmera2thirdwanniercourierstates")
fiosave(gmera2thirdstepdata[:minsvdforcourier], name="gmera2thirdminsvdforcourier")

rg2correlations = gmera2thirdstepdata[:couriercorrelations]
rg2H = gmera2thirdstepdata[:courierH]
rg2blockedcorrelations,rg2blockedH,rg2blocker = blocking(rg2correlations,rg2H,2)
fiosave(rg2blockedcorrelations, name="rg2blockedcorrelations")
fiosave(rg2blockedH, name="rg2blockedH")
couriercomposemapgmera12 = couriercomposemapgmera12*rg2blocker'

gmera3firststepdata = gmerafirststepbycount(rg2blockedcorrelations,rg2blockedH,18,1)
gmera3secondstepdata = gmerasecondstepbycount(gmera3firststepdata[:couriercorrelations],gmera3firststepdata[:courierH],12,3)
gmera3thirdstepdata = gmerathirdstepbycount(gmera3secondstepdata[:couriercorrelations],gmera3secondstepdata[:courierH],6,2)

gmera3thirdstepterminatedata = gmerathirdstepbycount(gmera3secondstepdata[:couriercorrelations],gmera3secondstepdata[:courierH],0,2)

gmera3firstapproximation = (couriercomposemapgmera12*gmera3firststepdata[:globalemptyisometry])*(couriercomposemapgmera12*gmera3firststepdata[:globalemptyisometry])'
gmera3secondapproximation = (couriercomposemapgmera12*gmera3firststepdata[:wanniercourierisometry]*gmera3secondstepdata[:globalemptyisometry])*(couriercomposemapgmera12*gmera3firststepdata[:wanniercourierisometry]*gmera3secondstepdata[:globalemptyisometry])'
gmera3thirdapproximation = (couriercomposemapgmera12*gmera3firststepdata[:wanniercourierisometry]*gmera3secondstepdata[:wanniercourierisometry]*gmera3thirdstepdata[:globalemptyisometry])*(couriercomposemapgmera12*gmera3firststepdata[:wanniercourierisometry]*gmera3secondstepdata[:wanniercourierisometry]*gmera3thirdstepdata[:globalemptyisometry])'
gmera3thirdterminateapproximation = (couriercomposemapgmera12*gmera3firststepdata[:wanniercourierisometry]*gmera3secondstepdata[:wanniercourierisometry]*gmera3thirdstepterminatedata[:globalemptyisometry])*(couriercomposemapgmera12*gmera3firststepdata[:wanniercourierisometry]*gmera3secondstepdata[:wanniercourierisometry]*gmera3thirdstepterminatedata[:globalemptyisometry])'

couriercomposemapgmera123 = couriercomposemapgmera12*(gmera3firststepdata[:wanniercourierisometry]*gmera3secondstepdata[:wanniercourierisometry]*gmera3thirdstepdata[:wanniercourierisometry])

gmera3approximatecorrelation = gmera3firstapproximation + gmera3secondapproximation  + gmera3thirdapproximation
gmera123approximatecorrelation = gmera1approximatecorrelation+gmera2approximatecorrelation+gmera3approximatecorrelation
gmera123terminateapproximatecorrelation = gmera12approximatecorrelation+gmera3firstapproximation + gmera3secondapproximation + gmera3thirdterminateapproximation
fiosave(gmera3approximatecorrelation, name="gmera3approximatecorrelation")
fiosave(gmera123terminateapproximatecorrelation, name="gmera123terminateapproximatecorrelation")

gmera123terminatediff = blockedcorrelations - gmera123terminateapproximatecorrelation
gmera123terminatediffnorm = focktraceL1norm(gmera123terminatediff,systemsize^2*6)
fiosave(gmera123terminatediffnorm, name="gmera123terminatetraceL1normofdiff")


fiosave(gmera3firststepdata[:localcorrelations], name="gmera3firstlocalcorrelations")
fiosave(gmera3firststepdata[:couriercorrelations], name="gmera3firstcouriercorrelations")
fiosave(gmera3firststepdata[:filledcorrelations], name="gmera3firstfilledcorrelations")
fiosave(gmera3firststepdata[:emptycorrelations], name="gmera3firstemptycorrelations")
fiosave(gmera3firststepdata[:courierH], name="gmera3firstcourierH")
fiosave(gmera3firststepdata[:filledH], name="gmera3firstfilledH")
fiosave(gmera3firststepdata[:emptyH], name="gmera3firstemptyH")
fiosave(gmera3firststepdata[:rawcouriercorrelations], name="gmera3firstrawcouriercorrelations")
fiosave(gmera3firststepdata[:wanniercourierisometry], name="gmera3firstwanniercourierisometry")
fiosave(gmera3firststepdata[:globalemptyisometry], name="gmera3firstglobalemptyisometry")
fiosave(gmera3firststepdata[:globalfilledisometry], name="gmera3firstglobalfilledisometry")
fiosave(gmera3firststepdata[:wanniercourierstates], name="gmera3firstwanniercourierstates")
fiosave(gmera3firststepdata[:minsvdforcourier], name="gmera3firstminsvdforcourier")

fiosave(gmera3secondstepdata[:localcorrelations], name="gmera3secondlocalcorrelations")
fiosave(gmera3secondstepdata[:couriercorrelations], name="gmera3secondcouriercorrelations")
fiosave(gmera3secondstepdata[:filledcorrelations], name="gmera3secondfilledcorrelations")
fiosave(gmera3secondstepdata[:emptycorrelations], name="gmera3secondemptycorrelations")
fiosave(gmera3secondstepdata[:courierH], name="gmera3secondcourierH")
fiosave(gmera3secondstepdata[:filledH], name="gmera3secondfilledH")
fiosave(gmera3secondstepdata[:emptyH], name="gmera3secondemptyH")
fiosave(gmera3secondstepdata[:rawcouriercorrelations], name="gmera3secondrawcouriercorrelations")
fiosave(gmera3secondstepdata[:wanniercourierisometry], name="gmera3secondwanniercourierisometry")
fiosave(gmera3secondstepdata[:globalemptyisometry], name="gmera3secondglobalemptyisometry")
fiosave(gmera3secondstepdata[:globalfilledisometry], name="gmera3secondglobalfilledisometry")
fiosave(gmera3secondstepdata[:wanniercourierstates], name="gmera3secondwanniercourierstates")
fiosave(gmera3secondstepdata[:minsvdforcourier], name="gmera3secondminsvdforcourier")

fiosave(gmera3thirdstepdata[:localcorrelations], name="gmera3thirdlocalcorrelations")
fiosave(gmera3thirdstepdata[:couriercorrelations], name="gmera3thirdcouriercorrelations")
fiosave(gmera3thirdstepdata[:filledcorrelations], name="gmera3thirdfilledcorrelations")
fiosave(gmera3thirdstepdata[:emptycorrelations], name="gmera3thirdemptycorrelations")
fiosave(gmera3thirdstepdata[:courierH], name="gmera3thirdcourierH")
fiosave(gmera3thirdstepdata[:filledH], name="gmera3thirdfilledH")
fiosave(gmera3thirdstepdata[:emptyH], name="gmera3thirdemptyH")
fiosave(gmera3thirdstepdata[:rawcouriercorrelations], name="gmera3thirdrawcouriercorrelations")
fiosave(gmera3thirdstepdata[:wanniercourierisometry], name="gmera3thirdwanniercourierisometry")
fiosave(gmera3thirdstepdata[:globalemptyisometry], name="gmera3thirdglobalemptyisometry")
fiosave(gmera3thirdstepdata[:globalfilledisometry], name="gmera3thirdglobalfilledisometry")
fiosave(gmera3thirdstepdata[:wanniercourierstates], name="gmera3thirdwanniercourierstates")
fiosave(gmera3thirdstepdata[:minsvdforcourier], name="gmera3thirdminsvdforcourier")

rg3correlations = gmera3thirdstepdata[:couriercorrelations]
rg3H = gmera3thirdstepdata[:courierH]
rg3blockedcorrelations,rg3blockedH,rg3blocker = blocking(rg3correlations,rg3H,2)
fiosave(rg3blockedcorrelations, name="rg3blockedcorrelations")
fiosave(rg3blockedH, name="rg3blockedH")
couriercomposemapgmera123 = couriercomposemapgmera123*rg3blocker'

gmera4firststepdata = gmerafirststepbycount(rg3blockedcorrelations,rg3blockedH,18,1)
gmera4secondstepdata = gmerasecondstepbycount(gmera4firststepdata[:couriercorrelations],gmera4firststepdata[:courierH],12,3)
gmera4thirdstepdata = gmerathirdstepbycount(gmera4secondstepdata[:couriercorrelations],gmera4secondstepdata[:courierH],6,2)

gmera4thirdstepterminatedata = gmerathirdstepbycount(gmera4secondstepdata[:couriercorrelations],gmera4secondstepdata[:courierH],0,2)

gmera4firstapproximation = (couriercomposemapgmera123*gmera4firststepdata[:globalemptyisometry])*(couriercomposemapgmera123*gmera4firststepdata[:globalemptyisometry])'
gmera4secondapproximation = (couriercomposemapgmera123*gmera4firststepdata[:wanniercourierisometry]*gmera4secondstepdata[:globalemptyisometry])*(couriercomposemapgmera123*gmera4firststepdata[:wanniercourierisometry]*gmera4secondstepdata[:globalemptyisometry])'
gmera4thirdapproximation = (couriercomposemapgmera123*gmera4firststepdata[:wanniercourierisometry]*gmera4secondstepdata[:wanniercourierisometry]*gmera4thirdstepdata[:globalemptyisometry])*(couriercomposemapgmera123*gmera4firststepdata[:wanniercourierisometry]*gmera4secondstepdata[:wanniercourierisometry]*gmera4thirdstepdata[:globalemptyisometry])'
gmera4thirdterminateapproximation = (couriercomposemapgmera123*gmera4firststepdata[:wanniercourierisometry]*gmera4secondstepdata[:wanniercourierisometry]*gmera4thirdstepterminatedata[:globalemptyisometry])*(couriercomposemapgmera123*gmera4firststepdata[:wanniercourierisometry]*gmera4secondstepdata[:wanniercourierisometry]*gmera4thirdstepterminatedata[:globalemptyisometry])'

couriercomposemapgmera1234 = couriercomposemapgmera123*(gmera4firststepdata[:wanniercourierisometry]*gmera4secondstepdata[:wanniercourierisometry]*gmera4thirdstepdata[:wanniercourierisometry])

gmera4approximatecorrelation = gmera4firstapproximation + gmera4secondapproximation  + gmera4thirdapproximation
gmera1234approximatecorrelation = gmera1approximatecorrelation+gmera2approximatecorrelation+gmera3approximatecorrelation+gmera4approximatecorrelation
gmera1234terminateapproximatecorrelation = gmera123approximatecorrelation+gmera4firstapproximation + gmera4secondapproximation + gmera4thirdterminateapproximation
fiosave(gmera4approximatecorrelation, name="gmera4approximatecorrelation")
fiosave(gmera1234terminateapproximatecorrelation, name="gmera1234terminateapproximatecorrelation")

gmera1234terminatediff = blockedcorrelations - gmera1234terminateapproximatecorrelation
gmera1234terminatediffnorm = focktraceL1norm(gmera1234terminatediff,systemsize^2*6)
fiosave(gmera1234terminatediffnorm, name="gmera1234terminatetraceL1normofdiff")

fiosave(gmera4firststepdata[:localcorrelations], name="gmera4firstlocalcorrelations")
fiosave(gmera4firststepdata[:couriercorrelations], name="gmera4firstcouriercorrelations")
fiosave(gmera4firststepdata[:filledcorrelations], name="gmera4firstfilledcorrelations")
fiosave(gmera4firststepdata[:emptycorrelations], name="gmera4firstemptycorrelations")
fiosave(gmera4firststepdata[:courierH], name="gmera4firstcourierH")
fiosave(gmera4firststepdata[:filledH], name="gmera4firstfilledH")
fiosave(gmera4firststepdata[:emptyH], name="gmera4firstemptyH")
fiosave(gmera4firststepdata[:rawcouriercorrelations], name="gmera4firstrawcouriercorrelations")
fiosave(gmera4firststepdata[:wanniercourierisometry], name="gmera4firstwanniercourierisometry")
fiosave(gmera4firststepdata[:globalemptyisometry], name="gmera4firstglobalemptyisometry")
fiosave(gmera4firststepdata[:globalfilledisometry], name="gmera4firstglobalfilledisometry")
fiosave(gmera4firststepdata[:wanniercourierstates], name="gmera4firstwanniercourierstates")
fiosave(gmera4firststepdata[:minsvdforcourier], name="gmera4firstminsvdforcourier")

fiosave(gmera4secondstepdata[:localcorrelations], name="gmera4secondlocalcorrelations")
fiosave(gmera4secondstepdata[:couriercorrelations], name="gmera4secondcouriercorrelations")
fiosave(gmera4secondstepdata[:filledcorrelations], name="gmera4secondfilledcorrelations")
fiosave(gmera4secondstepdata[:emptycorrelations], name="gmera4secondemptycorrelations")
fiosave(gmera4secondstepdata[:courierH], name="gmera4secondcourierH")
fiosave(gmera4secondstepdata[:filledH], name="gmera4secondfilledH")
fiosave(gmera4secondstepdata[:emptyH], name="gmera4secondemptyH")
fiosave(gmera4secondstepdata[:rawcouriercorrelations], name="gmera4secondrawcouriercorrelations")
fiosave(gmera4secondstepdata[:wanniercourierisometry], name="gmera4secondwanniercourierisometry")
fiosave(gmera4secondstepdata[:globalemptyisometry], name="gmera4secondglobalemptyisometry")
fiosave(gmera4secondstepdata[:globalfilledisometry], name="gmera4secondglobalfilledisometry")
fiosave(gmera4secondstepdata[:wanniercourierstates], name="gmera4secondwanniercourierstates")
fiosave(gmera4secondstepdata[:minsvdforcourier], name="gmera4secondminsvdforcourier")

fiosave(gmera4thirdstepdata[:localcorrelations], name="gmera4thirdlocalcorrelations")
fiosave(gmera4thirdstepdata[:couriercorrelations], name="gmera4thirdcouriercorrelations")
fiosave(gmera4thirdstepdata[:filledcorrelations], name="gmera4thirdfilledcorrelations")
fiosave(gmera4thirdstepdata[:emptycorrelations], name="gmera4thirdemptycorrelations")
fiosave(gmera4thirdstepdata[:courierH], name="gmera4thirdcourierH")
fiosave(gmera4thirdstepdata[:filledH], name="gmera4thirdfilledH")
fiosave(gmera4thirdstepdata[:emptyH], name="gmera4thirdemptyH")
fiosave(gmera4thirdstepdata[:rawcouriercorrelations], name="gmera4thirdrawcouriercorrelations")
fiosave(gmera4thirdstepdata[:wanniercourierisometry], name="gmera4thirdwanniercourierisometry")
fiosave(gmera4thirdstepdata[:globalemptyisometry], name="gmera4thirdglobalemptyisometry")
fiosave(gmera4thirdstepdata[:globalfilledisometry], name="gmera4thirdglobalfilledisometry")
fiosave(gmera4thirdstepdata[:wanniercourierstates], name="gmera4thirdwanniercourierstates")
fiosave(gmera4thirdstepdata[:minsvdforcourier], name="gmera4thirdminsvdforcourier")


rg4correlations = gmera4thirdstepdata[:couriercorrelations]
rg4H = gmera4thirdstepdata[:courierH]
rg4blockedcorrelations,rg4blockedH,rg4blocker = blocking(rg4correlations,rg4H,2)
fiosave(rg4blockedcorrelations, name="rg4blockedcorrelations")
fiosave(rg4blockedH, name="rg4blockedH")
couriercomposemapgmera1234 = couriercomposemapgmera1234*rg4blocker'

gmera5finalstepdata = gmerafinalstep(rg4blockedcorrelations,rg4blockedH,1)
fiosave(gmera5finalstepdata[:localcorrelations], name="gmera5finallocalcorrelations")
fiosave(gmera5finalstepdata[:filledcorrelations], name="gmera5finalfilledcorrelations")
fiosave(gmera5finalstepdata[:emptycorrelations], name="gmera5finalemptycorrelations")
fiosave(gmera5finalstepdata[:filledH], name="gmera5finalfilledH")
fiosave(gmera5finalstepdata[:emptyH], name="gmera5finalemptyH")
fiosave(gmera5finalstepdata[:globalemptyisometry], name="gmera5finalglobalemptyisometry")
fiosave(gmera5finalstepdata[:globalfilledisometry], name="gmera5finalglobalfilledisometry")

gmera5finalapproximation = (couriercomposemapgmera1234*gmera5finalstepdata[:globalemptyisometry])*(couriercomposemapgmera1234*gmera5finalstepdata[:globalemptyisometry])'
fiosave(gmera5finalapproximation, name="gmera5finalapproximation")
gmera12345allapproximatecorrelation = gmera1approximatecorrelation+gmera2approximatecorrelation+gmera3approximatecorrelation+gmera4approximatecorrelation+gmera5finalapproximation
fiosave(gmera12345allapproximatecorrelation, name="gmera12345allapproximatecorrelation")

gmera12345alldiff = blockedcorrelations - gmera12345allapproximatecorrelation
gmera12345alldiffnorm = focktraceL1norm(gmera12345alldiff,systemsize^2*6)
fiosave(gmera12345alldiffnorm, name="gmera12345alldiffnorm")