# using Pkg
# Pkg.resolve()
# Pkg.instantiate()
# Pkg.activate("../../../../../../../")
# Pkg.resolve()
# Pkg.instantiate()

using Zipper
using LinearAlgebra

setmaxthreads(8)

function focktraceL1norm(fockmap,systemsize)
    modeevalpairs = fockmap|>eigspech|>geteigenvalues
    return(sum([abs(modeevalpair[2]) for modeevalpair in modeevalpairs])/systemsize)
end

for onsitepotential in range(-1,-0.1,10)
    fiodir("/Users/slwongag/Desktop/data/trivial/systemsize40/localsize4/scheme[54,18,0]/onsitepotential_$onsitepotential")
    systemsize=40
    correlations,H = generatesystem(  onsitepotential+0im, -onsitepotential+0im,-1 + 0im,0im,systemsize)
    crystalfock = correlations|>getoutspace

    scale = Scale([4 0; 0 4], crystalfock|>getcrystal|>getspace)
    @info("Performing rgblocking...",scale)
    @info("Generating rgblocking transformation...")
    blocker = @time scale * crystalfock
    @info("Performing rgblocking on correlations...")
    blockedcorrelations = @time blocker * correlations * blocker'

    gmera1firststepdata = gmerafirststepbycount(blockedcorrelations,54)
    gmera1secondstepdata = gmerasecondstepbycount(gmera1firststepdata[:couriercorrelations],18,9)
    gmera1thirdstepdata = gmerathirdstepbycount(gmera1secondstepdata[:couriercorrelations],0,3)

    fiosave(gmera1firststepdata[:couriercorrelations], name="gmera1firstcouriercorrelations")
    fiosave(gmera1firststepdata[:filledcorrelations], name="gmera1firstfilledcorrelations")
    fiosave(gmera1firststepdata[:emptycorrelations], name="gmera1firstemptycorrelations")
    fiosave(gmera1firststepdata[:rawcouriercorrelations], name="gmera1firstrawcouriercorrelations")
    fiosave(gmera1firststepdata[:wanniercourierisometry], name="gmera1firstwanniercourierisometry")
    fiosave(gmera1firststepdata[:globalemptyisometry], name="gmera1firstglobalemptyisometry")
    fiosave(gmera1firststepdata[:globalfilledisometry], name="gmera1firstglobalfilledisometry")
    fiosave(gmera1firststepdata[:wanniercourierstates], name="gmera1firstwanniercourierstates")
    fiosave(gmera1firststepdata[:minsvdforcourier], name="gmera1firstminsvdforcourier")

    fiosave(gmera1secondstepdata[:couriercorrelations], name="gmera1secondcouriercorrelations")
    fiosave(gmera1secondstepdata[:filledcorrelations], name="gmera1secondfilledcorrelations")
    fiosave(gmera1secondstepdata[:emptycorrelations], name="gmera1secondemptycorrelations")
    fiosave(gmera1secondstepdata[:rawcouriercorrelations], name="gmera1secondrawcouriercorrelations")
    fiosave(gmera1secondstepdata[:wanniercourierisometry], name="gmera1secondwanniercourierisometry")
    fiosave(gmera1secondstepdata[:globalemptyisometry], name="gmera1secondglobalemptyisometry")
    fiosave(gmera1secondstepdata[:globalfilledisometry], name="gmera1secondglobalfilledisometry")
    fiosave(gmera1secondstepdata[:wanniercourierstates], name="gmera1secondwanniercourierstates")
    fiosave(gmera1secondstepdata[:minsvdforcourier], name="gmera1secondminsvdforcourier")

    fiosave(gmera1thirdstepdata[:filledcorrelations], name="gmera1thirdfilledcorrelations")
    fiosave(gmera1thirdstepdata[:emptycorrelations], name="gmera1thirdemptycorrelations")
    fiosave(gmera1thirdstepdata[:globalemptyisometry], name="gmera1thirdglobalemptyisometry")
    fiosave(gmera1thirdstepdata[:globalfilledisometry], name="gmera1thirdglobalfilledisometry")

    firstapproximation = gmera1firststepdata[:globalemptyisometry]*gmera1firststepdata[:globalemptyisometry]'
    secondapproximation = (gmera1firststepdata[:wanniercourierisometry]*gmera1secondstepdata[:globalemptyisometry])*(gmera1firststepdata[:wanniercourierisometry]*gmera1secondstepdata[:globalemptyisometry])'
    thirdapproximation = (gmera1firststepdata[:wanniercourierisometry]*gmera1secondstepdata[:wanniercourierisometry]*gmera1thirdstepdata[:globalemptyisometry])*(gmera1firststepdata[:wanniercourierisometry]*gmera1secondstepdata[:wanniercourierisometry]*gmera1thirdstepdata[:globalemptyisometry])'

    approximatecorrelation = firstapproximation + secondapproximation + thirdapproximation

    diff = blockedcorrelations-approximatecorrelation
    diffnorm = focktraceL1norm(diff,systemsize^2*6)
    fiosave(diffnorm, name="tracenormofdiff")
end