using Zipper
using LinearAlgebra
using Plots

setmaxthreads(8)

function focktraceL1norm(fockmap,systemsize)
    modeevalpairs = fockmap|>eigspech|>geteigenvalues
    return(sum([abs(modeevalpair[2]) for modeevalpair in modeevalpairs])/systemsize)
end

resultlist2 = []
[r for r in range(-0.4,-0.1,4)]
for onsitepotential in range(-1,-0.5,6)
    systemsize=40
    correlations,H = generatesystem(  onsitepotential+0im, -onsitepotential+0im,-1 + 0im,0im,systemsize)
    crystalfock = correlations|>getoutspace

    scale = Scale([2 0; 0 2], crystalfock|>getcrystal|>getspace)
    @info("Performing rgblocking...",scale)
    @info("Generating rgblocking transformation...")
    blocker = @time scale * crystalfock
    @info("Performing rgblocking on correlations...")
    blockedcorrelations = @time blocker * correlations * blocker'

    firstgmeracorrelations,firstglobalisometryfilled,firstglobalisometryempty,firstwanniercourierisometry = firstgmerastepbycount(blockedcorrelations,12)
    secondgmeracorrelations,secondglobalisometryfilled,secondglobalisometryempty,secondwanniercourierisometry = secondgmerastepbycount(firstgmeracorrelations,6,2)
    thirdglobalisometryfilled,thirdglobalisometryempty = thirdgmerastepbycount(secondgmeracorrelations,0,1)

    firstapproximation = firstglobalisometryempty*firstglobalisometryempty'
    secondapproximation = (firstwanniercourierisometry*secondglobalisometryempty)*(firstwanniercourierisometry*secondglobalisometryempty)'
    thirdapproximation = (firstwanniercourierisometry*secondwanniercourierisometry*thirdglobalisometryempty)*(firstwanniercourierisometry*secondwanniercourierisometry*thirdglobalisometryempty)'

    approximatecorrelation = firstapproximation + secondapproximation + thirdapproximation

    diff = blockedcorrelations-approximatecorrelation
    diffnorm = focktraceL1norm(diff,9600)
    append!(resultlist2,diffnorm)
end

scatter(resultlist2)