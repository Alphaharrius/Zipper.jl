using Zipper
using LinearAlgebra
using Plots

setmaxthreads(8)

function focktraceL1norm(fockmap,systemsize)
    modeevalpairs = fockmap|>eigspech|>geteigenvalues
    return(sum([abs(modeevalpair[2]) for modeevalpair in modeevalpairs])/systemsize)
end

resultlist5 = []
[r for r in range(-1,-0.3,8)]

for onsitepotential in range(-1,-0.3,8)
    systemsize=40
    correlations,H = generatesystem(  onsitepotential+0im, -onsitepotential+0im,-1 + 0im,0im,systemsize)
    crystalfock = correlations|>getoutspace

    scale = Scale([5 0; 0 5], crystalfock|>getcrystal|>getspace)
    @info("Performing rgblocking...",scale)
    @info("Generating rgblocking transformation...")
    blocker = @time scale * crystalfock
    @info("Performing rgblocking on correlations...")
    blockedcorrelations = @time blocker * correlations * blocker'

    firstgmeracorrelations,firstglobalisometryfilled,firstglobalisometryempty,firstwanniercourierisometry = firstgmerastepbycount(blockedcorrelations,30)
    secondgmeracorrelations,secondglobalisometryfilled,secondglobalisometryempty,secondwanniercourierisometry = secondgmerastepbycount(firstgmeracorrelations,12,5)
    thirdglobalisometryfilled,thirdglobalisometryempty = thirdgmerastepbycount(secondgmeracorrelations,0,2)

    firstapproximation = firstglobalisometryempty*firstglobalisometryempty'
    secondapproximation = (firstwanniercourierisometry*secondglobalisometryempty)*(firstwanniercourierisometry*secondglobalisometryempty)'
    thirdapproximation = (firstwanniercourierisometry*secondwanniercourierisometry*thirdglobalisometryempty)*(firstwanniercourierisometry*secondwanniercourierisometry*thirdglobalisometryempty)'

    approximatecorrelation = firstapproximation + secondapproximation + thirdapproximation

    diff = blockedcorrelations-approximatecorrelation
    diffnorm = focktraceL1norm(diff,9600)
    append!(resultlist5,diffnorm)
end

scatter(resultlist5)
scatter!(resultlist[1:8])
scatter!(resultlist2[1:8])