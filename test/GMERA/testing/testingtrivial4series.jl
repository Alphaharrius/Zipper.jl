using Zipper
using LinearAlgebra
using Plots

setmaxthreads(8)

function focktraceL1norm(fockmap,systemsize)
    modeevalpairs = fockmap|>eigspech|>geteigenvalues
    return(sum([abs(modeevalpair[2]) for modeevalpair in modeevalpairs])/systemsize)
end

resultlist = []
[r for r in range(-0.4,-0.1,4)]

for onsitepotential in range(-0.4,-0.1,4)
    systemsize=40
    correlations,H = generatesystem(  onsitepotential+0im, -onsitepotential+0im,-1 + 0im,0im,systemsize)
    crystalfock = correlations|>getoutspace

    scale = Scale([4 0; 0 4], crystalfock|>getcrystal|>getspace)
    @info("Performing rgblocking...",scale)
    @info("Generating rgblocking transformation...")
    blocker = @time scale * crystalfock
    @info("Performing rgblocking on correlations...")
    blockedcorrelations = @time blocker * correlations * blocker'

    firstgmeracorrelations,firstglobalisometryfilled,firstglobalisometryempty,firstwanniercourierisometry = firstgmerastepbycount(blockedcorrelations,24)
    secondgmeracorrelations,secondglobalisometryfilled,secondglobalisometryempty,secondwanniercourierisometry = secondgmerastepbycount(firstgmeracorrelations,12,4)
    thirdglobalisometryfilled,thirdglobalisometryempty = thirdgmerastepbycount(secondgmeracorrelations,0,2)

    firstapproximation = firstglobalisometryempty*firstglobalisometryempty'
    secondapproximation = (firstwanniercourierisometry*secondglobalisometryempty)*(firstwanniercourierisometry*secondglobalisometryempty)'
    thirdapproximation = (firstwanniercourierisometry*secondwanniercourierisometry*thirdglobalisometryempty)*(firstwanniercourierisometry*secondwanniercourierisometry*thirdglobalisometryempty)'

    approximatecorrelation = firstapproximation + secondapproximation + thirdapproximation

    diff = blockedcorrelations-approximatecorrelation
    diffnorm = focktraceL1norm(diff,9600)
    append!(resultlist,diffnorm)
end
resultlist
scatter!(resultlist)
data = [[resultlist2[r],resultlist[r],resultlist5[r]] for r in range(1,length(resultlist5))]
scatter([2,4,5],log.(data[1]))
scatter!([2,4,5],log.(data[2]))
scatter!([2,4,5],log.(data[3]))
scatter!([2,4,5],log.(data[4]))
scatter!([2,4,5],log.(data[5]))
scatter!([2,4,5],log.(data[6]))
scatter!([2,4,5],log.(data[7]))
scatter!([2,4,5],log.(data[8]))


for p in data[2:end]
    scatter!(p)
end