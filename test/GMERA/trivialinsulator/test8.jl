using Zipper
using LinearAlgebra
using Plots

setmaxthreads(8)

function focktraceL1norm(fockmap,systemsize)
    modeevalpairs = fockmap|>eigspech|>geteigenvalues
    return(sum([abs(modeevalpair[2]) for modeevalpair in modeevalpairs])/systemsize)
end

function focktraceL2norm(fockmap,volume)
    @info("Calculating L2norm...")
    return (real(sqrt(tr(fockmap*fockmap'|>rep)/volume)))
end

systemsize=32
correlations,H = generatesystem( -0.3+0im, 0.3+0im,-1 + 0im,0im,systemsize)
crystalfock = correlations|>getoutspace

scale = Scale([8 0; 0 8], crystalfock|>getcrystal|>getspace)
@info("Performing rgblocking...",scale)
@info("Generating rgblocking transformation...")
blocker = @time scale * crystalfock
@info("Performing rgblocking on correlations...")
blockedcorrelations = @time blocker * correlations * blocker'

rg1firstgmerastepdata = firstgmerastepbythreshold(blockedcorrelations,0.004)
rg1secondgmerastepdata = secondgmerastepbythreshold(rg1firstgmerastepdata[:couriercorrelations],0.004,div(rg1firstgmerastepdata[:nooflocalcouriermodes],6))
rg1thirdgmerastepdata = thirdgmerastepbythreshold(rg1secondgmerastepdata[:couriercorrelations],0.004,div(rg1secondgmerastepdata[:nooflocalcouriermodes],6))

firstapproximation = rg1firstgmerastepdata[:globalemptyisometry]*rg1firstgmerastepdata[:globalemptyisometry]'
secondapproximation = (rg1firstgmerastepdata[:wanniercourierisometry]*rg1secondgmerastepdata[:globalemptyisometry])*(rg1firstgmerastepdata[:wanniercourierisometry]*rg1secondgmerastepdata[:globalemptyisometry])'
thirdapproximation = (rg1firstgmerastepdata[:wanniercourierisometry]*rg1secondgmerastepdata[:wanniercourierisometry]*rg1thirdgmerastepdata[:globalemptyisometry])*(rg1firstgmerastepdata[:wanniercourierisometry]*rg1secondgmerastepdata[:wanniercourierisometry]*rg1thirdgmerastepdata[:globalemptyisometry])'

arrpoximatecorrelationeightprime = firstapproximation + secondapproximation + thirdapproximation

diffeightprime = blockedcorrelations-arrpoximatecorrelationeightprime
diffeightnormprime = focktraceL1norm(diffeightprime,6144)
focktraceL2norm(diffeightprime,6144)
diffeightnormprime

errorlist = log.([difftwonorm,difffournorm,difffivenorm,diffeightnorm])

# radius48 = ([8/(2^2*2^2*2),24/(2^2*3^2*2),48/(2^2*4^2*2),120/(2^2*6^2*2),224/(2^2*8^2*2),528/(2^2*12^2*2)])
radius = ([2,4,5,8])
Plots.scatter(radius,errorlist ,mode="markers")