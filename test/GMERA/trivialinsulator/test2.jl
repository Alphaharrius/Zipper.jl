using Zipper
using LinearAlgebra

setmaxthreads(8)

function focktraceL1norm(fockmap,systemsize)
    modeevalpairs = fockmap|>eigspech|>geteigenvalues
    return(sum([abs(modeevalpair[2]) for modeevalpair in modeevalpairs])/systemsize)
end

systemsize=40
correlations,H = generatesystem( -0.3+0im, 0.3+0im,-1 + 0im,0im,systemsize)
crystalfock = correlations|>getoutspace

scale = Scale([2 0; 0 2], crystalfock|>getcrystal|>getspace)
@info("Performing rgblocking...",scale)
@info("Generating rgblocking transformation...")
blocker = @time scale * crystalfock
@info("Performing rgblocking on correlations...")
blockedcorrelations = @time blocker * correlations * blocker'

rg1firstgmerastepdata = firstgmerastepbythreshold(blockedcorrelations,0.05)
rg1secondgmerastepdata = secondgmerastepbythreshold(rg1firstgmerastepdata[:couriercorrelations],0.05,div(rg1firstgmerastepdata[:nooflocalcouriermodes],6))
rg1thirdgmerastepdata = thirdgmerastepbythreshold(rg1secondgmerastepdata[:couriercorrelations],0.05,div(rg1secondgmerastepdata[:nooflocalcouriermodes],6))

firstapproximation = rg1firstgmerastepdata[:globalemptyisometry]*rg1firstgmerastepdata[:globalemptyisometry]'
secondapproximation = (rg1firstgmerastepdata[:wanniercourierisometry]*rg1secondgmerastepdata[:globalemptyisometry])*(rg1firstgmerastepdata[:wanniercourierisometry]*rg1secondgmerastepdata[:globalemptyisometry])'
thirdapproximation = (rg1firstgmerastepdata[:wanniercourierisometry]*rg1secondgmerastepdata[:wanniercourierisometry]*rg1thirdgmerastepdata[:globalemptyisometry])*(rg1firstgmerastepdata[:wanniercourierisometry]*rg1secondgmerastepdata[:wanniercourierisometry]*rg1thirdgmerastepdata[:globalemptyisometry])'

arrpoximatecorrelationtwo = firstapproximation + secondapproximation + thirdapproximation

difftwo = blockedcorrelations-arrpoximatecorrelationtwo
difftwonorm = focktraceL1norm(difftwo,9600)