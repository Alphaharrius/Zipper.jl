using Zipper
using LinearAlgebra

setmaxthreads(8)

function focktraceL1norm(fockmap,systemsize)
    modeevalpairs = fockmap|>eigspech|>geteigenvalues
    return(sum([abs(modeevalpair[2]) for modeevalpair in modeevalpairs])/systemsize)
end

systemsize=40
correlations,H = generatesystem(  -0.3+0im, 0.3+0im,-1 + 0im,0im,systemsize)
crystalfock = correlations|>getoutspace

scale = Scale([4 0; 0 4], crystalfock|>getcrystal|>getspace)
@info("Performing rgblocking...",scale)
@info("Generating rgblocking transformation...")
blocker = @time scale * crystalfock
@info("Performing rgblocking on correlations...")
blockedcorrelations = @time blocker * correlations * blocker'

results = Dict()
for threshold in range(0.002,0.05,49)
    @info("implementing GMERA step for threshold ", threshold)
    rg1firstgmerastepdata = firstgmerastepbythreshold(blockedcorrelations,threshold)
    rg1secondgmerastepdata = secondgmerastepbythreshold(rg1firstgmerastepdata[:couriercorrelations],threshold,div(rg1firstgmerastepdata[:nooflocalcouriermodes],6))
    rg1thirdgmerastepdata = thirdgmerastepbythreshold(rg1secondgmerastepdata[:couriercorrelations],threshold,div(rg1secondgmerastepdata[:nooflocalcouriermodes],6))
    if haskey(rg1thirdgmerastepdata,:nooflocalcouriermodes)
        @info("doesn't terminate")
    else
        firstapproximation = rg1firstgmerastepdata[:globalemptyisometry]*rg1firstgmerastepdata[:globalemptyisometry]'
        secondapproximation = (rg1firstgmerastepdata[:wanniercourierisometry]*rg1secondgmerastepdata[:globalemptyisometry])*(rg1firstgmerastepdata[:wanniercourierisometry]*rg1secondgmerastepdata[:globalemptyisometry])'
        thirdapproximation = (rg1firstgmerastepdata[:wanniercourierisometry]*rg1secondgmerastepdata[:wanniercourierisometry]*rg1thirdgmerastepdata[:globalemptyisometry])*(rg1firstgmerastepdata[:wanniercourierisometry]*rg1secondgmerastepdata[:wanniercourierisometry]*rg1thirdgmerastepdata[:globalemptyisometry])'

        arrpoximatecorrelationfour = firstapproximation+secondapproximation+thirdapproximation

        difffour = blockedcorrelations-arrpoximatecorrelationfour
        difffournorm = focktraceL1norm(difffour,9600)
        nooflocalcouriermodespersteps = [rg1firstgmerastepdata[:nooflocalcouriermodes],rg1secondgmerastepdata[:nooflocalcouriermodes]]
        results[threshold]=>(difffournorm,nooflocalcouriermodespersteps)
    end
end

haskey(rg1thirdgmerastepdata,:nooflocalcouriermodes)
firstapproximation = rg1firstgmerastepdata[:globalemptyisometry]*rg1firstgmerastepdata[:globalemptyisometry]'
secondapproximation = (rg1firstgmerastepdata[:wanniercourierisometry]*rg1secondgmerastepdata[:globalemptyisometry])*(rg1firstgmerastepdata[:wanniercourierisometry]*rg1secondgmerastepdata[:globalemptyisometry])'
thirdapproximation = (rg1firstgmerastepdata[:wanniercourierisometry]*rg1secondgmerastepdata[:wanniercourierisometry]*rg1thirdgmerastepdata[:globalemptyisometry])*(rg1firstgmerastepdata[:wanniercourierisometry]*rg1secondgmerastepdata[:wanniercourierisometry]*rg1thirdgmerastepdata[:globalemptyisometry])'

arrpoximatecorrelationfour = firstapproximation+secondapproximation+thirdapproximation

difffour = blockedcorrelations-arrpoximatecorrelationfour
difffournorm = focktraceL1norm(difffour,9600)

difffournorm