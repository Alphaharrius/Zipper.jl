using Zipper
using LinearAlgebra

setmaxthreads(8)

function focktraceL1norm(fockmap,systemsize)
    modeevalpairs = fockmap|>eigspech|>geteigenvalues
    return(sum([abs(modeevalpair[2]) for modeevalpair in modeevalpairs])/systemsize)
end

systemsize=64
correlations,H = generatesystem( 0im, 0im,-1 + 0im,0im,systemsize)
crystalfock = correlations|>getoutspace

scale = Scale([4 0; 0 4], crystalfock|>getcrystal|>getspace)
@info("Performing rgblocking...",scale)
@info("Generating rgblocking transformation...")
blocker = @time scale * crystalfock
@info("Performing rgblocking on correlations...")
blockedcorrelations = @time blocker * correlations * blocker'

rg1firstgmerastepdata = firstgmerastepbythreshold(blockedcorrelations,0.0002)
rg1firstgmerastepdata[:nooflocalcouriermodes]
rg1secondgmerastepdata = secondgmerastepbythreshold(rg1firstgmerastepdata[:couriercorrelations],0.002,div(rg1firstgmerastepdata[:nooflocalcouriermodes],6))
thirdgmeracorrelations,thirdglobalisometryfilled,thirdglobalisometryempty,thirdwanniercourierisometry = thirdgmerastepbycount(rg1secondgmerastepdata[:couriercorrelations],6,div(rg1secondgmerastepdata[:nooflocalcouriermodes],6))

rgedcrystalfock = thirdgmeracorrelations|>getoutspace
scale = Scale([4 0; 0 4], rgedcrystalfock |>getcrystal|>getspace)
@info("Performing rgblocking...",scale)
@info("Generating rgblocking transformation...")
rgedblocker = @time scale * rgedcrystalfock
@info("Performing rgblocking on correlations...")
rgedblockedcorrelations = @time rgedblocker * thirdgmeracorrelations * rgedblocker'

rgedfirstgmeracorrelations,rgedfirstglobalisometryfilled,rgedfirstglobalisometryempty,rgedfirstwanniercourierisometry = firstgmerastep(rgedblockedcorrelations,48)
rgedsecondgmeracorrelations,rgedsecondglobalisometryfilled,rgedsecondglobalisometryempty,rgedsecondwanniercourierisometry = secondgmerastep(rgedfirstgmeracorrelations,12,8)
rgedthirdgmeracorrelations,rgedthirdglobalisometryfilled,rgedthirdglobalisometryempty,rgedthirdwanniercourierisometry = thirdgmerastep(rgedsecondgmeracorrelations,6,2)