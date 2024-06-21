using Zipper
using LinearAlgebra

systemsize=32
correlations,H = generatesystem( -0.3+0im, 0.3+0im,-1 + 0im,0im,systemsize)
crystalfock = correlations|>getoutspace

scale = Scale([4 0; 0 4], crystalfock|>getcrystal|>getspace)
@info("Performing rgblocking...",scale)
@info("Generating rgblocking transformation...")
blocker = @time scale * crystalfock
@info("Performing rgblocking on correlations...")
blockedcorrelations = @time blocker * correlations * blocker'

firstgmeracorrelations,firstglobalisometryfilled,firstglobalisometryempty = firstgmerastep(blockedcorrelations,48)
secondgmeracorrelations,secondglobalisometryfilled,secondglobalisometryempty = secondgmerastep(firstgmeracorrelations,12,8)
thirdglobalisometryfilled,thirdglobalisometryempty = thirdgmerastep(secondgmeracorrelations,0,2)

rgedcrystalfock = thirdgmeracorrelations|>getoutspace

scale = Scale([4 0; 0 4], rgedcrystalfock|>getcrystal|>getspace)
@info("Performing rgblocking...",scale)
@info("Generating rgblocking transformation...")
blocker = @time scale * rgedcrystalfock
@info("Performing rgblocking on correlations...")
blockedrgedcorrelations = @time blocker * thirdgmeracorrelations * blocker'

firstgmerargedcorrelations = firstgmerastep(blockedrgedcorrelations,48)
secondgmerargedcorrelations = secondgmerastep(firstgmerargedcorrelations,12,8)
thirdgmerargedcorrelations = thirdgmerastep(secondgmerargedcorrelations,6,2)