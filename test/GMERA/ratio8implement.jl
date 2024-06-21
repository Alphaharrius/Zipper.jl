using Zipper
using LinearAlgebra

systemsize=32
correlations,H = generatesystem( 0im, 0im,-1 + 0im,0.4im,systemsize)
crystalfock = correlations|>getoutspace

scale = Scale([8 0; 0 8], crystalfock|>getcrystal|>getspace)
@info("Performing rgblocking...",scale)
@info("Generating rgblocking transformation...")
blocker = @time scale * crystalfock
@info("Performing rgblocking on correlations...")
blockedcorrelations = @time blocker * correlations * blocker'

firstgmeracorrelations = firstgmerastep(blockedcorrelations,96)
secondgmeracorrelations = secondgmerastep(firstgmeracorrelations,24,16)
thirdgmeracorrelations = thirdgmerastep(secondgmeracorrelations,6,4)

rgedcrystalfock = thirdgmeracorrelations|>getoutspace

scale = Scale([2 0; 0 2], rgedcrystalfock|>getcrystal|>getspace)
@info("Performing rgblocking...",scale)
@info("Generating rgblocking transformation...")
blocker = @time scale * rgedcrystalfock
@info("Performing rgblocking on correlations...")
blockedrgedcorrelations = @time blocker * thirdgmeracorrelations * blocker'

firstgmerargedcorrelations = firstgmerastep(blockedrgedcorrelations,18)
secondgmerargedcorrelations = secondgmerastep(firstgmerargedcorrelations,12,3)
thirdgmerargedcorrelations = thirdgmerastep(secondgmerargedcorrelations,6,2)

rgedcrystalfock = thirdgmerargedcorrelations|>getoutspace

scale = Scale([2 0; 0 2], rgedcrystalfock|>getcrystal|>getspace)
@info("Performing rgblocking...",scale)
@info("Generating rgblocking transformation...")
blocker = @time scale * rgedcrystalfock
@info("Performing rgblocking on correlations...")
blockedrgedcorrelations = @time blocker * thirdgmerargedcorrelations * blocker'

