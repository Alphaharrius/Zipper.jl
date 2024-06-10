using Zipper
using LinearAlgebra

correlations,H = generatesystem(-0.3 + 0im,0.3 + 0im,-1 + 0im,-(5/(3*sqrt(3)))im,32)

crystalfock = correlations|>getoutspace

scale = Scale([2 0; 0 2], crystalfock|>getcrystal|>getspace)
@info("Performing rgblocking...",scale)
@info("Generating rgblocking transformation...")
blocker = @time scale * crystalfock
@info("Performing rgblocking on correlations...")
blockedcorrelations = @time blocker * correlations * blocker'
blockedcrystal = blockedcorrelations|>getoutspace|>getcrystal

Asublatticeunitcellmodes,Bsublatticeunitcellmodes = generaterefAandBsublatticeunitcellmodes(blockedcrystal)
Asublatticeoffsets = Subset((mode|>getattr(:b)) for mode in Asublatticeunitcellmodes)
Bsublatticeoffsets = Subset((mode|>getattr(:b)) for mode in Bsublatticeunitcellmodes)
visualize(Asublatticeoffsets)
visualize(Bsublatticeoffsets)

firstrgedcorrelations = firstgmerastep(blockedcorrelations,0.3,6,Asublatticeunitcellmodes,Bsublatticeunitcellmodes)

secondrgedcorrelations = secondgmerastep(firstrgedcorrelations,0.3,6,Asublatticeunitcellmodes,Bsublatticeunitcellmodes)

thirdrgedcorrelations = thirdgmerastep(secondrgedcorrelations,0.3,6,Asublatticeunitcellmodes,Bsublatticeunitcellmodes)
