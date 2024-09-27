using Zipper
using LinearAlgebra,Plots

function offsetofmodes(modes::Subset)
    return Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in modes)
end


function eecontour(localcorrelations::SparseFockMap)
    localspectrum = localcorrelations|>eigspec
    emodewifevals = localspectrum|>geteigenvalues
    evectors = localspectrum|>geteigenvectors
    rsfock = localcorrelations|>getoutspace
    result::Dict{Mode, Number} = Dict()
    for rsmode in rsfock
        save = []
        for (emode,eval) in emodewifevals
            append!(save,norm(evectors[rsmode,emode])^2*ee(norm(eval)))
        end
        result[rsmode] = sum(save)
    end
    return result
end

setmaxthreads(setmaxthreads(Threads.nthreads()))
usecrystaldensemap()
plotlyjs()

power = 5
onsitepotential = 0
nnhopping = 0.2im
systemsize=2^power
correlations,H = generatesystem(-onsitepotential+0im, onsitepotential+0im,-1 + 0im,-nnhopping,systemsize)
crystalfock = correlations|>getoutspace
crystal = crystalfock|>getcrystal

scaling = 2

blockedcorrelations,blockedH,blocker = blocking(correlations,H,scaling)

blockedcrystalfock = blockedcorrelations|>getoutspace
blockedcrystal = blockedcorrelations|>getoutspace|>getcrystal
blockedspace::RealSpace = blockedcrystal|>getspace

refrot = inv([2/3 -1/3; -1/3 2/3]')
c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
c3 = c6^2

noofflavourpermode=1

firstcenter = [0,0] ∈ blockedspace
firsthexagonalregion = gethexagonalregion(rot = refrot,crystal=blockedcrystal, center=firstcenter, metricspace=blockedspace)
firsthexagonalregionfock = quantize(firsthexagonalregion,noofflavourpermode)
firsthexagonalregion|>visualize

# @info "Computing local correlations..."
localrestrict = fourier(blockedcrystalfock, firsthexagonalregionfock) / (blockedcrystal|>vol|>sqrt)
localcorrelations = localrestrict'*blockedcorrelations*localrestrict
localspectrum = localcorrelations|>eigspec

