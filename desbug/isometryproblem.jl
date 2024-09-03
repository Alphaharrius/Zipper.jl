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
sortgroupdictwifdist(firsteecontour,true,3)[4][2]|>offsetofmodes|>visualize
firsteecontour = eecontour(localcorrelations)
ee([eval for (m,eval) in  localspectrum|>geteigenvalues][1])

localcorrelations|>getoutspace|>orderedmodes
offsetofmodes(localcorrelations|>getoutspace|>orderedmodes)|>visualize

setmaxthreads(setmaxthreads(Threads.nthreads()))
usecrystaldensemap()
plotlyjs()

power = 5
onsitepotential = 0
nnhopping = 0
systemsize=2^power
bonds,correlations,H = generatesystem(-onsitepotential+0im, onsitepotential+0im,-1 + 0im,-nnhopping,systemsize)
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

firstinnerhexagonalregion = rankedandgroupoffsets(firsthexagonalregion,1)
offstlist = [offset for offset in firstinnerhexagonalregion]
firstinnerhexagonalregionhalf1 = Subset([offstlist[1],c3*offstlist[1],(c3)^2*offstlist[1]])
firstinnerhexagonalregionhalf2 = Subset([c6*offstlist[1],c6*c3*offstlist[1],c6*(c3)^2*offstlist[1]])
firstinnerhexagonalregionfock = quantize(firstinnerhexagonalregion,1)
firstinnerhexagonalregionhalf1fock = quantize(firstinnerhexagonalregionhalf1,1)
firstinnerhexagonalregionhalf2fock = quantize(firstinnerhexagonalregionhalf2,1)

firstboundaryhexagonalregion = firsthexagonalregion-firstinnerhexagonalregion
firstboundaryhexagonalregionfock = quantize(firstboundaryhexagonalregion,1)

# firstinnerhexagonalregion = rankedandgroupoffsets(firsthexagonalregion,2)
# offstlist = [offset for offset in firstinnerhexagonalregion]
# firstinnerhexagonalregionhalf1 = Subset([offstlist[end],c3*offstlist[end],(c3)^2*offstlist[end]])
# firstinnerhexagonalregionhalf2 = Subset([c6*offstlist[end],c6*c3*offstlist[end],c6*(c3)^2*offstlist[end]])
# firstinnerhexagonalregionfock = quantize(firstinnerhexagonalregion,1)
# firstinnerhexagonalregionhalf1fock = quantize(firstinnerhexagonalregionhalf1,1)
# firstinnerhexagonalregionhalf2fock = quantize(firstinnerhexagonalregionhalf2,1)

# firstboundaryhexagonalregion = firsthexagonalregion-firstinnerhexagonalregionhalf1-firstinnerhexagonalregionhalf2
# firstboundaryhexagonalregionfock = quantize(firstboundaryhexagonalregion,1)

# @info "Computing local correlations..."
localrestrict = fourier(blockedcrystalfock, firsthexagonalregionfock) / (blockedcrystal|>vol|>sqrt)
localcorrelations = localrestrict'*blockedcorrelations*localrestrict
localspectrum = localcorrelations|>eigspec
display(localspectrum|>visualize)

localstates = getregionstates(localcorrelations=localcorrelations, grouping=[3,18,3])
localcourierisometry = localstates[2]|>FockMap
localfilledisometry = localstates[1]|>FockMap
localemptyisometry = localstates[3]|>FockMap

iden = idmap(localcorrelations|>getoutspace, localcorrelations|>getinspace)
localfilledseeds = iden[:,firstinnerhexagonalregionhalf1fock]
localemptyseeds = iden[:,firstinnerhexagonalregionhalf2fock]
localcourierseeds = iden[:,firstboundaryhexagonalregionfock]
Ufilled,svdvals,Vdfilled = svd(localfilledseeds'*localfilledisometry)
Uempty,svdvals,Vdempty = svd(localemptyseeds'*localemptyisometry)

localrotfilled = Ufilled*Vdfilled
localwannierfilled = localfilledisometry*localrotfilled'
localrotempty = Uempty*Vdempty
localwannierempty = localemptyisometry*localrotempty'
localwanniercourier,minsvd = localwannierization(localcourierisometry, localcourierseeds)
localdisentangler = localwanniercourier+localwannierfilled+localwannierempty
redefinelocaldisentangler = FockMap(localdisentangler|>getoutspace,localdisentangler|>getinspace|>RegionFock,localdisentangler|>rep)

newblockedcrystalfock = getcrystalfock(redefinelocaldisentangler|>getinspace|>unitcellfock, blockedcrystal)
newlocalrestrict = fourier(newblockedcrystalfock, redefinelocaldisentangler|>getinspace) * (blockedcrystal|>vol|>sqrt)

globaldisentangler = broadcast(*,(localrestrict*redefinelocaldisentangler), newlocalrestrict')
rotcorrelations = globaldisentangler'*blockedcorrelations*globaldisentangler
rotH = globaldisentangler'*blockedH*globaldisentangler

firstgmeracrystalfock = rotcorrelations|>getoutspace
firstgmeracrystal::Crystal = firstgmeracrystalfock|>getcrystal
firstgmeraspace::RealSpace = firstgmeracrystal|>getspace

refrot = inv([2/3 -1/3; -1/3 2/3]')
c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
c3 = c6^2

secondcenter = [2/3,-1/3] ∈ firstgmeraspace
secondhexagonalregion = gethexagonalregion(rot=refrot,crystal=firstgmeracrystal, center=secondcenter, metricspace=firstgmeraspace)
secondhexagonalregionfock = quantize(secondhexagonalregion,noofflavourpermode)

secondinnerhexagonalregion = rankedandgroupoffsets(secondhexagonalregion,1)
secondoffstlist = [offset for offset in secondinnerhexagonalregion]
firstinnerhexagonalregionhalf1 = Subset([secondoffstlist[1],c3*secondoffstlist[1],(c3)^2*secondoffstlist[1]])
firstinnerhexagonalregionhalf2 = Subset([c6*offstlist[1],c6*c3*offstlist[1],c6*(c3)^2*offstlist[1]])
firstinnerhexagonalregionfock = quantize(firstinnerhexagonalregion,1)
firstinnerhexagonalregionhalf1fock = quantize(firstinnerhexagonalregionhalf1,1)
firstinnerhexagonalregionhalf2fock = quantize(firstinnerhexagonalregionhalf2,1)


secondboundaryhexagonalregion = secondhexagonalregion-secondinnerhexagonalregion
secondboundaryhexagonalregionfock = quantize(secondboundaryhexagonalregion,1)

c6recenter = recenter(c6,secondcenter)
c3recenter = recenter(c3,secondcenter)

localcorrelations = regioncorrelations(rotcorrelations,secondhexagonalregionfock)
localspectrum = localcorrelations|>eigspech
display(localspectrum|>visualize)
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[6, 12, 6])
localcourierisometry = localstates[2]|>FockMap
localfilledisometry = localstates[1]|>FockMap
localemptyisometry = localstates[3]|>FockMap