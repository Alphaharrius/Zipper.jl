using Zipper
using LinearAlgebra,Plots

function computeenergyspectrum(bonds::FockMap; crystal::Crystal)::CrystalSpectrum
    # function computeenergyspectrum(bonds::FockMap; crystal::Crystal)
        bondmodes::Subset{Mode} = bonds|>getoutspace|>orderedmodes
        homefock::NormalFock = bonds|>getoutspace|>RegionFock|>unitcellfock
        transform = fourier(getcrystalfock(homefock, crystal), bondmodes|>RegionFock)
        function compute(k)
            ktransform = transform[k, :]
            return k=>(ktransform*bonds*ktransform')
        end
        momentumhamiltonians = paralleltasks(
            name="computeenergyspectrum",
            tasks=(()->compute(k) for k in crystal|>brillouinzone),
            count=crystal|>vol)|>parallel
        return crystalspectrum(momentumhamiltonians, crystal=crystal)
    end


blockedhomefock = blockedcrystalfock|>unitcellfock
refblockedcrystalfock = getcrystalfock(blockedhomefock, blockedcrystal)
newblockedhomefock = localwannieriso|>getoutspace
newunitcell = Subset((m|>getattr(:b))+(m|>getattr(:r)) for m in newblockedhomefock|>orderedmodes)
newblockedcrystal = Crystal(newunitcell, blockedcrystal|>size)
refnewblockedcrystalfock = getcrystalfock(newblockedhomefock, newblockedcrystal)
refnewblockedcrystalfock|>getcrystal|>getunitcell|>visualize

fourier(refblockedcrystalfock, newblockedhomefock)
fourier(refnewblockedcrystalfock, newblockedhomefock)
[m for m in refnewblockedcrystalfock|>unitcellfock]
[m for m in newblockedhomefock|>orderedmodes]

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
rgshiftedcenter1 = [2/3,-1/3] ∈ blockedspace
firstrgshiftedhexagonalregion1 = firsthexagonalregion.+rgshiftedcenter1
rgshiftedcenter2 = [1/3,1/3] ∈ blockedspace
firstrgshiftedhexagonalregion2 = firsthexagonalregion.+rgshiftedcenter2

firstinnerhexagonalregion = rankedandgroupoffsets(firsthexagonalregion,1)
offstlist = [offset for offset in firstinnerhexagonalregion]
firstinnerhexagonalregionhalf1 = Subset([offstlist[1],c3*offstlist[1],(c3)^2*offstlist[1]])
firstinnerhexagonalregionhalf2 = Subset([c6*offstlist[1],c6*c3*offstlist[1],c6*(c3)^2*offstlist[1]])
firstinnerhexagonalregionfock = quantize(firstinnerhexagonalregion,1)
firstinnerhexagonalregionhalf1fock = quantize(firstinnerhexagonalregionhalf1,1)
firstinnerhexagonalregionhalf2fock = quantize(firstinnerhexagonalregionhalf2,1)

firstboundaryhexagonalregion = firsthexagonalregion-firstinnerhexagonalregion
firstboundaryhexagonalregionfock = quantize(firstboundaryhexagonalregion,1)

localcorrelations = regioncorrelations(blockedcorrelations,firsthexagonalregionfock)
localspectrum = localcorrelations|>eigspec
display(localspectrum|>visualize)

localstates = getregionstates(localcorrelations=localcorrelations, grouping=[3,18,3])
localcourierisometry = localstates[2]|>FockMap
localfilledisometry = localstates[1]|>FockMap
localemptyisometry = localstates[3]|>FockMap
localfrozenisometry = (localstates[1]+localstates[3])|>FockMap

@info "Computing local correlations..."
localrestrict = fourier(blockedcrystalfock, firsthexagonalregionfock) / (blockedcrystal|>vol|>sqrt)
localcorrelations = localrestrict' * blockedcorrelations * localrestrict
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[3, 18, 3])


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
localwannieriso = localwanniercourier+localwannierfilled+localwannierempty

localwannierresults =  Dict(:localwanniercourier => localwanniercourier,:localwannierfilled => localwannierfilled,:localwannierempty => localwannierempty, :localwannieriso => localwannieriso, :localcourierseeds => localcourierseeds)
wannierinfos =  Dict(firsthexagonalregionfock=>localwannierresults)
localwannierresults[:localwannieriso]
crystalisometry(localisometry=localwannierresults[:localwannieriso],crystalfock=(blockedcorrelations|>getinspace))

disentanglerisometry = globalwannierfunction(blockedcorrelations,localwannierresults[:localwannieriso])
rotcorrelations = disentanglerisometry' * blockedcorrelations * disentanglerisometry
rotH = disentanglerisometry' * blockedH * disentanglerisometry

firstgmeracrystalfock = rotcorrelations|>getoutspace
firstgmeracrystal::Crystal = firstgmeracrystalfock|>getcrystal
firstgmeraspace::RealSpace = firstgmeracrystal|>getspace

refrot = inv([2/3 -1/3; -1/3 2/3]')
c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
c3 = c6^2

secondcenter = [2/3,-1/3] ∈ firstgmeraspace
secondhexagonalregion = gethexagonalregion(rot=refrot,crystal=firstgmeracrystal, center=secondcenter, metricspace=firstgmeraspace)
secondhexagonalregionfock = quantize(secondhexagonalregion,noofflavourpermode)
rgshiftedcenter1 = [2/3,-1/3] ∈ firstgmeraspace
secondrgshiftedhexagonalregion1 = secondhexagonalregion.+rgshiftedcenter1
rgshiftedcenter2 = [1/3,1/3] ∈ firstgmeraspace
secondrgshiftedhexagonalregion2 = secondhexagonalregion.+rgshiftedcenter2

secondinnerhexagonalregion = rankedandgroupoffsets(secondhexagonalregion,1)
secondinnerhexagonalregionfock = quantize(secondinnerhexagonalregion,1)

secondboundaryhexagonalregion = secondhexagonalregion-secondinnerhexagonalregion
secondboundaryhexagonalregionfock = quantize(secondboundaryhexagonalregion,1)

c6recenter = recenter(c6,secondcenter)
c3recenter = recenter(c3,secondcenter)

localcorrelations = regioncorrelations(rotcorrelations,secondhexagonalregionfock)
localspectrum = localcorrelations|>eigspech
display(localspectrum|>visualize)
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[3, 18, 3])
localcourierisometry = localstates[2]|>FockMap
localfilledisometry = localstates[1]|>FockMap
localemptyisometry = localstates[3]|>FockMap
localfrozenisometry = (localstates[1]+localstates[3])|>FockMap