using Zipper
using LinearAlgebra,Plots

setmaxthreads(setmaxthreads(Threads.nthreads()))
usecrystaldensemap()
plotlyjs()

power = 5
onsitepotential = 0
nnhopping = 0
systemsize=2^power
correlations,H = generatesystem(-onsitepotential+0im, onsitepotential+0im,-1 + 0im,-nnhopping,systemsize)
crystalfock = correlations|>getoutspace

scaling = 2

blockedcorrelations,blockedH,blocker = blocking(correlations,H,scaling)
# noofmodesinlocalreg = (6*(scaling^2))|>Int
# noofdistillablemodes = ((1/4)*noofmodesinlocalreg)|>Int
# noofcouriermodesinfirststep = noofmodesinlocalreg - noofdistillablemodes
# noofcouriermodes = noofcouriermodesinfirststep
# noofcouriermodesinsecondstep = noofcouriermodesinfirststep - noofdistillablemodes
# noofcouriermodesinthirdstep = noofcouriermodesinsecondstep - noofdistillablemodes
noofflavourpermode = 1

blockedcrystal = blockedcorrelations|>getoutspace|>getcrystal
blockedspace::RealSpace = blockedcrystal|>getspace
    
refrot = inv([2/3 -1/3; -1/3 2/3]')
c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
c3 = c6^2

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

identity = idmap(localcorrelations|>getoutspace, localcorrelations|>getinspace)
localfilledseeds = identity[:,firstinnerhexagonalregionhalf1fock]
localemptyseeds = identity[:,firstinnerhexagonalregionhalf2fock]
localcourierseeds = identity[:,firstboundaryhexagonalregionfock]
Ufilled,svdvals,Vdfilled = svd(localfilledseeds'*localfilledisometry)
Uempty,svdvals,Vdempty = svd(localemptyseeds'*localemptyisometry)

# rmlist = [m for m in Uempty|>getoutspace|>orderedmodes]
# filledrsmodes = Subset(rmlist[1:3])|>FockSpace
# emptyrsmodes = Subset(rmlist[4:6])|>FockSpace

localrotfilled = Ufilled*Vdfilled
localwannierfilled = localfilledisometry*localrotfilled'
localrotempty = Uempty*Vdempty
localwannierempty = localemptyisometry*localrotempty'

[(localwannierfilled'*localcorrelations*localwannierfilled)[m,m]|>rep for m in (localwannierfilled'*localcorrelations*localwannierfilled)|>getinspace|>orderedmodes]
(localwannierfilled'*localcorrelations*localwannierfilled)|>eigspec|>visualize
localwanniercourier,minsvd = localwannierization(localcourierisometry, localcourierseeds)
# localwannierfrozen,minsvd = localwannierization(localfrozenisometry, localfrozenseeds)
localwannieriso = localwanniercourier+localwannierfilled+localwannierempty

localwannierresults =  Dict(:localwannieriso => localwannieriso, :localcourierseeds => localcourierseeds)
wannierinfos =  Dict(firsthexagonalregionfock=>localwannierresults)

shiftedfirstcenter1 = [1,0] ∈ blockedspace
shiftedfirstcenter2 = [0,1] ∈ blockedspace
shiftedfirstcenter3 = [1,1] ∈ blockedspace
    
firstshiftedhexagonalcenters = [shiftedfirstcenter1, shiftedfirstcenter2, shiftedfirstcenter3]

for hexagonalcenter in firstshiftedhexagonalcenters
    # hexagonalcenter = shiftedfirstcenter1
    shiftedfirsthexagonalregion = firsthexagonalregion.+hexagonalcenter
    shiftedfirsthexagonalregionfock = quantize(shiftedfirsthexagonalregion,noofflavourpermode)
    c6recenter = recenter(c6,hexagonalcenter)
    c3recenter = recenter(c3,hexagonalcenter)
    
    shiftedfirstinnerhexagonalregion = rankedandgroupoffsets(shiftedfirsthexagonalregion,1)
    shiftedoffstlist = [offset for offset in shiftedfirstinnerhexagonalregion]
    shiftedfirstinnerhexagonalregionhalf1 = Subset([shiftedoffstlist[1],c3recenter*shiftedoffstlist[1],(c3recenter)^2*shiftedoffstlist[1]])
    shiftedfirstinnerhexagonalregionhalf2 = Subset([c6recenter*shiftedoffstlist[1],c6recenter*c3recenter*shiftedoffstlist[1],c6recenter*(c3recenter)^2*shiftedoffstlist[1]])
    shiftedfirstinnerhexagonalregionfock = quantize(shiftedfirstinnerhexagonalregion,1)
    shiftedfirstinnerhexagonalregionhalf1fock = quantize(shiftedfirstinnerhexagonalregionhalf1,1)
    shiftedfirstinnerhexagonalregionhalf2fock = quantize(shiftedfirstinnerhexagonalregionhalf2,1)

    shiftedfirstboundaryhexagonalregion = shiftedfirsthexagonalregion-shiftedfirstinnerhexagonalregion
    shiftedfirstboundaryhexagonalregionfock = quantize(shiftedfirstboundaryhexagonalregion,1)

    shiftedlocalcorrelations = regioncorrelations(blockedcorrelations,shiftedfirsthexagonalregionfock)
    shiftedlocalstates = getregionstates(localcorrelations=shiftedlocalcorrelations, grouping=[3, 18, 3])
    shiftedlocalcourierisometry = shiftedlocalstates[2]|>FockMap
    shiftedlocalfilledisometry = shiftedlocalstates[1]|>FockMap
    shiftedlocalemptyisometry = shiftedlocalstates[3]|>FockMap
    shiftedlocalfrozenisometry = (shiftedlocalstates[1]+shiftedlocalstates[3])|>FockMap
    
    shiftedidentity = idmap(shiftedlocalcorrelations|>getoutspace, shiftedlocalcorrelations|>getinspace)
    # shiftedlocalfrozenseeds = shiftedidentity[:,shiftedfirstinnerhexagonalregionfock]
    shiftedlocalfilledseeds = shiftedidentity[:,shiftedfirstinnerhexagonalregionhalf1fock]
    shiftedlocalemptyseeds = shiftedidentity[:,shiftedfirstinnerhexagonalregionhalf2fock]
    shiftedlocalcourierseeds = shiftedidentity[:,shiftedfirstboundaryhexagonalregionfock]

    shiftedUfilled,shiftedsvdvals,shiftedVdfilled = svd(shiftedlocalfilledseeds'*shiftedlocalfilledisometry)
    shiftedUempty,shiftedsvdvals,shiftedVdempty = svd(shiftedlocalemptyseeds'*shiftedlocalemptyisometry)

    # shiftedrmlist = [m for m in shiftedUempty|>getoutspace|>orderedmodes]
    # shiftedfilledrsmodes = Subset(shiftedrmlist[1:3])|>FockSpace
    # shiftedemptyrsmodes = Subset(shiftedrmlist[4:6])|>FockSpace

    shiftedlocalrotfilled = shiftedUfilled*shiftedVdfilled
    shiftedlocalwannierfilled = shiftedlocalfilledisometry*shiftedlocalrotfilled'
    shiftedlocalrotempty = shiftedUempty*shiftedVdempty
    shiftedlocalwannierempty = shiftedlocalemptyisometry*shiftedlocalrotempty'

    shiftedlocalwanniercourier,minsvd = localwannierization(shiftedlocalcourierisometry, shiftedlocalcourierseeds)
    # shiftedlocalwannierfrozen,minsvd = localwannierization(shiftedlocalfrozenisometry, shiftedlocalfrozenseeds)
    shiftedlocalwannieriso = shiftedlocalwanniercourier + shiftedlocalwannierfilled + shiftedlocalwannierempty

    shiftedlocalwannierresults =  Dict(:localwannieriso => shiftedlocalwannieriso, :localcourierseeds => shiftedlocalcourierseeds)
    wannierinfos[shiftedfirsthexagonalregionfock] =  shiftedlocalwannierresults
end

ref = [quantize(firsthexagonalregion.+hexagonalcenter,noofflavourpermode) for hexagonalcenter in firstshiftedhexagonalcenters]
firsthexagonalregionfocklist = [firsthexagonalregionfock,ref...]

extendedwannierizediso =  sum(wannierinfos[regionfock][:localwannieriso] for regionfock in firsthexagonalregionfocklist)
                
origin = [0, 0] ∈ blockedspace
refunictcellfock = FockSpace(Subset(mode for mode in extendedwannierizediso |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))

disentanglerisometry = globalwannierfunction(blockedcorrelations,extendedwannierizediso[:,refunictcellfock])
disentangledidentity = disentanglerisometry'*disentanglerisometry


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

identity = idmap(localcorrelations|>getoutspace, localcorrelations|>getinspace)
localfrozenseeds = identity[:,secondinnerhexagonalregionfock]
localcourierseeds = identity[:,secondboundaryhexagonalregionfock]

localwanniercourier,minsvd = localwannierization(localcourierisometry, localcourierseeds)
localwannierfrozen,minsvd = localwannierization(localfrozenisometry, localfrozenseeds)
localwannieriso = localwanniercourier+localwannierfrozen

localwannierresults =  Dict(:localwannieriso => localwannieriso, :localcourierseeds => localcourierseeds, :localfrozenseeds => localfrozenseeds)
wannierinfos =  Dict(secondhexagonalregionfock=>localwannierresults)

shiftedsecondcenter1 = [0,1] ∈ firstgmeraspace
shiftedsecondcenter2 = [-1,1] ∈ firstgmeraspace

secondshiftedhexagonalcenters = [shiftedsecondcenter1, shiftedsecondcenter2]

for hexagonalcenter in secondshiftedhexagonalcenterlist
    shiftedsecondhexagonalregion = secondhexagonalregion.+hexagonalcenter
    shiftedsecondhexagonalregionfock = quantize(shiftedsecondhexagonalregion,noofflavourpermode)
    c6recenter = recenter(c6,secondcenter+hexagonalcenter)
    c3recenter = recenter(c3,secondcenter+hexagonalcenter)

    shiftedsecondinnerhexagonalregion = rankedandgroupoffsets(shiftedsecondhexagonalregion,1)
    shiftedsecondinnerhexagonalregionfock = quantize(shiftedsecondinnerhexagonalregion,1)

    shiftedsecondboundaryhexagonalregion = shiftedsecondhexagonalregion-shiftedsecondinnerhexagonalregion
    shiftedsecondboundaryhexagonalregionfock = quantize(shiftedsecondboundaryhexagonalregion,1)

    shiftedlocalcorrelations = regioncorrelations(rotcorrelations,shiftedsecondhexagonalregionfock)
    shiftedlocalstates = getregionstates(localcorrelations=shiftedlocalcorrelations, grouping=[3, 18, 3])
    shiftedlocalcourierisometry = shiftedlocalstates[2]|>FockMap
    shiftedlocalfilledisometry = shiftedlocalstates[1]|>FockMap
    shiftedlocalemptyisometry = shiftedlocalstates[3]|>FockMap
    shiftedlocalfrozenisometry = (shiftedlocalstates[1]+shiftedlocalstates[3])|>FockMap
    
    shiftedidentity = idmap(shiftedlocalcorrelations|>getoutspace, shiftedlocalcorrelations|>getinspace)
    shiftedlocalfrozenseeds = shiftedidentity[:,shiftedsecondinnerhexagonalregionfock]
    shiftedlocalcourierseeds = shiftedidentity[:,shiftedsecondboundaryhexagonalregionfock]

    shiftedlocalwanniercourier,minsvd = localwannierization(shiftedlocalcourierisometry, shiftedlocalcourierseeds)
    shiftedlocalwannierfrozen,minsvd = localwannierization(shiftedlocalfrozenisometry, shiftedlocalfrozenseeds)
    shiftedlocalwannieriso = shiftedlocalwannierfrozen + shiftedlocalwanniercourier

    
    shiftedlocalwannierresults =  Dict(:localwannieriso => shiftedlocalwannieriso, :localcourierseeds => shiftedlocalcourierseeds, :localfrozenseeds => shiftedlocalfrozenseeds)
    wannierinfos[shiftedsecondhexagonalregionfock] =  shiftedlocalwannierresults
end

ref = [quantize(secondhexagonalregion.+hexagonalcenter,noofflavourpermode) for hexagonalcenter in secondshiftedhexagonalcenters]
secondhexagonalregionfocklist = [secondhexagonalregionfock,ref...]

extendedwannierizediso =  sum(wannierinfos[regionfock][:localwannieriso] for regionfock in secondhexagonalregionfocklist)
                
origin = [0, 0] ∈ blockedspace
refunictcellfock = FockSpace(Subset(mode for mode in extendedwannierizediso |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))

disentanglerisometry2 = globalwannierfunction(blockedcorrelations,extendedwannierizediso[:,refunictcellfock])


rotrotcorrelations = disentanglerisometry2' * rotcorrelations * disentanglerisometry2
rotrotH = disentanglerisometry2' * rotH * disentanglerisometry2

secondgmeracrystalfock = rotrotcorrelations|>getoutspace
secondgmeracrystal::Crystal = secondgmeracrystalfock|>getcrystal
secondgmeraspace::RealSpace = secondgmeracrystal|>getspace

refrot = inv([2/3 -1/3; -1/3 2/3]')
c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
c3 = c6^2

thirdcenter = [1/3,1/3] ∈ secondgmeraspace
thirdhexagonalregion = gethexagonalregion(rot = refrot,crystal=secondgmeracrystal, center=thirdcenter, metricspace=secondgmeraspace)
thirdhexagonalregionfock = quantize(thirdhexagonalregion,noofflavourpermode)
rgshiftedcenter1 = [2/3,-1/3] ∈ secondgmeraspace
thirdrgshiftedhexagonalregion1 = thirdhexagonalregion.+rgshiftedcenter1
rgshiftedcenter2 = [1/3,1/3] ∈ secondgmeraspace
thirdrgshiftedhexagonalregion2 = thirdhexagonalregion.+rgshiftedcenter2

c6recenter = recenter(c6,thirdcenter)
c3recenter = recenter(c3,thirdcenter)

localcorrelations = regioncorrelations(rotrotcorrelations,thirdhexagonalregionfock)
localspectrum = localcorrelations|>eigspech
display(localspectrum|>visualize)
