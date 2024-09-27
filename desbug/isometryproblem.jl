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
crystal = crystalfock|>getcrystal

scaling = 4

blockedcorrelations,blockedH,blocker = blocking(correlations,H,scaling)

blockedcrystalfock = blockedcorrelations|>getoutspace
blockedcrystal = blockedcorrelations|>getoutspace|>getcrystal
blockedspace::RealSpace = blockedcrystal|>getspace
    
refrot = inv([2/3 -1/3; -1/3 2/3]')
c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
c3 = c6^2

noofflavourpermode=1

firstcenter = [1/3,1/3] ∈ blockedspace
firsthexagonalregion = gethexagonalregion(rot = refrot,crystal=blockedcrystal, center=firstcenter, metricspace=blockedspace)
firsthexagonalregionfock = quantize(firsthexagonalregion,noofflavourpermode)
firsthexagonalmodes = Subset(mode for mode in firsthexagonalregionfock)

# @info "Computing local correlations..."
localrestrict = fourier(blockedcrystalfock, firsthexagonalregionfock) / (blockedcrystal|>vol|>sqrt)
localcorrelations = localrestrict'*blockedcorrelations*localrestrict
potentialfilledseeds = Subset([mode for mode in firsthexagonalregionfock if (localcorrelations[mode,mode]|>rep)[1,1]|>abs<0.5])
potentialemptyseeds = Subset([mode for mode in firsthexagonalregionfock if (localcorrelations[mode,mode]|>rep)[1,1]|>abs>0.5])

localspectrum = localcorrelations|>eigspec
localspectrum|>visualize
# emodewifevals = localspectrum|>geteigenvalues
# # chosenfilledevecmodes = append!(sort([(m,abs(val)) for (m,val) in emodewifevals],by=x->x[2])[1:12],sort([(m,abs(val)) for (m,val) in emodewifevals],by=x->x[2])[end-11:end])
# chosenfilledevecmodes = sort([(m,abs(val)) for (m,val) in emodewifevals],by=x->x[2])[37:end-36]
# # chosenfilledevecmodes = sort([(m,abs(val)) for (m,val) in emodewifevals],by=x->x[2])[13:end-12]
# evectors = localspectrum|>geteigenvectors
# rsfock = localcorrelations|>getoutspace
# result::Dict{Mode, Number} = Dict()
# for rsmode in rsfock
#     save = []
#     for (emode,eval) in chosenfilledevecmodes
#         append!(save,norm(evectors[rsmode,emode])^2*ee(norm(eval)))
#     end
#     result[rsmode] = sum(save)
# end
# sortgroupdictwifvalue(result,true)
# ((sortgroupdictwifvalue(result,true)[1][2]|>offsetofmodes)+(sortgroupdictwifvalue(result,true)[2][2]|>offsetofmodes))|>visualize
# intersect(sortgroupdictwifvalue(result,true)[5][2],potentialfilledseeds)
# sortgroupdictwifvalue(result,true)[end][2]|>offsetofmodes|>visualize

entanglementcontourinfo = entanglementcontour(localcorrelations)
groupedandsortedmodeswifecontour = sortgroupdictwifvalue(entanglementcontourinfo,false)
# nooffrozenmodes = div(firsthexagonalregionfock|>length,4)*1
nooffrozenmodes = 96
frozenseeds = groupedandsortedmodeswifecontour[1][2]
for (entanglementval,modes) in groupedandsortedmodeswifecontour[2:end]
    if (frozenseeds|>length)<nooffrozenmodes
        @info "include more seeds"
        frozenseeds = frozenseeds+modes
    elseif (frozenseeds|>length)==nooffrozenmodes
        @info "just enough frozenseeds"
        break
    else 
        @error "too many frozenseeds, sth wrong"
    end
end
frozenseeds|>offsetofmodes|>visualize
potentialfilledseeds|>offsetofmodes|>visualize
filledseeds = intersect(frozenseeds,potentialfilledseeds)
emptyseeds = intersect(frozenseeds,potentialemptyseeds)
courierseeds = firsthexagonalmodes - frozenseeds
# courierseeds|>offsetofmodes|>visualize
filledseedsfock = quantize(filledseeds|>offsetofmodes,noofflavourpermode)
emptyseedsfock = quantize(emptyseeds|>offsetofmodes,noofflavourpermode)
# courierseedsfock = quantize(courierseeds|>offsetofmodes,noofflavourpermode)


# firstinnerhexagonalregion = rankedandgroupoffsets(firsthexagonalregion,1)
# offstlist = [offset for offset in firstinnerhexagonalregion]
# firstinnerhexagonalregionhalf1 = Subset([offstlist[1],c3*offstlist[1],(c3)^2*offstlist[1]])
# firstinnerhexagonalregionhalf2 = Subset([c6*offstlist[1],c6*c3*offstlist[1],c6*(c3)^2*offstlist[1]])
# firstinnerhexagonalregionfock = quantize(firstinnerhexagonalregion,1)
# firstinnerhexagonalregionhalf1fock = quantize(firstinnerhexagonalregionhalf1,1)
# firstinnerhexagonalregionhalf2fock = quantize(firstinnerhexagonalregionhalf2,1)

# firstboundaryhexagonalregion = firsthexagonalregion-firstinnerhexagonalregion
# firstboundaryhexagonalregionfock = quantize(firstboundaryhexagonalregion,1)

localstates = getregionstates(localcorrelations=localcorrelations, grouping=[48,48 ])
# localcourierisometry = localstates[2]|>FockMap
localfilledisometry = localstates[1]|>FockMap
localemptyisometry = localstates[2]|>FockMap

iden = idmap(localcorrelations|>getoutspace, localcorrelations|>getinspace)
localfilledseedsisometry = iden[:,filledseedsfock]
localemptyseedsisometry = iden[:,emptyseedsfock]
# localcourierseedsisometry = iden[:,courierseedsfock]

localwannierfilled,minsvdempty = localwannierization(localfilledisometry, localfilledseedsisometry)
localwannierempty,minsvdfilled = localwannierization(localemptyisometry, localemptyseedsisometry)
# localwanniercourier,minsvdcourier = localwannierization(localcourierisometry, localcourierseedsisometry)
localdisentangler = localwannierfilled+localwannierempty
# localdisentangler = localwannierfilled+localwannierempty
redefinelocaldisentangler = FockMap(localdisentangler|>getoutspace,localdisentangler|>getinspace|>RegionFock,localdisentangler|>rep)

# inspacemodes = Subset(m for m in redefinelocaldisentangler|>getinspace)
# inspaceoffsets = inspacemodes|>offsetofmodes
# rotlocalcorrelations = redefinelocaldisentangler*localcorrelations*redefinelocaldisentangler'
# chosenfilledmodes = Subset([m for m in inspacemodes if (rotlocalcorrelations[m,m]|>rep|>diag)[1,1]|>abs<0.2])
# chosenfilledoffsets = chosenfilledmodes|>offsetofmodes
# chosenemptymodes = Subset([m for m in inspacemodes if (rotlocalcorrelations[m,m]|>rep|>diag)[1,1]|>abs>1-0.2])
# chosenemptyoffsets = chosenemptymodes|>offsetofmodes
# (inspaceoffsets-chosenfilledoffsets-chosenemptyoffsets)|>visualize
# (chosenfilledoffsets+chosenemptyoffsets)|>visualize
# chosenregionfock = (chosenfilledmodes+chosenemptymodes)|>RegionFock
# restredefinelocaldisentangler = redefinelocaldisentangler[chosenregionfock,:]
# restredefinelocaldisentangler*localcorrelations*restredefinelocaldisentangler'
# # [ for val in (redefinelocaldisentangler*localcorrelations*redefinelocaldisentangler')|>rep|>diag if abs(val)<0.03] 
# scatter(sort([abs(val) for val in (redefinelocaldisentangler*localcorrelations*redefinelocaldisentangler')|>rep|>diag]))
# scatter(sort([abs(val) for val in (restredefinelocaldisentangler*localcorrelations*restredefinelocaldisentangler')|>rep|>diag]))

newblockedcrystalfock = getcrystalfock(redefinelocaldisentangler|>getinspace|>unitcellfock, blockedcrystal)
# newcouriercrystalfock = getcrystalfock(localwanniercourier|>getinspace|>unitcellfock, Crystal(courierseeds|>offsetofmodes,blockedcrystal.sizes))
# newfilledcrystalfock = getcrystalfock(localwannierfilled|>getinspace|>unitcellfock, Crystal(filledseeds|>offsetofmodes,blockedcrystal.sizes))
# newemptycrystalfock = getcrystalfock(localwannierempty|>getinspace|>unitcellfock, Crystal(emptyseeds|>offsetofmodes,blockedcrystal.sizes))

newlocalrestrict = fourier(newblockedcrystalfock, redefinelocaldisentangler|>getinspace) * (blockedcrystal|>vol|>sqrt)
# newlocalcourierrestrict = fourier(newcouriercrystalfock, localwanniercourier|>getinspace) * (blockedcrystal|>vol|>sqrt)
# newlocalfilledrestrict = fourier(newfilledcrystalfock, localwannierfilled|>getinspace) * (blockedcrystal|>vol|>sqrt)
# newlocalemptyrestrict = fourier(newemptycrystalfock, localwannierempty|>getinspace) * (blockedcrystal|>vol|>sqrt)

globaldisentangler = broadcast(*,(localrestrict*redefinelocaldisentangler), newlocalrestrict')
# globalcourierisometry = broadcast(*,(localrestrict*localwanniercourier), newlocalcourierrestrict')
# globalfilledisometry = broadcast(*,(localrestrict*localwannierfilled), newlocalfilledrestrict')
# globalemptyisometry = broadcast(*,(localrestrict*localwannierempty), newlocalemptyrestrict')

rotcorrelations = globaldisentangler'*blockedcorrelations*globaldisentangler
rotH = globaldisentangler'*blockedH*globaldisentangler
# couriercorrelations = globalcourierisometry'*blockedcorrelations*globalcourierisometry
# couriercorrelationspectrum = couriercorrelations|>crystalspectrum
# # couriercorrelationspectrum|>visualize
# purifiedcouriercorrelations = roundingpurification(couriercorrelationspectrum)|>crystalfockmap

firstgmeracrystalfock = rotcorrelations|>getoutspace
firstgmeracrystal::Crystal = firstgmeracrystalfock|>getcrystal
firstgmeraspace::RealSpace = firstgmeracrystal|>getspace

refrot = inv([2/3 -1/3; -1/3 2/3]')
c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
c3 = c6^2

secondcenter = [2/3,-1/3] ∈ firstgmeraspace
secondhexagonalregion = gethexagonalregion(rot=refrot,crystal=firstgmeracrystal, center=secondcenter, metricspace=firstgmeraspace)
secondhexagonalregionfock = quantize(secondhexagonalregion,noofflavourpermode)
secondhexagonalmodes = Subset(mode for mode in secondhexagonalregionfock)

# secondinnerhexagonalregion = rankedandgroupoffsets(secondhexagonalregion,1)
# secondoffstlist = [offset for offset in secondinnerhexagonalregion]
# firstinnerhexagonalregionhalf1 = Subset([secondoffstlist[1],c3*secondoffstlist[1],(c3)^2*secondoffstlist[1]])
# firstinnerhexagonalregionhalf2 = Subset([c6*offstlist[1],c6*c3*offstlist[1],c6*(c3)^2*offstlist[1]])
# firstinnerhexagonalregionfock = quantize(firstinnerhexagonalregion,1)
# firstinnerhexagonalregionhalf1fock = quantize(firstinnerhexagonalregionhalf1,1)
# firstinnerhexagonalregionhalf2fock = quantize(firstinnerhexagonalregionhalf2,1)


# secondboundaryhexagonalregion = secondhexagonalregion-secondinnerhexagonalregion
# secondboundaryhexagonalregionfock = quantize(secondboundaryhexagonalregion,1)

c6recenter = recenter(c6,secondcenter)
c3recenter = recenter(c3,secondcenter)

# @info "Computing local correlations..."
localrestrict = fourier(firstgmeracrystalfock, secondhexagonalregionfock) / (firstgmeracrystal|>vol|>sqrt)
localcorrelations = localrestrict'*rotcorrelations*localrestrict
potentialfilledseeds2 = Subset([mode for mode in secondhexagonalregionfock if (localcorrelations[mode,mode]|>rep)[1,1]|>abs<0.5])
potentialemptyseeds2 = Subset([mode for mode in secondhexagonalregionfock if (localcorrelations[mode,mode]|>rep)[1,1]|>abs>0.5])

localspectrum = localcorrelations|>eigspech
localspectrum|>visualize
entanglementcontourinfo = entanglementcontour(localcorrelations)
groupedandsortedmodeswifecontour = sortgroupdictwifvalue(entanglementcontourinfo,false)
# nooffrozenmodes = div(firsthexagonalregionfock|>length,4)*1
nooffrozenmodes = 96
frozenseeds = groupedandsortedmodeswifecontour[1][2]
for (entanglementval,modes) in groupedandsortedmodeswifecontour[2:end]
    if (frozenseeds|>length)<nooffrozenmodes
        @info "include more seeds"
        frozenseeds = frozenseeds+modes
    elseif (frozenseeds|>length)==nooffrozenmodes
        @info "just enough frozenseeds"
        break
    else 
        @error "too many frozenseeds, sth wrong"
    end
end
filledseeds = intersect(frozenseeds,potentialfilledseeds2)
emptyseeds = intersect(frozenseeds,potentialemptyseeds2)
courierseeds = secondhexagonalmodes - frozenseeds
# courierseeds|>offsetofmodes|>visualize
frozenseeds|>offsetofmodes|>visualize
filledseedsfock = quantize(filledseeds|>offsetofmodes,noofflavourpermode)
emptyseedsfock = quantize(emptyseeds|>offsetofmodes,noofflavourpermode)
# courierseedsfock = quantize(courierseeds|>offsetofmodes,noofflavourpermode)

localstates = getregionstates(localcorrelations=localcorrelations, grouping=[48,48])
# localcourierisometry = localstates[2]|>FockMap
localfilledisometry = localstates[1]|>FockMap
localemptyisometry = localstates[2]|>FockMap

iden = idmap(localcorrelations|>getoutspace, localcorrelations|>getinspace)
localfilledseedsisometry = iden[:,filledseedsfock]
localemptyseedsisometry = iden[:,emptyseedsfock]
# localcourierseedsisometry = iden[:,courierseedsfock]

localwannierfilled,minsvdempty = localwannierization(localfilledisometry, localfilledseedsisometry)
localwannierempty,minsvdfilled = localwannierization(localemptyisometry, localemptyseedsisometry)
# localwanniercourier,minsvdcourier = localwannierization(localcourierisometry, localcourierseedsisometry)
localdisentangler = localwannierfilled+localwannierempty
redefinelocaldisentangler = FockMap(localdisentangler|>getoutspace,localdisentangler|>getinspace|>RegionFock,localdisentangler|>rep)

newblockedcrystalfock = getcrystalfock(redefinelocaldisentangler|>getinspace|>unitcellfock, firstgmeracrystal)
# newcouriercrystalfock = getcrystalfock(localwanniercourier|>getinspace|>unitcellfock, Crystal(courierseeds|>offsetofmodes,blockedcrystal.sizes))

newlocalrestrict = fourier(newblockedcrystalfock, redefinelocaldisentangler|>getinspace) * (firstgmeracrystal|>vol|>sqrt)
# newlocalcourierrestrict = fourier(newcouriercrystalfock, localwanniercourier|>getinspace) * (firstgmeracrystal|>vol|>sqrt)

globaldisentangler = broadcast(*,(localrestrict*redefinelocaldisentangler), newlocalrestrict')
# globalcourierisometry = broadcast(*,(localrestrict*localwanniercourier), newlocalcourierrestrict')

rotcorrelations2 = globaldisentangler'*rotcorrelations*globaldisentangler
# couriercorrelations2 = globalcourierisometry'*purifiedcouriercorrelations*globalcourierisometry

# couriercorrelationspectrum2 = couriercorrelations2|>crystalspectrum
# couriercorrelationspectrum|>visualize
# purifiedcouriercorrelations2 = roundingpurification(couriercorrelationspectrum2)|>crystalfockmap

rotH2 = globaldisentangler'*rotH*globaldisentangler

secondgmeracrystalfock = rotcorrelations2 |>getoutspace
secondgmeracrystal::Crystal = secondgmeracrystalfock|>getcrystal
secondgmeraspace::RealSpace = secondgmeracrystal|>getspace

refrot = inv([2/3 -1/3; -1/3 2/3]')
c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
c3 = c6^2

thirdcenter = [0,0] ∈ secondgmeraspace
thirdhexagonalregion = gethexagonalregion(rot = refrot,crystal=secondgmeracrystal, center=thirdcenter, metricspace=secondgmeraspace)
thirdhexagonalregionfock = quantize(thirdhexagonalregion,noofflavourpermode)
thirdhexagonalmodes = Subset(mode for mode in thirdhexagonalregionfock)

c6recenter = recenter(c6,thirdcenter)
c3recenter = recenter(c3,thirdcenter)

# @info "Computing local correlations..."
localrestrict = fourier(secondgmeracrystalfock, thirdhexagonalregionfock) / (secondgmeracrystal|>vol|>sqrt)
localcorrelations = localrestrict'*rotcorrelations2 *localrestrict
potentialfilledseeds3 = Subset([mode for mode in thirdhexagonalregionfock if (localcorrelations[mode,mode]|>rep)[1,1]|>abs<0.5])
potentialemptyseeds3 = Subset([mode for mode in thirdhexagonalregionfock if (localcorrelations[mode,mode]|>rep)[1,1]|>abs>0.5])

sort([abs(val) for val in localcorrelations|>rep|>diag])|>scatter

inspacemodes = Subset(m for m in localcorrelations|>getinspace)
inspaceoffsets = inspacemodes|>offsetofmodes
chosenfilledmodes = Subset([m for m in inspacemodes if (localcorrelations[m,m]|>rep|>diag)[1,1]|>abs<0.03])
chosenfilledoffsets = chosenfilledmodes|>offsetofmodes
chosenemptymodes = Subset([m for m in inspacemodes if (localcorrelations[m,m]|>rep|>diag)[1,1]|>abs>1-0.03])
chosenemptyoffsets = chosenemptymodes|>offsetofmodes
(inspaceoffsets-chosenfilledoffsets-chosenemptyoffsets)|>visualize
(chosenfilledoffsets+chosenemptyoffsets)|>visualize


localspectrum = localcorrelations|>eigspech
localspectrum|>visualize

emodewifevals = localspectrum|>geteigenvalues
emodewifevalspair = sort([(m,eval) for (m,eval) in emodewifevals],by=x->x[2])
rankedemodewifevals = Subset((m,eval) for (m,eval) in emodewifevalspair[1:end])
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
trial = sortgroupdictwifvalue(result,false)[end][2]+sortgroupdictwifvalue(result,false)[end-1][2]+sortgroupdictwifvalue(result,false)[end-2][2]+sortgroupdictwifvalue(result,false)[end-3][2]
sortgroupdictwifvalue(result,false)[end-14][2]|>offsetofmodes|>visualize

entanglementcontourinfo = entanglementcontour(localcorrelations)
groupedandsortedmodeswifecontour = sortgroupdictwifvalue(entanglementcontourinfo,false)
# nooffrozenmodes = div(firsthexagonalregionfock|>length,4)*3
groupedandsortedmodeswifecontour
nooffrozenmodes = 18
frozenseeds = groupedandsortedmodeswifecontour[1][2]
for (entanglementval,modes) in groupedandsortedmodeswifecontour[2:end]
    if (frozenseeds|>length)<nooffrozenmodes
        @info "include more seeds"
        frozenseeds = frozenseeds+modes
    elseif (frozenseeds|>length)==nooffrozenmodes
        @info "just enough frozenseeds"
        break
    else 
        @error "too many frozenseeds, sth wrong"
    end
end
filledseeds = intersect(frozenseeds,potentialfilledseeds3)
emptyseeds = intersect(frozenseeds,potentialemptyseeds3)
courierseeds = sortgroupdictwifvalue(result,false)[end-14][2]
frozenseeds|>offsetofmodes|>visualize
courierseeds|>offsetofmodes|>visualize


filledseedsfock = quantize(filledseeds|>offsetofmodes,noofflavourpermode)
emptyseedsfock = quantize(emptyseeds|>offsetofmodes,noofflavourpermode)
courierseedsfock = quantize(courierseeds|>offsetofmodes,noofflavourpermode)

localstates = getregionstates(localcorrelations=localcorrelations, grouping=[45,6,45])
localcourierisometry = localstates[2]|>FockMap
localfilledisometry = localstates[1]|>FockMap
localemptyisometry = localstates[2]|>FockMap

iden = idmap(localcorrelations|>getoutspace, localcorrelations|>getinspace)
localfilledseedsisometry = iden[:,filledseedsfock]
localemptyseedsisometry = iden[:,emptyseedsfock]
localcourierseedsisometry = iden[:,courierseedsfock]

localwannierfilled,minsvdempty = localwannierization(localfilledisometry, localfilledseedsisometry)
localwannierempty,minsvdfilled = localwannierization(localemptyisometry, localemptyseedsisometry)
localwanniercourier,minsvdcourier = localwannierization(localcourierisometry, localcourierseedsisometry)
localdisentangler = localwannierfilled+localwannierempty
redefinelocaldisentangler = FockMap(localdisentangler|>getoutspace,localdisentangler|>getinspace|>RegionFock,localdisentangler|>rep)

rotlocalcorrelations = redefinelocaldisentangler*localcorrelations*redefinelocaldisentangler'
sort([abs(val) for val in rotlocalcorrelations|>rep|>diag])|>scatter

inspacemodes = Subset(m for m in rotlocalcorrelations|>getinspace)
inspaceoffsets = inspacemodes|>offsetofmodes
chosenfilledmodes = Subset([m for m in inspacemodes if (rotlocalcorrelations[m,m]|>rep|>diag)[1,1]|>abs<0.05])
chosenfilledoffsets = chosenfilledmodes|>offsetofmodes
chosenemptymodes = Subset([m for m in inspacemodes if (rotlocalcorrelations[m,m]|>rep|>diag)[1,1]|>abs>1-0.05])
chosenemptyoffsets = chosenemptymodes|>offsetofmodes
(inspaceoffsets-chosenfilledoffsets-chosenemptyoffsets)|>visualize
(chosenfilledoffsets+chosenemptyoffsets)|>visualize

newblockedcrystalfock = getcrystalfock(redefinelocaldisentangler|>getinspace|>unitcellfock, secondgmeracrystal)
# newcouriercrystalfock = getcrystalfock(localwanniercourier|>getinspace|>unitcellfock, Crystal(courierseeds|>offsetofmodes,secondgmeracrystal.sizes))
# newfilledcrystalfock = getcrystalfock(localwannierfilled|>getinspace|>unitcellfock, Crystal(filledseeds|>offsetofmodes,secondgmeracrystal.sizes))
# newemptycrystalfock = getcrystalfock(localwannierempty|>getinspace|>unitcellfock, Crystal(emptyseeds|>offsetofmodes,secondgmeracrystal.sizes))

newlocalrestrict = fourier(newblockedcrystalfock, redefinelocaldisentangler|>getinspace) * (secondgmeracrystal|>vol|>sqrt)
# newlocalcourierrestrict = fourier(newcouriercrystalfock, localwanniercourier|>getinspace) * (secondgmeracrystal|>vol|>sqrt)
# newlocalfilledrestrict = fourier(newfilledcrystalfock, localwannierfilled|>getinspace) * (secondgmeracrystal|>vol|>sqrt)
# newlocalemptyrestrict = fourier(newemptycrystalfock, localwannierempty|>getinspace) * (secondgmeracrystal|>vol|>sqrt)

globaldisentangler = broadcast(*,(localrestrict*redefinelocaldisentangler), newlocalrestrict')
globalcourierisometry = broadcast(*,(localrestrict*localwanniercourier), newlocalcourierrestrict')
globalfilledisometry = broadcast(*,(localrestrict*localwannierfilled), newlocalfilledrestrict')
globalemptyisometry = broadcast(*,(localrestrict*localwannierempty), newlocalemptyrestrict')

rotcorrelations3 = globaldisentangler'*rotcorrelations2*globaldisentangler
rotcorrelations3|>crystalspectrum|>visualize
rotH3 = globaldisentangler'*rotH2*globaldisentangler

thirdgmeracrystalfock = rotcorrelations3|>getoutspace
thirdgmeracrystal::Crystal = thirdgmeracrystalfock|>getcrystal
thirdgmeraspace::RealSpace = thirdgmeracrystal|>getspace

thirdcenter = [0,0] ∈ thirdgmeraspace
thirdhexagonalregion = gethexagonalregion(rot = refrot,crystal=secondgmeracrystal, center=thirdcenter, metricspace=secondgmeraspace)
thirdhexagonalregionfock = quantize(thirdhexagonalregion,noofflavourpermode)
thirdhexagonalmodes = Subset(mode for mode in thirdhexagonalregionfock)

# @info "Computing local correlations..."
localrestrict = fourier(thirdgmeracrystalfock, thirdhexagonalregionfock) / (thirdgmeracrystal|>vol|>sqrt)
localcorrelations = localrestrict'*rotcorrelations2*localrestrict

entanglementcontourinfo = entanglementcontour(localcorrelations)
groupedandsortedmodeswifecontour = sortgroupdictwifvalue(entanglementcontourinfo,false)

groupedandsortedmodeswifecontour[1][2]|>offsetofmodes|>visualize

localcorrelations|>eigspech|>visualize






couriercorrelations = globalcourierisometry'*rotH2*globalcourierisometry
couriercorrelations|>getinspace|>unitcellfock
couriercorrelations|>getinspace|>getcrystal|>getunitcell|>visualize

blockedH|>crystalspectrum|>visualize