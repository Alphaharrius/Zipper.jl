using Zipper
using LinearAlgebra

function generatescaledAsublatticeunitcellmodes(scale::Scale,size::Number)
    triangular = RealSpace([1 0; 1/2 sqrt(3)/2]')
    pa = [1/3, 0] ∈ triangular
    pc = [0, 2/3] ∈ triangular
    pe = [2/3, 1/3] ∈ triangular
    unitcell = Subset(pa, pc, pe)
    crystal = Crystal(unitcell, [size, size])
    scaledcrystal = scale*crystal
    return sum([Subset(m|>removeattr(:r)) for m in quantize(scaledcrystal|>getunitcell,1)])
end

function generatescaledBsublatticeunitcellmodes(scale::Scale,size::Number)
    triangular = RealSpace([1 0; 1/2 sqrt(3)/2]')
    pb = [0, 1/3] ∈ triangular
    pd = [1/3, 2/3] ∈ triangular
    pf = [2/3, 0] ∈ triangular
    unitcell = Subset(pb, pd, pf)
    crystal = Crystal(unitcell, [size, size])
    scaledcrystal = scale*crystal
    return sum([Subset(m|>removeattr(:r)) for m in quantize(scaledcrystal|>getunitcell,1)])
end

function distinguishABsublatticemodesforhexagon(modes::Subset{Mode},Arefmodes::Subset{Mode},Brefmodes::Subset{Mode})
    modesdict = Dict(m|>removeattr(:r)=> m for m in modes )
    Amodes = sum([Subset(modesdict[am]) for am in Arefmodes if haskey(modesdict,am)])
    Bmodes = sum([Subset(modesdict[bm]) for bm in Brefmodes if haskey(modesdict,bm)]) 
    return Amodes,Bmodes
end

function sortgroupdictwifvalue(dict::Dict{Mode, Number},rev::Bool)
    sortedtupledata = sort([(value,key) for (key,value) in dict],rev=rev,by=first)
    refvalue = sortedtupledata[1][2]|>getattr(:eigenindex)
    result = []
    evals = []
    subresult = Subset(sortedtupledata[1][2])
    for pair in sortedtupledata
        if pair[2]|>getattr(:eigenindex)==refvalue
            subresult = subresult + Subset(pair[2])
            append!(evals,pair[1])
        else
            append!(result,[subresult])
            refvalue = pair[2]|>getattr(:eigenindex)
            evals = []
            subresult = Subset(pair[2])
        end
    end
    append!(result,[subresult])

    return result
end

triangular = RealSpace([1 0; 1/2 sqrt(3)/2]')
kspace = convert(MomentumSpace, triangular)

pa = [1/3, 0] ∈ triangular
pb = [0, 1/3] ∈ triangular
pc = [0, 2/3] ∈ triangular
pd = [1/3, 2/3] ∈ triangular
pe = [2/3, 1/3] ∈ triangular
pf = [2/3, 0] ∈ triangular

pg = (pa + pb + pc + pd + pe + pf) / 6
spatialsnappingcalibration((pa, pb, pc, pd, pe, pf, pg))

unitcell = Subset(pa, pb, pc, pd, pe, pf)
visualize(unitcell)
# visualize(Subset( pd, pe, pf))
crystal = Crystal(unitcell, [20, 20])
reciprocalhashcalibration(crystal.sizes)
# visualize(crystal|>sitepoints)

modes::Subset{Mode} = quantize(unitcell, 1)|>orderedmodes
m0, m1, m2, m3, m4, m5 = members(modes)
m5|>getattr(:b)

t_a = -0.3 + 0im
t_b = 0.3 + 0im
tₙ = -1 + 0im
tₕ = -(0.25/(3*sqrt(3)))im

onsite = [
    (m0, m0) => t_a,
    (m1, m1) => t_b,
    (m2, m2) => t_a,
    (m3, m3) => t_b,
    (m4, m4) => t_a,
    (m5, m5) => t_b
]

nearestneighbor = [
    (m1, m0) => tₙ,
    (m1, m4|> setattr(:r => [-1, 0] ∈ triangular)) => tₙ,
    (m1, m2) => tₙ,
    (m3, m2) => tₙ,
    (m3, m0|> setattr(:r => [0, 1] ∈ triangular)) => tₙ,
    (m3, m4) => tₙ,
    (m5, m4) => tₙ,
    (m5, m2|> setattr(:r => [1, -1] ∈ triangular)) => tₙ,
    (m5, m0) => tₙ]

haldaneterm = [
    (m1, m3|> setattr(:r => [-1, 0] ∈ triangular)) => tₕ,
    (m1, m3) => tₕ,
    (m1, m3|> setattr(:r => [0, -1] ∈ triangular)) => tₕ,
    (m3, m5|> setattr(:r => [-1, 1] ∈ triangular)) => tₕ,
    (m3, m5|> setattr(:r => [0, 1] ∈ triangular)) => tₕ,
    (m3, m5) => tₕ,
    (m5, m1) => tₕ,
    (m5, m1|> setattr(:r => [1, 0] ∈ triangular)) => tₕ,
    (m5, m1|> setattr(:r => [1, -1] ∈ triangular)) => tₕ,
    (m0, m4|> setattr(:r => [-1, 0] ∈ triangular)) => -tₕ,
    (m0, m4) => -tₕ,
    (m0, m4|> setattr(:r => [0, -1] ∈ triangular)) => -tₕ,
    (m2, m0|> setattr(:r => [-1, 1] ∈ triangular)) => -tₕ,
    (m2, m0|> setattr(:r => [0, 1] ∈ triangular)) => -tₕ,
    (m2, m0) => -tₕ,
    (m4, m2) => -tₕ,
    (m4, m2|> setattr(:r => [1, 0] ∈ triangular)) => -tₕ,
    (m4, m2|> setattr(:r => [1, -1] ∈ triangular)) => -tₕ]

bonds::FockMap = bondmap([onsite..., nearestneighbor...,haldaneterm...])
energyspectrum = @time computeenergyspectrum(bonds, crystal=crystal)
# energyspectrum|>visualize

groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=3)
groundstateprojector = groundstates|>crystalprojector
correlations = idmap(groundstateprojector|>getoutspace) - groundstateprojector
H = CrystalFockMap(energyspectrum)

@info("Starting RG...")
crystalfock = correlations|>getoutspace

scale = Scale([5 0; 0 5], crystalfock|>getcrystal|>getspace)
@info("Performing rgblocking...")
@info("Generating rgblocking transformation...")
blocker = @time scale * crystalfock
@info("Performing rgblocking on correlations...")
blockedcorrelations = @time blocker * correlations * blocker'
blockedcrystalfock = blockedcorrelations|>getoutspace
blockedcrystal::Crystal = blockedcrystalfock|>getcrystal
blockedspace::RealSpace = blockedcrystal|>getspace
blockedcrystal|>getunitcell|>visualize

@info("Computing local correlations...")
firstcenter = [0,0] ∈ blockedspace
firsthexagonalregion = gethexagonalregion(crystal=blockedcrystal, center=firstcenter, metricspace=blockedspace)
visualize(firsthexagonalregion)
firsthexagonalregionfock = quantize(firsthexagonalregion,1)

Asublatticeunitcellmodes = generatescaledAsublatticeunitcellmodes(scale,20)
Bsublatticeunitcellmodes = generatescaledBsublatticeunitcellmodes(scale,20)

firstAsublatticethexagonalregionmodes,firstBsublatticethexagonalregionmodes =  distinguishABsublatticemodesforhexagon(firsthexagonalregionfock|>orderedmodes,Asublatticeunitcellmodes,Bsublatticeunitcellmodes)
firstAsublatticethexagonalregionfock = RegionFock(firstAsublatticethexagonalregionmodes)
firstBsublatticethexagonalregionfock = RegionFock(firstBsublatticethexagonalregionmodes)

shiftedfirstcenter1 = [1,0] ∈ blockedspace
shiftedfirsthexagonalregion1 = firsthexagonalregion.+shiftedfirstcenter1
shiftedfirsthexagonalregion1fock = quantize(shiftedfirsthexagonalregion1,1)
shiftedfirstAsublatticethexagonalregionmodes1,shiftedfirstBsublatticethexagonalregionmodes1 =  distinguishABsublatticemodesforhexagon(shiftedfirsthexagonalregion1fock|>orderedmodes,Asublatticeunitcellmodes,Bsublatticeunitcellmodes)
shiftedfirstAsublatticethexagonalregionfock1 = RegionFock(shiftedfirstAsublatticethexagonalregionmodes1)
shiftedfirstBsublatticethexagonalregionfock1 = RegionFock(shiftedfirstBsublatticethexagonalregionmodes1)

shiftedfirstcenter2 = [0,1] ∈ blockedspace
shiftedfirsthexagonalregion2 = firsthexagonalregion.+shiftedfirstcenter2
shiftedfirsthexagonalregion2fock = quantize(shiftedfirsthexagonalregion2,1)
shiftedfirstAsublatticethexagonalregionmodes2,shiftedfirstBsublatticethexagonalregionmodes2 =  distinguishABsublatticemodesforhexagon(shiftedfirsthexagonalregion2fock|>orderedmodes,Asublatticeunitcellmodes,Bsublatticeunitcellmodes)
shiftedfirstAsublatticethexagonalregionfock2 = RegionFock(shiftedfirstAsublatticethexagonalregionmodes2)
shiftedfirstBsublatticethexagonalregionfock2 = RegionFock(shiftedfirstBsublatticethexagonalregionmodes2)

shiftedfirstcenter3 = [1,1] ∈ blockedspace
shiftedfirsthexagonalregion3 = firsthexagonalregion.+shiftedfirstcenter3
shiftedfirsthexagonalregion3fock = quantize(shiftedfirsthexagonalregion3,1)
shiftedfirstAsublatticethexagonalregionmodes3,shiftedfirstBsublatticethexagonalregionmodes3 =  distinguishABsublatticemodesforhexagon(shiftedfirsthexagonalregion3fock|>orderedmodes,Asublatticeunitcellmodes,Bsublatticeunitcellmodes)
shiftedfirstAsublatticethexagonalregionfock3 = RegionFock(shiftedfirstAsublatticethexagonalregionmodes3)
shiftedfirstBsublatticethexagonalregionfock3 = RegionFock(shiftedfirstBsublatticethexagonalregionmodes3)

allregion = firsthexagonalregion+shiftedfirsthexagonalregion1+shiftedfirsthexagonalregion2+shiftedfirsthexagonalregion3
visualize(allregion)
if intersect(allregion,blockedcrystal|>getunitcell) == blockedcrystal|>getunitcell
    @info("cover the whole unitcell")
else
    @error("the allregion cannot cover the whole unitcell ")
end

localcorrelations = regioncorrelations(blockedcorrelations,firsthexagonalregionfock)
localspectrum = localcorrelations|>eigspech
filleddict,emptydict = separatefilledandemptymodes(localspectrum)

rankedandgroupedfilledmodes = sortgroupdictwifvalue(filleddict,false)
rankedandgroupedemptymodes = sortgroupdictwifvalue(emptydict,true)

chosenfilledmodes = sum(rankedandgroupedfilledmodes[1:24])
chosenemptymodes = sum(rankedandgroupedemptymodes[1:24])

rsmodedict = Dict(mode=>norm(euclidean((mode|>getattr(:b))+(mode|>getattr(:r))-firstcenter)) for mode in localspectrum |> geteigenvectors|>getoutspace|>orderedmodes)
sortedrsmode = sortgroupdictwifdist(rsmodedict,false)
frozenrsmodes = sum([pair[2] for pair in sortedrsmode[1:6]])
frozenrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in frozenrsmodes)
visualize(frozenrsmodesoffset)
filledrsmodes,emptyrsmodes = distinguishABsublatticemodesforhexagon(frozenrsmodes,Asublatticeunitcellmodes,Bsublatticeunitcellmodes)

courierrsmodes = sum([pair[2] for pair in sortedrsmode[7:end]])
courierrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in courierrsmodes)
visualize(courierrsmodesoffset)

courierseeds = idmap(firsthexagonalregionfock, firsthexagonalregionfock)[:,FockSpace(mode for mode in courierrsmodes)]
filledseeds = idmap(firsthexagonalregionfock, firsthexagonalregionfock)[:,FockSpace(mode for mode in filledrsmodes)]
emptyseeds = idmap(firsthexagonalregionfock, firsthexagonalregionfock)[:,FockSpace(mode for mode in emptyrsmodes)]
frozenseeds = idmap(firsthexagonalregionfock, firsthexagonalregionfock)[:,FockSpace(mode for mode in frozenrsmodes)]

# rankedandgroupedoverlapwiffilledmodes = sort([calculateavgoverlapofmodes(modes,localspectrum,frozenrsmodes) for modes in rankedandgroupedfilledmodes],rev=true,by=first)
# rankedandgroupedoverlapwifemptymodes = sort([calculateavgoverlapofmodes(modes,localspectrum,frozenrsmodes) for modes in rankedandgroupedemptymodes],rev=true,by=first)

# chosenfilledmodes = sum([pair[2] for pair in rankedandgroupedoverlapwiffilledmodes[1:24]])
# chosenemptymodes = sum([pair[2] for pair in rankedandgroupedoverlapwifemptymodes[1:24]])
chosenfrozenmodes = chosenfilledmodes+chosenemptymodes
chosencouriermodes = (localspectrum |> geteigenmodes)-chosenfrozenmodes

localisofilled = columns(localspectrum |> geteigenvectors, FockSpace(chosenfilledmodes))
localisoempty = columns(localspectrum |> geteigenvectors, FockSpace(chosenemptymodes))
localisoforzen = columns(localspectrum |> geteigenvectors, FockSpace(chosenfrozenmodes))
localisocourier = columns(localspectrum |> geteigenvectors, FockSpace(chosencouriermodes))

wannierfilled = localwannierization(localisofilled, filledseeds)
wannierempty = localwannierization(localisoempty, emptyseeds)
wannierfrozen = localwannierization(localisoforzen, frozenseeds)
wanniercourier = localwannierization(localisocourier, courierseeds)

visualize(wannierfilled'*localcorrelations*wannierfilled|>eigspech)
visualize(wannierempty'*localcorrelations*wannierempty|>eigspech)
visualize(wanniercourier'*localcorrelations*wanniercourier|>eigspech)
visualize(wannierfrozen'*localcorrelations*wannierfrozen|>eigspech)

localwannierresults =  Dict(:wannierfilled => wannierfilled,:wannierempty => wannierempty,
                    :wannierfrozen => wannierfrozen, :wanniercourier => wanniercourier,
                    :filledseeds => filledseeds,:emptyseeds => emptyseeds,
                    :frozenseeds => frozenseeds,:courierseeds => courierseeds)

wannierinfos =  Dict(firsthexagonalregionfock=>localwannierresults)

firstshiftedhexagonalregionfocklist = [(shiftedfirstcenter1,shiftedfirsthexagonalregion1fock), (shiftedfirstcenter2,shiftedfirsthexagonalregion2fock), (shiftedfirstcenter3,shiftedfirsthexagonalregion3fock)]

for (shifted,hexagonalregionfock) in firstshiftedhexagonalregionfocklist
    localcorrelations = regioncorrelations(blockedcorrelations,hexagonalregionfock)
    localspectrum = localcorrelations|>eigspech
    filleddict,emptydict = separatefilledandemptymodes(localspectrum)
    rankedandgroupedfilledmodes = sortgroupdictwifvalue(filleddict,false)
    rankedandgroupedemptymodes = sortgroupdictwifvalue(emptydict,true)

    chosenfilledmodes = sum(rankedandgroupedfilledmodes[1:24])
    chosenemptymodes = sum(rankedandgroupedemptymodes[1:24])
    
    rsmodedict = Dict(mode=>norm(euclidean((mode|>getattr(:b))+(mode|>getattr(:r))-(firstcenter+shifted))) for mode in localspectrum |> geteigenvectors|>getoutspace|>orderedmodes)
    sortedrsmode = sortgroupdictwifdist(rsmodedict,false)
    frozenrsmodes = sum([pair[2] for pair in sortedrsmode[1:6]])
    frozenrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in frozenrsmodes)
    filledrsmodes,emptyrsmodes = distinguishABsublatticemodesforhexagon(frozenrsmodes,Asublatticeunitcellmodes,Bsublatticeunitcellmodes)
    courierrsmodes = sum([pair[2] for pair in sortedrsmode[7:end]])
    courierrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in courierrsmodes)

    courierseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in courierrsmodes)]
    filledseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in filledrsmodes)]
    emptyseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in emptyrsmodes)]
    frozenseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in frozenrsmodes)]


    # rankedandgroupedoverlapwiffilledmodes = sort([calculateavgoverlapofmodes(modes,localspectrum,frozenrsmodes) for modes in rankedandgroupedfilledmodes],rev=true,by=first)
    # rankedandgroupedoverlapwifemptymodes = sort([calculateavgoverlapofmodes(modes,localspectrum,frozenrsmodes) for modes in rankedandgroupedemptymodes],rev=true,by=first)
    
    # chosenfilledmodes = sum([pair[2] for pair in rankedandgroupedoverlapwiffilledmodes[1:24]])
    # chosenemptymodes = sum([pair[2] for pair in rankedandgroupedoverlapwifemptymodes[1:24]])
    chosenfrozenmodes = chosenfilledmodes+chosenemptymodes
    chosencouriermodes = (localspectrum |> geteigenmodes)-chosenfrozenmodes
    
    localisofilled = columns(localspectrum |> geteigenvectors, FockSpace(chosenfilledmodes))
    localisoempty = columns(localspectrum |> geteigenvectors, FockSpace(chosenemptymodes))
    localisoforzen = columns(localspectrum |> geteigenvectors, FockSpace(chosenfrozenmodes))
    localisocourier = columns(localspectrum |> geteigenvectors, FockSpace(chosencouriermodes))

    wannierfilled = localwannierization(localisofilled, filledseeds)
    wannierempty = localwannierization(localisoempty, emptyseeds)
    wannierfrozen = localwannierization(localisoforzen, frozenseeds)
    wanniercourier = localwannierization(localisocourier, courierseeds)
    
    shiftedlocalwannierresults =  Dict(:wannierfilled => wannierfilled,:wannierempty => wannierempty,
                                :wannierfrozen => wannierfrozen, :wanniercourier => wanniercourier,
                                :filledseeds => filledseeds,:emptyseeds => emptyseeds,
                                :frozenseeds => frozenseeds,:courierseeds => courierseeds)

    wannierinfos[hexagonalregionfock] =  shiftedlocalwannierresults
end 

firsthexagonalregionfocklist = [firsthexagonalregionfock, shiftedfirsthexagonalregion1fock, shiftedfirsthexagonalregion2fock, shiftedfirsthexagonalregion3fock]

extendedwannierizedfilled =  sum(wannierinfos[regionfock][:wannierfilled] for regionfock in firsthexagonalregionfocklist)
extendedwannierizedempty =  sum(wannierinfos[regionfock][:wannierempty] for regionfock in firsthexagonalregionfocklist)
extendedwannierizedfrozen =  sum(wannierinfos[regionfock][:wannierfrozen] for regionfock in firsthexagonalregionfocklist)
extendedwannierizedcourier =  sum(wannierinfos[regionfock][:wanniercourier] for regionfock in firsthexagonalregionfocklist)
            
origin = [0, 0] ∈ blockedspace
refunictcellfockfilled = FockSpace(Subset(mode for mode in extendedwannierizedfilled |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
refunictcellfockempty = FockSpace(Subset(mode for mode in extendedwannierizedempty |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
refunictcellfockfrozen = FockSpace(Subset(mode for mode in extendedwannierizedfrozen |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
refunictcellfockcourier = FockSpace(Subset(mode for mode in extendedwannierizedcourier |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))

wannierfilledisometry = globalwannierfunction(blockedcorrelations,extendedwannierizedfilled[:,refunictcellfockfilled])
wannieremptyisometry = globalwannierfunction(blockedcorrelations,extendedwannierizedempty[:,refunictcellfockempty])
wannierforzenisometry = globalwannierfunction(blockedcorrelations,extendedwannierizedfrozen[:,refunictcellfockfrozen])
wanniercourierisometry = globalwannierfunction(blockedcorrelations,extendedwannierizedcourier[:,refunictcellfockcourier])

couriercorrelations = wanniercourierisometry' * blockedcorrelations * wanniercourierisometry
couriercorrelations|>eigspech|>visualize
couriercorrelationspectrum = couriercorrelations |> crystalspectrum
purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
purifiedcouriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap

firstrgedcorrelations = purifiedcouriercorrelations
firstrgedcrystalfock = purifiedcouriercorrelations|>getoutspace
firstrgedcrystal::Crystal = firstrgedcrystalfock|>getcrystal
firstrgedspace::RealSpace = firstrgedcrystal|>getspace

@info("Computing local correlations...")
secondcenter = [2/3,-1/3] ∈ firstrgedspace
secondhexagonalregion = gethexagonalregion(crystal=firstrgedcrystal, center=secondcenter, metricspace=firstrgedspace)
visualize(secondhexagonalregion)
secondhexagonalregionfock = quantize(secondhexagonalregion,1)

shiftedsecondcenter1 = [0,1] ∈ firstrgedspace
shiftedsecondhexagonalregion1 = secondhexagonalregion.+shiftedsecondcenter1
shiftedsecondhexagonalregion1fock = quantize(shiftedsecondhexagonalregion1,1)

shiftedsecondcenter2 = [-1,1] ∈ firstrgedspace
shiftedsecondhexagonalregion2 = secondhexagonalregion.+shiftedsecondcenter2
shiftedsecondhexagonalregion2fock = quantize(shiftedsecondhexagonalregion2,1)


allregion2 = secondhexagonalregion+shiftedsecondhexagonalregion1+shiftedsecondhexagonalregion2
visualize(allregion2)
if intersect(allregion2,firstrgedcrystal|>getunitcell) == firstrgedcrystal|>getunitcell
    @info("cover the whole unitcell")
else
    @error("the allregion cannot cover the whole unitcell ")
end

localcorrelations = regioncorrelations(firstrgedcorrelations,secondhexagonalregionfock)
localspectrum = localcorrelations|>eigspech
display(visualize(localspectrum))
filleddict,emptydict = separatefilledandemptymodes(localspectrum)
rankedandgroupedfilledmodes = sortgroupdictwifvalue(filleddict,false)
rankedandgroupedemptymodes = sortgroupdictwifvalue(emptydict,true)

rsmodedict = Dict(mode=>norm(euclidean((mode|>getattr(:b))+(mode|>getattr(:r))-secondcenter)) for mode in localspectrum |> geteigenvectors|>getoutspace|>orderedmodes)
sortedrsmode = sortgroupdictwifdist(rsmodedict,false)
frozenrsmodes = sum([pair[2] for pair in sortedrsmode[1:8]])
frozenrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in frozenrsmodes)
visualize(frozenrsmodesoffset)
filledrsmodes,emptyrsmodes = distinguishABsublatticemodesforhexagon(frozenrsmodes,Asublatticeunitcellmodes,Bsublatticeunitcellmodes)

courierrsmodes = sum([pair[2] for pair in sortedrsmode[9:end]])
courierrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in courierrsmodes)
visualize(courierrsmodesoffset)

courierseeds = idmap(secondhexagonalregionfock, secondhexagonalregionfock)[:,FockSpace(mode for mode in courierrsmodes)]
filledseeds = idmap(secondhexagonalregionfock, secondhexagonalregionfock)[:,FockSpace(mode for mode in filledrsmodes)]
emptyseeds = idmap(secondhexagonalregionfock, secondhexagonalregionfock)[:,FockSpace(mode for mode in emptyrsmodes)]
frozenseeds = idmap(secondhexagonalregionfock, secondhexagonalregionfock)[:,FockSpace(mode for mode in frozenrsmodes)]

rankedandgroupedoverlapwiffilledmodes = sort([calculateavgoverlapofmodes(modes,localspectrum,frozenrsmodes) for modes in rankedandgroupedfilledmodes],rev=true,by=first)
rankedandgroupedoverlapwifemptymodes = sort([calculateavgoverlapofmodes(modes,localspectrum,frozenrsmodes) for modes in rankedandgroupedemptymodes],rev=true,by=first)

# for r in range(1,48)
#     println([(localspectrum|>geteigenvalues)[m] for m in rankedandgroupedoverlapwiffilledmodes[r:r][1][2]])
# end
# [(localspectrum|>geteigenvalues)[m] for m in rankedandgroupedoverlapwiffilledmodes[32:32][1][2]]

# chosenfilledmodes = sum([pair[2] for pair in rankedandgroupedoverlapwiffilledmodes[1:23]])+sum([pair[2] for pair in rankedandgroupedoverlapwiffilledmodes[25:25]])
chosenfilledmodes = sum([pair[2] for pair in rankedandgroupedoverlapwiffilledmodes[1:23]])+ref
ref = sum([pair[2] for pair in rankedandgroupedoverlapwiffilledmodes[25:25]])
# calculateavgoverlapofmodes(chosenfilledmodes,localspectrum,filledrsmodes)
# calculateavgoverlapofmodes(chosenfilledmodes,localspectrum,emptyrsmodes)
# sort([calculateavgoverlapofmodes(modes,localspectrum,frozenrsmodes) for modes in rankedandgroupedfilledmodes],rev=true,by=first)

chosenemptymodes = sum([pair[2] for pair in rankedandgroupedoverlapwifemptymodes[1:23]])+sum([pair[2] for pair in rankedandgroupedoverlapwifemptymodes[25:25]])
chosenfrozenmodes = chosenfilledmodes+chosenemptymodes
chosencouriermodes = (localspectrum |> geteigenmodes)-chosenfrozenmodes
columns(localspectrum |> geteigenvectors, FockSpace(ref))'*filledseeds
localisofilled = columns(localspectrum |> geteigenvectors, FockSpace(chosenfilledmodes))

U, Σ, Vt = svd(filledseeds'*localisofilled)
Σ

norm((U*U')*(columns(localspectrum |> geteigenvectors, FockSpace(ref))'*filledseeds)'-(columns(localspectrum |> geteigenvectors, FockSpace(ref))'*filledseeds)')
(columns(localspectrum |> geteigenvectors, FockSpace(ref))'*filledseeds)'
norm(((localisofilled'*filledseeds)'*(localisofilled'*filledseeds))*(columns(localspectrum |> geteigenvectors, FockSpace(ref))'*filledseeds)'-(columns(localspectrum |> geteigenvectors, FockSpace(ref))'*filledseeds)')
norm((localisofilled'*filledseeds)*(columns(localspectrum |> geteigenvectors, FockSpace(ref))'*filledseeds)')
localisoempty = columns(localspectrum |> geteigenvectors, FockSpace(chosenemptymodes))
localisofrozen = columns(localspectrum |> geteigenvectors, FockSpace(chosenfrozenmodes))
localisocourier = columns(localspectrum |> geteigenvectors, FockSpace(chosencouriermodes))

wannierfilled = localwannierization(localisofilled, filledseeds)
wannierempty = localwannierization(localisoempty, emptyseeds)
wannierfrozen = localwannierization(localisofrozen, frozenseeds)
wanniercourier = localwannierization(localisocourier, courierseeds)

chosen = sum([Subset(pair[2]) for pair in sort([(norm(columns((localspectrum|>geteigenvectors),FockSpace(m))'*wannierfilledref)^2,m) for m in chosenfilledmodes],rev=true,by=first)[1:24]])
localisofilled = columns(localspectrum |> geteigenvectors, FockSpace(chosen))
wannierfilledref = localwannierization(localisofilled, filledseeds)

(wannierfilled'*localcorrelations*wannierfilled)|>eigspech|>visualize
(wannierfrozen'*localcorrelations*wannierfrozen)|>eigspech|>visualize
(wanniercourier'*localcorrelations*wanniercourier)|>eigspech|>visualize

localwannierresults =  Dict(:wannierfrozen => wannierfrozen, :wanniercourier => wanniercourier,
                    :frozenseeds => frozenseeds,:courierseeds => courierseeds)

wannierinfos =  Dict(secondhexagonalregionfock=>localwannierresults)

secondshiftedhexagonalregionfocklist = [(shiftedsecondcenter1,shiftedsecondhexagonalregion1fock), (shiftedsecondcenter2,shiftedsecondhexagonalregion2fock)]

for (shifted,hexagonalregionfock) in secondshiftedhexagonalregionfocklist
    localcorrelations = regioncorrelations(firstrgedcorrelations,hexagonalregionfock)
    localspectrum = localcorrelations|>eigspech
    filleddict,emptydict = separatefilledandemptymodes(localspectrum)
    rankedandgroupedfilledmodes = sortgroupdictwifvalue(filleddict,false)
    rankedandgroupedemptymodes = sortgroupdictwifvalue(emptydict,true)
    
    rsmodedict = Dict(mode=>norm(euclidean((mode|>getattr(:b))+(mode|>getattr(:r))-(secondcenter+shifted))) for mode in localspectrum |> geteigenvectors|>getoutspace|>orderedmodes)
    sortedrsmode = sortgroupdictwifdist(rsmodedict,false)
    frozenrsmodes = sum([pair[2] for pair in sortedrsmode[1:8]])
    frozenrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in frozenrsmodes)
    courierrsmodes = sum([pair[2] for pair in sortedrsmode[9:end]])
    courierrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in courierrsmodes)

    courierseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in courierrsmodes)]
    frozenseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in frozenrsmodes)]


    rankedandgroupedoverlapwifmodes = sort([calculateavgoverlapofmodes(pair[1]+pair[2],localspectrum,frozenrsmodes) for pair in zip(rankedandgroupedfilledmodes,rankedandgroupedemptymodes)],rev=true,by=first)
    chosenfrozenmodesref = sum([pair[2] for pair in rankedandgroupedoverlapwifmodes[1:27]])
    localisoforzenref = columns(localspectrum |> geteigenvectors, FockSpace(chosenfrozenmodesref))
    wannierfrozenref = localwannierization(localisoforzenref, frozenseeds)
    overlap = columns(localspectrum|>geteigenvectors,FockSpace(chosenfrozenmodesref))'*wannierfrozenref
    sorted = sort([(norm(rows(overlap,FockSpace(mode)))^2,mode) for mode in overlap|>getoutspace],by=first,rev=true)

    chosenfrozenmodes = sum([Subset(data[2]) for data in sorted[1:48]])
    chosencouriermodes = (localspectrum |> geteigenmodes)-chosenfrozenmodes
    
    localisoforzen = columns(localspectrum |> geteigenvectors, FockSpace(chosenfrozenmodes))
    localisocourier = columns(localspectrum |> geteigenvectors, FockSpace(chosencouriermodes))
    wannierfrozen = localwannierization(localisoforzen, frozenseeds)
    wanniercourier = localwannierization(localisocourier, courierseeds)
    
    shiftedlocalwannierresults =  Dict(:wannierfrozen => wannierfrozen, :wanniercourier => wanniercourier,
                    :frozenseeds => frozenseeds,:courierseeds => courierseeds)

    wannierinfos[hexagonalregionfock] =  shiftedlocalwannierresults
end 

secondhexagonalregionfocklist = [secondhexagonalregionfock, shiftedsecondhexagonalregion1fock, shiftedsecondhexagonalregion2fock]

extendedwannierizedfrozen =  sum(wannierinfos[regionfock][:wannierfrozen] for regionfock in secondhexagonalregionfocklist)
extendedwannierizedcourier =  sum(wannierinfos[regionfock][:wanniercourier] for regionfock in secondhexagonalregionfocklist)
            
origin = [0, 0] ∈ blockedspace
refunictcellfockfrozen = FockSpace(Subset(mode for mode in extendedwannierizedfrozen |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
refunictcellfockcourier = FockSpace(Subset(mode for mode in extendedwannierizedcourier |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))

wannierforzenisometry = globalwannierfunction(firstrgedcorrelations,extendedwannierizedfrozen[:,refunictcellfockfrozen])
wanniercourierisometry = globalwannierfunction(firstrgedcorrelations,extendedwannierizedcourier[:,refunictcellfockcourier])

frozencorrelations = wannierforzenisometry' * firstrgedcorrelations * wannierforzenisometry
frozencorrelations|>eigspech|>visualize

couriercorrelations2 = wanniercourierisometry' * firstrgedcorrelations * wanniercourierisometry
couriercorrelations2|>eigspech|>visualize

couriercorrelationspectrum2 = couriercorrelations2 |> crystalspectrum
purifiedcorrelationspectrum2 = couriercorrelationspectrum2 |> roundingpurification
purifiedcouriercorrelations2 = purifiedcorrelationspectrum2 |> CrystalFockMap

secondrgedcorrelations = purifiedcouriercorrelations2
secondrgedcrystalfock = purifiedcouriercorrelations2|>getoutspace
secondrgedcrystal::Crystal = secondrgedcrystalfock|>getcrystal
secondrgedspace::RealSpace = secondrgedcrystal|>getspace
secondrgedcrystal|>getunitcell|>visualize

@info("Computing local correlations...")
thirdcenter = [1/3,1/3] ∈ blockedspace
thirdhexagonalregion = gethexagonalregion(crystal=secondrgedcrystal, center=thirdcenter, metricspace=secondrgedspace)
thirdhexagonalregionfock = quantize(thirdhexagonalregion,1)
visualize(thirdhexagonalregion)

allregion3 = thirdhexagonalregion
visualize(allregion3)
visualize(secondrgedcrystal|>getunitcell)
if intersect(allregion3,firstrgedcrystal|>getunitcell) == secondrgedcrystal|>getunitcell
    @info("cover the whole unitcell")
else
    @error("the allregion cannot cover the whole unitcell ")
end

localcorrelations = regioncorrelations(secondrgedcorrelations,thirdhexagonalregionfock)
localspectrum = localcorrelations|>eigspech
display(visualize(localspectrum))
filleddict,emptydict = separatefilledandemptymodes(localspectrum)
rankedandgroupedfilledmodes = sortgroupdictwifvalue(filleddict,false)
rankedandgroupedemptymodes = sortgroupdictwifvalue(emptydict,true)

rsmodedict = Dict(mode=>norm(euclidean((mode|>getattr(:b))+(mode|>getattr(:r))-thirdcenter)) for mode in localspectrum |> geteigenvectors|>getoutspace|>orderedmodes)
sortedrsmode = sortgroupdictwifdist(rsmodedict,false)
frozenrsmodes = sum([pair[2] for pair in sortedrsmode[1:7]])
frozenrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in frozenrsmodes)
visualize(frozenrsmodesoffset)
courierrsmodes = sum([pair[2] for pair in sortedrsmode[8:end]])
courierrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in courierrsmodes)
visualize(courierrsmodesoffset)
courierseeds = idmap(thirdhexagonalregionfock, thirdhexagonalregionfock)[:,FockSpace(mode for mode in courierrsmodes)]
frozenseeds = idmap(thirdhexagonalregionfock, thirdhexagonalregionfock)[:,FockSpace(mode for mode in frozenrsmodes)]

rankedandgroupedoverlapwifmodes = sort([calculateavgoverlapofmodes(pair[1]+pair[2],localspectrum,frozenrsmodes) for pair in zip(rankedandgroupedfilledmodes,rankedandgroupedemptymodes)],rev=true,by=first)

chosenfrozenmodes = sum([pair[2] for pair in rankedandgroupedoverlapwifmodes[1:21]])
chosencouriermodes = (localspectrum |> geteigenmodes)-chosenfrozenmodes


localisoforzen = columns(localspectrum |> geteigenvectors, FockSpace(chosenfrozenmodes))
localisocourier = columns(localspectrum |> geteigenvectors, FockSpace(chosencouriermodes))
wannierfrozen = localwannierization(localisoforzen, frozenseeds)
wanniercourier = localwannierization(localisocourier, courierseeds)

visualize(wannierfrozen'*localcorrelations*wannierfrozen|>eigspech)

localwannierresults =  Dict(:wannierfrozen => wannierfrozen, :wanniercourier => wanniercourier,
                :frozenseeds => frozenseeds,:courierseeds => courierseeds)

wannierinfos =  Dict(thirdhexagonalregionfock=>localwannierresults)

wannierforzenisometry = globalwannierfunction(secondrgedcorrelations,wannierinfos[thirdhexagonalregionfock][:wannierfrozen])
wanniercourierisometry = globalwannierfunction(secondrgedcorrelations,wannierinfos[thirdhexagonalregionfock][:wanniercourier])

frozencorrelations = wannierforzenisometry' * secondrgedcorrelations * wannierforzenisometry
frozencorrelations|>eigspech|>visualize

couriercorrelations = wanniercourierisometry' * secondrgedcorrelations * wanniercourierisometry
couriercorrelations|>eigspech|>visualize
