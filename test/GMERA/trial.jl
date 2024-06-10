using Zipper
using LinearAlgebra

function checkingdiffwifchosen(truncatedeigenvec, chosenmodes::Subset{Mode},newmodes::Subset{Mode})
    chosencolumns = columns(truncatedeigenvec, FockSpace(chosenmodes))
    U, Σ, Vt = svd(chosencolumns)
    # minsvdvalue::Number = minimum(v for (_, v) in Σ)
    # println("min svdvalue during checking before", minsvdvalue)
    projector = U*U'
    candidate = columns(truncatedeigenvec, FockSpace(newmodes))
    W, Σ, Xt = svd(candidate)
    diff = projector*W-W
    normlist = [norm(columns(diff,FockSpace(m))) for m in diff|>getinspace]
    return sum(normlist)/length(normlist)
end

function choosemodesforwannierizationboth(rankedandgroupedmodes1,rankedandgroupedmodes2,truncatedeigenvec,threshold::Float64,noofmodes::Number)
    chosenmodes = rankedandgroupedmodes[1]
    for modes in rankedandgroupedmodes[2:end]
        if length(chosenmodes)<noofmodes
            if length(chosenmodes)+length(modes)>noofmodes
                @info("more than need")
            else
                diff = checkingdiffwifchosen(truncatedeigenvec, chosenmodes,modes)
                if diff>threshold
                    chosenmodes = chosenmodes+modes
                    chosencolumns = columns(truncatedeigenvec, FockSpace(chosenmodes))
                    # U, Σ, Vt = svd(chosencolumns)
                    # minsvdvalue::Number = minimum(v for (_, v) in Σ)
                    # println("min svdvalue during checking after", minsvdvalue)
                else
                    @info("reject this candidate due to small diff", diff)
                end
            end
        else
            return chosenmodes
        end
    end
    @error("threshold too large")
end

function growchosenmodes(rankedandgroupedmodes,truncatedeigenvec,noofmodes::Number)
    chosenmodes = rankedandgroupedmodes[1]
    remainingmodes = sum(rankedandgroupedmodes)-chosenmodes
    while length(chosenmodes)<noofmodes
        modewifminimumsvdvalsdict = Dict()
        minimumsvdvalslist = []
        for m in remainingmodes
            u,s,vt = svd(columns(truncatedeigenvec,FockSpace(chosenmodes+Subset(m))))
            append!(minimumsvdvalslist,minimum(v for (_, v) in s))
            modewifminimumsvdvalsdict[minimum(v for (_, v) in s)]= m
        end
        chosenmodes = chosenmodes + Subset(modewifminimumsvdvalsdict[maximum(minimumsvdvalslist)])
        remainingmodes = remainingmodes-Subset(modewifminimumsvdvalsdict[maximum(minimumsvdvalslist)])
    end
    return chosenmodes
end

systemsize=32

correlations,H = generatesystem(-0.3 + 0im,0.3 + 0im,-1 + 0im,+(0/(3*sqrt(3)))im,systemsize)
# correlations,H = generatesystem(0im,0im,-1 + 0im,0im,32)

crystalfock = correlations|>getoutspace

scale = Scale([8 0; 0 8], crystalfock|>getcrystal|>getspace)
@info("Performing rgblocking...",scale)
@info("Generating rgblocking transformation...")
blocker = @time scale * crystalfock
@info("Performing rgblocking on correlations...")
blockedcorrelations = @time blocker * correlations * blocker'
blockedcrystal = blockedcorrelations|>getoutspace|>getcrystal

Asublatticeunitcellmodes,Bsublatticeunitcellmodes = generaterefAandBsublatticeunitcellmodes(blockedcrystal)
Asublatticeoffsets = Subset((mode|>getattr(:b)) for mode in Asublatticeunitcellmodes)
Bsublatticeoffsets = Subset((mode|>getattr(:b)) for mode in Bsublatticeunitcellmodes)

blockedcrystalfock = blockedcorrelations|>getoutspace
blockedspace::RealSpace = blockedcrystal|>getspace

refrot = inv([2/3 -1/3; -1/3 2/3]')
c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
c3 = c6^2


@info("Computing local correlations...")
firstcenter = [0,0] ∈ blockedspace
firsthexagonalregion = gethexagonalregion(rot = refrot,crystal=blockedcrystal, center=firstcenter, metricspace=blockedspace)
secondcenter = [2/3,-1/3] ∈ firstrgedspace
secondhexagonalregion = firsthexagonalregion.+secondcenter 
thirdcenter = [1/3,1/3] ∈ secondrgedspace
thirdhexagonalregion = firsthexagonalregion.+thirdcenter 
ref = intersect(intersect(firsthexagonalregion,secondhexagonalregion),thirdhexagonalregion)
visualize(ref)
quantize(ref,1)
localcorrelations = regioncorrelations(blockedcorrelations,quantize((ref+c6*ref+(c6)^2*ref+(c6)^3*ref+(c6)^4*ref+(c6)^5*ref),1))
localspectrum = localcorrelations|>eigspech
localspectrum|>visualize

(ref+c6*ref+(c6)^2*ref+(c6)^3*ref+(c6)^4*ref+(c6)^5*ref)|>visualize

visualize(firsthexagonalregion)
firsthexagonalregionfock = quantize(firsthexagonalregion,1)
    
shiftedfirstcenter1 = [1,0] ∈ blockedspace
# refrot*shiftedfirstcenter1.pos
# shiftedfirsthexagonalregion1 = gethexagonalregion(rot = refrot,crystal=blockedcrystal, center=shiftedfirstcenter1, metricspace=blockedspace)
shiftedfirsthexagonalregion1 = firsthexagonalregion.+shiftedfirstcenter1
# visualize(shiftedfirsthexagonalregion1)
shiftedfirsthexagonalregion1fock = quantize(shiftedfirsthexagonalregion1,1)
    
shiftedfirstcenter2 = [0,1] ∈ blockedspace
shiftedfirsthexagonalregion2 = firsthexagonalregion.+shiftedfirstcenter2
shiftedfirsthexagonalregion2fock = quantize(shiftedfirsthexagonalregion2,1)
    
shiftedfirstcenter3 = [1,1] ∈ blockedspace
shiftedfirsthexagonalregion3 = firsthexagonalregion.+shiftedfirstcenter3
shiftedfirsthexagonalregion3fock = quantize(shiftedfirsthexagonalregion3,1)
    
allregion = firsthexagonalregion+shiftedfirsthexagonalregion1+shiftedfirsthexagonalregion2+shiftedfirsthexagonalregion3
if intersect(allregion,blockedcrystal|>getunitcell) == blockedcrystal|>getunitcell
    @info("cover the whole unitcell")
else
    @error("the allregion cannot cover the whole unitcell ")
end

localcorrelations = regioncorrelations(blockedcorrelations,firsthexagonalregionfock)
localspectrum = localcorrelations|>eigspech
display(localspectrum|>visualize)
filleddict,emptydict = separatefilledandemptymodes(localspectrum)
rankedandgroupedfilledmodes = sortgroupdictwifvaluefilled(filleddict, 0.00001)
rankedandgroupedemptymodes = sortgroupdictwifvalueempty(emptydict, 0.00001)
for data in rankedandgroupedemptymodes
    println(data)
end
for data in rankedandgroupedfilledmodes
    println(data)
end
rankedandgroupedfilledwifemptymodes = [sum(pair) for pair in zip(rankedandgroupedfilledmodes,rankedandgroupedemptymodes)]

rsmodedict = Dict(mode=>norm(euclidean((mode|>getattr(:b))+(mode|>getattr(:r))-firstcenter)) for mode in localspectrum |> geteigenvectors|>getoutspace|>orderedmodes)
sortedrsmode = sortgroupdictwifdist(rsmodedict,false)
frozenrsmodes, courierrsmodes = findfrozenandcourierrsmode(sortedrsmode,138)
courieroffsets = Subset([(mode|>getattr(:b))+(mode|>getattr(:r)) for mode in courierrsmodes])
frozenoffsets = Subset([(mode|>getattr(:b))+(mode|>getattr(:r)) for mode in frozenrsmodes])
visualize(courieroffsets)

frozenseeds = idmap(firsthexagonalregionfock, firsthexagonalregionfock)[:,FockSpace(mode for mode in frozenrsmodes)]
courierseeds = idmap(firsthexagonalregionfock, firsthexagonalregionfock)[:,FockSpace(mode for mode in courierrsmodes)]

truncatedeigenvectofrozen= rows((localspectrum|>geteigenvectors),FockSpace(frozenrsmodes))

# svd(columns(truncatedeigenvectofrozen,FockSpace(sum(rankedandgroupedfilledmodes[1:17])+sum(rankedandgroupedfilledmodes[19:19])+sum(rankedandgroupedemptymodes[1:17])+sum(rankedandgroupedemptymodes[19:19]))))

# checkingdiffwifchosen(truncatedeigenvectofrozen, sum(rankedandgroupedfilledmodes[1:17]),rankedandgroupedemptymodes[20])

chosenfrozenmodes  = growchosenmodes(rankedandgroupedfilledwifemptymodes,truncatedeigenvectofrozen,138)
chosencouriermodes = (localspectrum |> geteigenmodes)-chosenfrozenmodes

localisofrozen = columns(localspectrum |> geteigenvectors, FockSpace(chosenfrozenmodes))
localisocourier = columns(localspectrum |> geteigenvectors, FockSpace(chosencouriermodes))

wannierfrozen = localwannierization(localisofrozen, frozenseeds)
wanniercourier = localwannierization(localisocourier, courierseeds)

(wannierfrozen'*localcorrelations*wannierfrozen)|>eigspec|>visualize|>display

localwannierresults =  Dict(:wannierfrozen => wannierfrozen, :wanniercourier => wanniercourier,
                            :frozenseeds => frozenseeds,:courierseeds => courierseeds)

wannierinfos =  Dict(firsthexagonalregionfock=>localwannierresults)

firstshiftedhexagonalregionfocklist = [(shiftedfirstcenter1,shiftedfirsthexagonalregion1fock), (shiftedfirstcenter2,shiftedfirsthexagonalregion2fock), (shiftedfirstcenter3,shiftedfirsthexagonalregion3fock)]

for (shifted,hexagonalregionfock) in firstshiftedhexagonalregionfocklist
    localcorrelations = regioncorrelations(blockedcorrelations,hexagonalregionfock)
    localspectrum = localcorrelations|>eigspech
    filleddict,emptydict = separatefilledandemptymodes(localspectrum)
    rankedandgroupedfilledmodes =  sortgroupdictwifvaluefilled(filleddict, 0.001)
    rankedandgroupedemptymodes = sortgroupdictwifvalueempty(emptydict, 0.001)
    rankedandgroupedfilledwifemptymodes = [sum(pair) for pair in zip(rankedandgroupedfilledmodes,rankedandgroupedemptymodes)]
    
    rsmodedict = Dict(mode=>norm(euclidean((mode|>getattr(:b))+(mode|>getattr(:r))-(firstcenter+shifted))) for mode in localspectrum |> geteigenvectors|>getoutspace|>orderedmodes)
    sortedrsmode = sortgroupdictwifdist(rsmodedict,false)
    frozenrsmodes, courierrsmodes = findfrozenandcourierrsmode(sortedrsmode,138)
    filledrsmodes,emptyrsmodes = distinguishABsublatticemodesforhexagon(frozenrsmodes,Asublatticeunitcellmodes,Bsublatticeunitcellmodes)

    courierseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in courierrsmodes)]
    frozenseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in frozenrsmodes)]

    truncatedeigenvectofrozen= rows((localspectrum|>geteigenvectors),FockSpace(frozenrsmodes))


    chosenfrozenmodes = growchosenmodes(rankedandgroupedfilledwifemptymodes,truncatedeigenvectofrozen,138)
    chosencouriermodes = (localspectrum |> geteigenmodes)-chosenfrozenmodes

    localisoforzen = columns(localspectrum |> geteigenvectors, FockSpace(chosenfrozenmodes))
    localisocourier = columns(localspectrum |> geteigenvectors, FockSpace(chosencouriermodes))

    wannierfrozen = localwannierization(localisoforzen, frozenseeds)
    wanniercourier = localwannierization(localisocourier, courierseeds)
    
    shiftedlocalwannierresults =  Dict(:wannierfrozen => wannierfrozen, :wanniercourier => wanniercourier,
                                    :frozenseeds => frozenseeds,:courierseeds => courierseeds)

    wannierinfos[hexagonalregionfock] =  shiftedlocalwannierresults
end 

firsthexagonalregionfocklist = [firsthexagonalregionfock, shiftedfirsthexagonalregion1fock, shiftedfirsthexagonalregion2fock, shiftedfirsthexagonalregion3fock]

extendedwannierizedfrozen =  sum(wannierinfos[regionfock][:wannierfrozen] for regionfock in firsthexagonalregionfocklist)
extendedwannierizedcourier =  sum(wannierinfos[regionfock][:wanniercourier] for regionfock in firsthexagonalregionfocklist)
                
origin = [0, 0] ∈ blockedspace
refunictcellfockfrozen = FockSpace(Subset(mode for mode in extendedwannierizedfrozen |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
refunictcellfockcourier = FockSpace(Subset(mode for mode in extendedwannierizedcourier |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))

wannierforzenisometry = globalwannierfunction(blockedcorrelations,extendedwannierizedfrozen[:,refunictcellfockfrozen])
wanniercourierisometry = globalwannierfunction(blockedcorrelations,extendedwannierizedcourier[:,refunictcellfockcourier])

couriercorrelations = wanniercourierisometry' * blockedcorrelations * wanniercourierisometry
couriercorrelations|>eigspech|>visualize|>display
couriercorrelationspectrum = couriercorrelations |> crystalspectrum
purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
purifiedcouriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap

firstrgedcorrelations = purifiedcouriercorrelations
firstrgedcrystalfock = firstrgedcorrelations|>getoutspace
firstrgedcrystal::Crystal = firstrgedcrystalfock|>getcrystal
firstrgedspace::RealSpace = firstrgedcrystal|>getspace
firstrgedcrystal|>getunitcell|>visualize

@info("Computing local correlations...")
secondcenter = [2/3,-1/3] ∈ firstrgedspace
secondhexagonalregion = gethexagonalregion(rot=refrot,crystal=firstrgedcrystal, center=secondcenter, metricspace=firstrgedspace)
secondhexagonalregionfock = quantize(secondhexagonalregion,1)
visualize(secondhexagonalregion )

shiftedsecondcenter1 = [0,1] ∈ firstrgedspace
shiftedsecondhexagonalregion1 = secondhexagonalregion.+shiftedsecondcenter1
shiftedsecondhexagonalregion1fock = quantize(shiftedsecondhexagonalregion1,1)

shiftedsecondcenter2 = [-1,1] ∈ firstrgedspace
shiftedsecondhexagonalregion2 = secondhexagonalregion.+shiftedsecondcenter2
shiftedsecondhexagonalregion2fock = quantize(shiftedsecondhexagonalregion2,1)

allregion2 = secondhexagonalregion+shiftedsecondhexagonalregion1+shiftedsecondhexagonalregion2
# visualize(allregion2)
    if intersect(allregion2,firstrgedcrystal|>getunitcell) == firstrgedcrystal|>getunitcell
        @info("cover the whole unitcell")
    else
        @error("the allregion cannot cover the whole unitcell ")
    end

localcorrelations = regioncorrelations(firstrgedcorrelations,secondhexagonalregionfock)
localspectrum = localcorrelations|>eigspech
display(localspectrum|>visualize)
filleddict,emptydict = separatefilledandemptymodes(localspectrum)
rankedandgroupedfilledmodes = sortgroupdictwifvaluefilled(filleddict, 0.005)
rankedandgroupedemptymodes = sortgroupdictwifvalueempty(emptydict, 0.005)
rankedandgroupedfilledwifemptymodes = [sum(pair) for pair in zip(rankedandgroupedfilledmodes,rankedandgroupedemptymodes)]


rsmodedict = Dict(mode=>norm(euclidean((mode|>getattr(:b))+(mode|>getattr(:r))-secondcenter)) for mode in localspectrum |> geteigenvectors|>getoutspace|>orderedmodes)
sortedrsmode = sortgroupdictwifdist(rsmodedict,false)
frozenrsmodes, courierrsmodes = findfrozenandcourierrsmode(sortedrsmode,126)
courieroffsets = Subset([(mode|>getattr(:b))+(mode|>getattr(:r)) for mode in courierrsmodes])
frozenoffsets = Subset([(mode|>getattr(:b))+(mode|>getattr(:r)) for mode in frozenrsmodes])
visualize(courieroffsets)

frozenseeds = idmap(secondhexagonalregionfock, secondhexagonalregionfock)[:,FockSpace(mode for mode in frozenrsmodes)]
courierseeds = idmap(secondhexagonalregionfock, secondhexagonalregionfock)[:,FockSpace(mode for mode in courierrsmodes)]


truncatedeigenvectofrozen= rows((localspectrum|>geteigenvectors),FockSpace(frozenrsmodes))


chosenfrozenmodes = growchosenmodes(rankedandgroupedfilledwifemptymodes,truncatedeigenvectofrozen,126)
chosencouriermodes = (localspectrum |> geteigenmodes)-chosenfrozenmodes

localisofrozen = columns(localspectrum |> geteigenvectors, FockSpace(chosenfrozenmodes))
localisocourier = columns(localspectrum |> geteigenvectors, FockSpace(chosencouriermodes))

wannierfrozen = localwannierization(localisofrozen, frozenseeds)
wanniercourier = localwannierization(localisocourier, courierseeds)

(wannierfrozen'*localcorrelations*wannierfrozen)|>eigspec|>visualize|>display

localwannierresults =  Dict(:wannierfrozen => wannierfrozen, :wanniercourier => wanniercourier,
                            :frozenseeds => frozenseeds,:courierseeds => courierseeds)

wannierinfos =  Dict(secondhexagonalregionfock=>localwannierresults)

secondshiftedhexagonalregionfocklist = [(shiftedsecondcenter1,shiftedsecondhexagonalregion1fock), (shiftedsecondcenter2,shiftedsecondhexagonalregion2fock)]
                    
    for (shifted,hexagonalregionfock) in secondshiftedhexagonalregionfocklist
        localcorrelations = regioncorrelations(firstrgedcorrelations,hexagonalregionfock)
        localspectrum = localcorrelations|>eigspech
        filleddict,emptydict = separatefilledandemptymodes(localspectrum)
        rankedandgroupedfilledmodes = sortgroupdictwifvaluefilled(filleddict, 0.005)
        rankedandgroupedemptymodes = sortgroupdictwifvalueempty(emptydict, 0.005)
        rankedandgroupedfilledwifemptymodes = [sum(pair) for pair in zip(rankedandgroupedfilledmodes,rankedandgroupedemptymodes)]
                        
        rsmodedict = Dict(mode=>norm(euclidean((mode|>getattr(:b))+(mode|>getattr(:r))-(secondcenter+shifted))) for mode in localspectrum |> geteigenvectors|>getoutspace|>orderedmodes)
        sortedrsmode = sortgroupdictwifdist(rsmodedict,false)
        frozenrsmodes, courierrsmodes = findfrozenandcourierrsmode(sortedrsmode,126)
        
        courierseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in courierrsmodes)]
        frozenseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in frozenrsmodes)]
                    
        truncatedeigenvectofrozen= rows((localspectrum|>geteigenvectors),FockSpace(frozenrsmodes))

        chosenfrozenmodes = growchosenmodes(rankedandgroupedfilledwifemptymodes,truncatedeigenvectofrozen,126)
        chosencouriermodes = (localspectrum |> geteigenmodes)-chosenfrozenmodes
                        
        localisofrozen = columns(localspectrum |> geteigenvectors, FockSpace(chosenfrozenmodes))
        localisocourier = columns(localspectrum |> geteigenvectors, FockSpace(chosencouriermodes))

        wannierfrozen = localwannierization(localisofrozen, frozenseeds)
        wanniercourier = localwannierization(localisocourier, courierseeds)

                            
        shiftedlocalwannierresults =  Dict(:wannierfrozen => wannierfrozen, :wanniercourier => wanniercourier,
                                    :frozenseeds => frozenseeds,:courierseeds => courierseeds)
                    
        wannierinfos[hexagonalregionfock] =  shiftedlocalwannierresults
    end 

secondhexagonalregionfocklist = [secondhexagonalregionfock, shiftedsecondhexagonalregion1fock, shiftedsecondhexagonalregion2fock]
                    
extendedwannierizedfrozen =  sum(wannierinfos[regionfock][:wannierfrozen] for regionfock in secondhexagonalregionfocklist)
extendedwannierizedcourier =  sum(wannierinfos[regionfock][:wanniercourier] for regionfock in secondhexagonalregionfocklist)
                                
origin = [0, 0] ∈ firstrgedspace
refunictcellfockfrozen = FockSpace(Subset(mode for mode in extendedwannierizedfrozen |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
refunictcellfockcourier = FockSpace(Subset(mode for mode in extendedwannierizedcourier |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
                    
wannierforzenisometry = globalwannierfunction(firstrgedcorrelations,extendedwannierizedfrozen[:,refunictcellfockfrozen])
wanniercourierisometry = globalwannierfunction(firstrgedcorrelations,extendedwannierizedcourier[:,refunictcellfockcourier])
                    
frozencorrelations = wannierforzenisometry' * firstrgedcorrelations * wannierforzenisometry              
couriercorrelations2 = wanniercourierisometry' * firstrgedcorrelations * wanniercourierisometry
                    
display(couriercorrelations2|>eigspech|>visualize)
couriercorrelationspectrum2 = couriercorrelations2 |> crystalspectrum
purifiedcorrelationspectrum2 = couriercorrelationspectrum2 |> roundingpurification
purifiedcouriercorrelations2 = purifiedcorrelationspectrum2 |> CrystalFockMap
                    
secondrgedcorrelations = purifiedcouriercorrelations2

secondrgedcrystalfock = secondrgedcorrelations|>getoutspace
secondrgedcrystal::Crystal = secondrgedcrystalfock|>getcrystal
secondrgedspace::RealSpace = secondrgedcrystal|>getspace
    
@info("Computing local correlations...")
thirdcenter = [1/3,1/3] ∈ secondrgedspace
thirdhexagonalregion = gethexagonalregion(rot = refrot,crystal=secondrgedcrystal, center=thirdcenter, metricspace=secondrgedspace)
visualize(thirdhexagonalregion)
thirdhexagonalregionfock = quantize(thirdhexagonalregion,1)

allregion3 = thirdhexagonalregion

    if intersect(allregion3,secondrgedcrystal|>getunitcell) == secondrgedcrystal|>getunitcell
        @info("cover the whole unitcell")
    else
        @error("the allregion cannot cover the whole unitcell ")
    end

localcorrelations = regioncorrelations(secondrgedcorrelations,thirdhexagonalregionfock)
localspectrum = localcorrelations|>eigspech
display(visualize(localspectrum))

filleddict,emptydict = separatefilledandemptymodes(localspectrum)
rankedandgroupedfilledmodes = sortgroupdictwifvaluefilled(filleddict, 0.15)
rankedandgroupedemptymodes = sortgroupdictwifvalueempty(emptydict, 0.15)

rankedandgroupedfilledwifemptymodes = [sum(pair) for pair in zip(rankedandgroupedfilledmodes,rankedandgroupedemptymodes)]

rsmodedict = Dict(mode=>norm(euclidean((mode|>getattr(:b))+(mode|>getattr(:r))-thirdcenter)) for mode in localspectrum |> geteigenvectors|>getoutspace|>orderedmodes)
sortedrsmode = sortgroupdictwifdist(rsmodedict,false)
frozenrsmodes, courierrsmodes = findfrozenandcourierrsmode(sortedrsmode,24)
courieroffsets = Subset([(mode|>getattr(:b))+(mode|>getattr(:r)) for mode in courierrsmodes])
frozenoffsets = Subset([(mode|>getattr(:b))+(mode|>getattr(:r)) for mode in frozenrsmodes])
visualize(courieroffsets)

courierseeds = idmap(thirdhexagonalregionfock, thirdhexagonalregionfock)[:,FockSpace(mode for mode in courierrsmodes)]
frozenseeds = idmap(thirdhexagonalregionfock, thirdhexagonalregionfock)[:,FockSpace(mode for mode in frozenrsmodes)]

truncatedeigenvectofrozen= rows((localspectrum|>geteigenvectors),FockSpace(frozenrsmodes))


chosenfrozenmodes = growchosenmodes(rankedandgroupedfilledwifemptymodes,truncatedeigenvectofrozen,24)
chosencouriermodes = (localspectrum |> geteigenmodes)-chosenfrozenmodes

localisofrozen = columns(localspectrum |> geteigenvectors, FockSpace(chosenfrozenmodes))
localisocourier = columns(localspectrum |> geteigenvectors, FockSpace(chosencouriermodes))

wannierfrozen = localwannierization(localisofrozen, frozenseeds)
wanniercourier = localwannierization(localisocourier, courierseeds)

(wannierfrozen'*localcorrelations*wannierfrozen)|>eigspec|>visualize|>display

localwannierresults =  Dict(:wannierfrozen => wannierfrozen, :wanniercourier => wanniercourier,
                            :frozenseeds => frozenseeds,:courierseeds => courierseeds)

wannierinfos =  Dict(thirdhexagonalregionfock=>localwannierresults)

wannierforzenisometry = globalwannierfunction(secondrgedcorrelations,wannierinfos[thirdhexagonalregionfock][:wannierfrozen])
wanniercourierisometry = globalwannierfunction(secondrgedcorrelations,wannierinfos[thirdhexagonalregionfock][:wanniercourier])

frozencorrelations = wannierforzenisometry' * secondrgedcorrelations * wannierforzenisometry
couriercorrelations3 = wanniercourierisometry' * secondrgedcorrelations * wanniercourierisometry

couriercorrelationspectrum3 = couriercorrelations3 |> crystalspectrum
display(couriercorrelations3|>eigspech|>visualize)
purifiedcorrelationspectrum3 = couriercorrelationspectrum3 |> roundingpurification
purifiedcouriercorrelations3 = purifiedcorrelationspectrum3 |> CrystalFockMap

thirdrgedcorrelations = purifiedcouriercorrelations3
thirdrgedcorrelations|>getoutspace|>getcrystal|>getunitcell|>visualize
recrystalfock = thirdrgedcorrelations|>getoutspace

# rescale = Scale([4 0; 0 4], recrystalfock|>getcrystal|>getspace)
# @info("Performing rgblocking...",rescale)
# @info("Generating rgblocking transformation...")
# reblocker = @time rescale * recrystalfock
# @info("Performing rgblocking on correlations...")
# reblockedcorrelations = @time reblocker * thirdrgedcorrelations * reblocker'
# reblockedcrystal = reblockedcorrelations|>getoutspace|>getcrystal

# reblockedcrystalfock = reblockedcorrelations|>getoutspace
# # blockedcrystal::Crystal = blockedcrystalfock|>getcrystal
# reblockedspace::RealSpace = reblockedcrystal|>getspace

# @info("Computing local correlations...")
# refirstcenter = [0,0] ∈ reblockedspace
# refirsthexagonalregion = gethexagonalregion(crystal=reblockedcrystal, center=refirstcenter, metricspace=reblockedspace)
# refirsthexagonalregionfock = quantize(refirsthexagonalregion,1)
    
# reshiftedfirstcenter1 = [1,0] ∈ reblockedspace
# reshiftedfirsthexagonalregion1 = refirsthexagonalregion.+reshiftedfirstcenter1
# reshiftedfirsthexagonalregion1fock = quantize(reshiftedfirsthexagonalregion1,1)
    
# reshiftedfirstcenter2 = [0,1] ∈ reblockedspace
# reshiftedfirsthexagonalregion2 = refirsthexagonalregion.+reshiftedfirstcenter2
# reshiftedfirsthexagonalregion2fock = quantize(reshiftedfirsthexagonalregion2,1)
    
# reshiftedfirstcenter3 = [1,1] ∈ reblockedspace
# reshiftedfirsthexagonalregion3 = refirsthexagonalregion.+reshiftedfirstcenter3
# reshiftedfirsthexagonalregion3fock = quantize(reshiftedfirsthexagonalregion3,1)
    
# reallregion = refirsthexagonalregion+reshiftedfirsthexagonalregion1+reshiftedfirsthexagonalregion2+reshiftedfirsthexagonalregion3
# reallregion|>visualize
# if intersect(reallregion,reblockedcrystal|>getunitcell) == reblockedcrystal|>getunitcell
#     @info("cover the whole unitcell")
# else
#     @error("the allregion cannot cover the whole unitcell ")
# end

# localcorrelations = regioncorrelations(blockedcorrelations,firsthexagonalregionfock)
# localspectrum = localcorrelations|>eigspech
# display(localspectrum|>visualize)