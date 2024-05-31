using Zipper
using LinearAlgebra
using IterTools

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
crystal = Crystal(unitcell, [32, 32])
reciprocalhashcalibration(crystal.sizes)

modes::Subset{Mode} = quantize(unitcell, 1)|>orderedmodes
m0, m1, m2, m3, m4, m5 = members(modes)
m5|>getattr(:b)

t_a = -0.7 + 0im
t_b = -0.3 + 0im
tₙ = -1 + 0im

onsite = [
    (m0, m0) => t_b,
    (m1, m1) => t_a,
    (m2, m2) => t_b,
    (m3, m3) => t_a,
    (m4, m4) => t_b,
    (m5, m5) => t_a
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

# bonds::FockMap = bondmap([onsite..., nearestneighbor...])
bonds::FockMap = bondmap([nearestneighbor...])

energyspectrum = @time computeenergyspectrum(bonds, crystal=crystal)
groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=3)


groundstateprojector = groundstates|>crystalprojector
correlations = idmap(groundstateprojector|>getoutspace) - groundstateprojector
H = CrystalFockMap(energyspectrum)

function sortgroupdictwifvaluethird(dict::Dict{Mode, Float64},rev::Bool)
    sortedtupledata = sort([(value,key) for (key,value) in dict],rev=rev,by=first)
    # return sortedtupledata
    refvalue = round(sortedtupledata[1][1],digits=3)
    result = []
    subresult = Subset(sortedtupledata[1][2])+Subset(sortedtupledata[2][2])
    for pair in sortedtupledata
        if round(pair[1],digits=3)==refvalue
            subresult = subresult+Subset(pair[2])
        else
            append!(result,tuple([refvalue,subresult]))
            refvalue = round(pair[1],digits=3)
            subresult = Subset(pair[2])
        end
    end
    append!(result,tuple([refvalue,subresult]))

    return result
end

function sortgroupdictwifvaluesec(dict::Dict{Mode, Float64},rev::Bool)
    sortedtupledata = sort([(value,key) for (key,value) in dict],rev=rev,by=first)
    # return sortedtupledata
    refvalue = round(sortedtupledata[1][1],digits=5)
    result = []
    subresult = Subset(sortedtupledata[1][2])+Subset(sortedtupledata[2][2])
    for pair in sortedtupledata
        if round(pair[1],digits=5)==refvalue
            subresult = subresult+Subset(pair[2])
        else
            append!(result,tuple([refvalue,subresult]))
            refvalue = round(pair[1],digits=5)
            subresult = Subset(pair[2])
        end
    end
    append!(result,tuple([refvalue,subresult]))

    return result
end

function sortgroupdictwifvaluefirst(dict::Dict{Mode, Float64},rev::Bool)
    sortedtupledata = sort([(value,key) for (key,value) in dict],rev=rev,by=first)
    # return sortedtupledata
    refvalue = round(sortedtupledata[1][1],digits=10)
    result = []
    subresult = Subset(sortedtupledata[1][2])+Subset(sortedtupledata[2][2])
    for pair in sortedtupledata
        if round(pair[1],digits=10)==refvalue
            subresult = subresult+Subset(pair[2])
        else
            append!(result,tuple([refvalue,subresult]))
            refvalue = round(pair[1],digits=10)
            subresult = Subset(pair[2])
        end
    end
    append!(result,tuple([refvalue,subresult]))

    return result
end

function calculateavgoverlapofmodes(evalwifmodes,spectrum::EigenSpectrum,rsmodes::Subset{Mode})
    overlaps = [(norm(rows(columns(spectrum |> geteigenvectors, FockSpace(mode)),FockSpace(rsmodes))))^2 for mode in evalwifmodes[3]]
    avgoverlap = sum(overlaps)/length(overlaps)
    return tuple(evalwifmodes[3],avgoverlap,evalwifmodes[1],evalwifmodes[2])
end

function separatefilledandemptymodes(spectrum::EigenSpectrum)
    evaldict = spectrum |> geteigenvalues
    filleddict = Dict( key=>value for (key,value) in evaldict if value<0.5)
    emptydict = Dict( key=>value for (key,value) in evaldict if value>0.5)
    return filleddict,emptydict
end

@info("Starting RG...")
    crystalfock = correlations|>getoutspace

    scale = Scale([4 0; 0 4], crystalfock|>getcrystal|>getspace)
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
    visualize(regioncorrelations(blockedcorrelations,firsthexagonalregionfock)|>eigspech)

    shiftedfirstcenter1 = [1,0] ∈ blockedspace
    shiftedfirsthexagonalregion1 = firsthexagonalregion.+shiftedfirstcenter1
    shiftedfirsthexagonalregion1fock = quantize(shiftedfirsthexagonalregion1,1)

    shiftedfirstcenter2 = [0,1] ∈ blockedspace
    shiftedfirsthexagonalregion2 = firsthexagonalregion.+shiftedfirstcenter2
    shiftedfirsthexagonalregion2fock = quantize(shiftedfirsthexagonalregion2,1)

    shiftedfirstcenter3 = [1,1] ∈ blockedspace
    shiftedfirsthexagonalregion3 = firsthexagonalregion.+shiftedfirstcenter3
    shiftedfirsthexagonalregion3fock = quantize(shiftedfirsthexagonalregion3,1)

    allregion = firsthexagonalregion+shiftedfirsthexagonalregion1+shiftedfirsthexagonalregion2+shiftedfirsthexagonalregion3
    visualize(allregion)
    if intersect(allregion,blockedcrystal|>getunitcell) == blockedcrystal|>getunitcell
        @info("cover the whole unitcell")
    else
        @error("the allregion cannot cover the whole unitcell ")
    end


    localcorrelations = regioncorrelations(blockedcorrelations,firsthexagonalregionfock)
    # (regioncorrelations(blockedcorrelations,firsthexagonalregionfock)|> rep |> Matrix |> Hermitian |> eigvals)[1:32]
    localspectrum = localcorrelations|>eigspech
    filleddict,emptydict = separatefilledandemptymodes(localspectrum)
    
    sortedfilled = sortgroupdictwifvaluefirst(filleddict,false)
    sortedempty = sortgroupdictwifvaluefirst(emptydict,true)
    rsmodedict = Dict(mode=>norm(euclidean((mode|>getattr(:b))+(mode|>getattr(:r))-firstcenter)) for mode in localspectrum |> geteigenvectors|>getoutspace|>orderedmodes)
    sortedrsmode = sortgroupdictwifvalue(rsmodedict,false)
    frozenrsmodes = sum([pair[2] for pair in sortedrsmode[1:4]])
    frozenrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in frozenrsmodes)
    visualize(frozenrsmodesoffset)
    courierrsmodes = sum([pair[2] for pair in sortedrsmode[5:end]])
    courierrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in courierrsmodes)
    visualize(courierrsmodesoffset)

    sortedwifevalfrozen = [(pair[1][1],pair[2][1],pair[1][2]+pair[2][2]) for pair in zip(sortedfilled,sortedempty)]
    frozenmodeswifoverlap = [calculateavgoverlapofmodes(frozen,localspectrum,frozenrsmodes) for frozen in sortedwifevalfrozen]
    sortedfrozenmodeswifoverlap = sort(frozenmodeswifoverlap,by=x->x[2],rev=true)
    sortedfrozenmodeswifoverlap[1:8]
    # filledmodeswifoverlap = [calculateavgoverlapofmodes(filled,localspectrum,frozenrsmodes) for filled in sortedfilled]
    # emptymodeswifoverlap = [calculateavgoverlapofmodes(empty,localspectrum,frozenrsmodes) for empty in sortedempty]
    # sortedfilledmodeswifoverlap = sort(filledmodeswifoverlap,by=x->x[2],rev=true)
    # sortedfilledmodeswifoverlap
    # sortedemptymodeswifoverlap = sort(emptymodeswifoverlap,by=x->x[2],rev=true)
    # sortedemptymodeswifoverlap
    # rankedfilledandemptymodes = [pair[1][1]+pair[2][1] for pair in zip(sortedfilledmodeswifoverlap,sortedemptymodeswifoverlap) if length(pair[1][1]) == length(pair[2][1])]
    chosenfrozenmodes = sum([data[1] for data in sortedfrozenmodeswifoverlap[1:8]])
    chosencouriermodes = (localspectrum |> geteigenmodes)-chosenfrozenmodes
    

    
    courierseeds = idmap(firsthexagonalregionfock, firsthexagonalregionfock)[:,FockSpace(mode for mode in courierrsmodes)]
    localisocourier = columns(localspectrum |> geteigenvectors, FockSpace(chosencouriermodes))
    wanniercourier = localwannierization(localisocourier, courierseeds)

    frozenseeds = idmap(firsthexagonalregionfock, firsthexagonalregionfock)[:,FockSpace(mode for mode in frozenrsmodes)]
    localisoforzen = columns(localspectrum |> geteigenvectors, FockSpace(chosenfrozenmodes))
    wannierfrozen = localwannierization(localisoforzen, frozenseeds)


    localwannierresults =  Dict(:wannierfrozen => wannierfrozen, :wanniercourier => wanniercourier,
                    :frozenseeds => frozenseeds,:courierseeds => courierseeds)

    wannierinfos =  Dict(firsthexagonalregionfock=>localwannierresults)

    firstshiftedhexagonalregionfocklist = [(shiftedfirstcenter1,shiftedfirsthexagonalregion1fock), (shiftedfirstcenter2,shiftedfirsthexagonalregion2fock), (shiftedfirstcenter3,shiftedfirsthexagonalregion3fock)]

    for (shifted,hexagonalregionfock) in firstshiftedhexagonalregionfocklist
        localcorrelations = regioncorrelations(blockedcorrelations,hexagonalregionfock)
        localspectrum = localcorrelations|>eigspech
        filleddict,emptydict = separatefilledandemptymodes(localspectrum)
        
        sortedfilled = sortgroupdictwifvaluefirst(filleddict,false)
        sortedempty = sortgroupdictwifvaluefirst(emptydict,true)
        rsmodedict = Dict(mode=>norm(euclidean((mode|>getattr(:b))+(mode|>getattr(:r))-(firstcenter+shifted))) for mode in localspectrum |> geteigenvectors|>getoutspace|>orderedmodes)
        sortedrsmode = sortgroupdictwifvalue(rsmodedict,false)
        frozenrsmodes = sum([pair[2] for pair in sortedrsmode[1:4]])
        frozenrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in frozenrsmodes)
        visualize(frozenrsmodesoffset)
        courierrsmodes = sum([pair[2] for pair in sortedrsmode[5:end]])
        courierrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in courierrsmodes)
        visualize(courierrsmodesoffset)

        sortedwifevalfrozen = [(pair[1][1],pair[2][1],pair[1][2]+pair[2][2]) for pair in zip(sortedfilled,sortedempty)]
        frozenmodeswifoverlap = [calculateavgoverlapofmodes(frozen,localspectrum,frozenrsmodes) for frozen in sortedwifevalfrozen]
        sortedfrozenmodeswifoverlap = sort(frozenmodeswifoverlap,by=x->x[2],rev=true)
        sortedfrozenmodeswifoverlap[1:8]
        # filledmodeswifoverlap = [calculateavgoverlapofmodes(filled,localspectrum,frozenrsmodes) for filled in sortedfilled]
        # emptymodeswifoverlap = [calculateavgoverlapofmodes(empty,localspectrum,frozenrsmodes) for empty in sortedempty]
        # sortedfilledmodeswifoverlap = sort(filledmodeswifoverlap,by=x->x[2],rev=true)
        # sortedfilledmodeswifoverlap
        # sortedemptymodeswifoverlap = sort(emptymodeswifoverlap,by=x->x[2],rev=true)
        # sortedemptymodeswifoverlap
        # rankedfilledandemptymodes = [pair[1][1]+pair[2][1] for pair in zip(sortedfilledmodeswifoverlap,sortedemptymodeswifoverlap) if length(pair[1][1]) == length(pair[2][1])]
        chosenfrozenmodes = sum([data[1] for data in sortedfrozenmodeswifoverlap[1:8]])
        chosencouriermodes = (localspectrum |> geteigenmodes)-chosenfrozenmodes
        

        
        courierseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in courierrsmodes)]
        localisocourier = columns(localspectrum |> geteigenvectors, FockSpace(chosencouriermodes))
        wanniercourier = localwannierization(localisocourier, courierseeds)

        frozenseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in frozenrsmodes)]
        localisoforzen = columns(localspectrum |> geteigenvectors, FockSpace(chosenfrozenmodes))
        wannierfrozen = localwannierization(localisoforzen, frozenseeds)

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
    couriercorrelations|>eigspech|>visualize
    couriercorrelationspectrum = couriercorrelations |> crystalspectrum
    purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
    purifiedcouriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap

    frozencorrelations = wannierforzenisometry' * blockedcorrelations * wannierforzenisometry
    frozencorrelations|>eigspech|>visualize

    visualize(purifiedcouriercorrelations|>getinspace|>getcrystal|>getunitcell)

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

    
    sortedfilled = sortgroupdictwifvaluesec(filleddict,false)
    sortedempty = sortgroupdictwifvaluesec(emptydict,true)
    for data in sortedfilled
        println(data)
    end
    for data in sortedempty
        println(data)
    end
    rsmodedict = Dict(mode=>norm(euclidean((mode|>getattr(:b))+(mode|>getattr(:r))-secondcenter)) for mode in localspectrum |> geteigenvectors|>getoutspace|>orderedmodes)
    sortedrsmode = sortgroupdictwifvalue(rsmodedict,false)
    sortedrsmode
    frozenrsmodes = sum([pair[2] for pair in sortedrsmode[1:4]])
    frozenrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in frozenrsmodes)
    visualize(frozenrsmodesoffset)
    courierrsmodes = sum([pair[2] for pair in sortedrsmode[5:end]])
    courierrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in courierrsmodes)
    visualize(courierrsmodesoffset)

    sortedwifevalfrozen = [(pair[1][1],pair[2][1],pair[1][2]+pair[2][2]) for pair in zip(sortedfilled,sortedempty)]
    frozenmodeswifoverlap = [calculateavgoverlapofmodes(frozen,localspectrum,frozenrsmodes) for frozen in sortedwifevalfrozen]
    sortedfrozenmodeswifoverlap = sort(frozenmodeswifoverlap,by=x->x[2],rev=true)
    sortedfrozenmodeswifoverlap[1:8]
    # filledmodeswifoverlap = [calculateavgoverlapofmodes(filled,localspectrum,frozenrsmodes) for filled in sortedfilled]
    # emptymodeswifoverlap = [calculateavgoverlapofmodes(empty,localspectrum,frozenrsmodes) for empty in sortedempty]
    # sortedfilledmodeswifoverlap = sort(filledmodeswifoverlap,by=x->x[2],rev=true)
    # sortedfilledmodeswifoverlap
    # sortedemptymodeswifoverlap = sort(emptymodeswifoverlap,by=x->x[2],rev=true)
    # sortedemptymodeswifoverlap
    # rankedfilledandemptymodes = [pair[1][1]+pair[2][1] for pair in zip(sortedfilledmodeswifoverlap,sortedemptymodeswifoverlap) if length(pair[1][1]) == length(pair[2][1])]
    chosenfrozenmodes = sum([data[1] for data in sortedfrozenmodeswifoverlap[1:10]])
    chosencouriermodes = (localspectrum |> geteigenmodes)-chosenfrozenmodes

    courierseeds = idmap(secondhexagonalregionfock, secondhexagonalregionfock)[:,FockSpace(mode for mode in courierrsmodes)]
    localisocourier = columns(localspectrum |> geteigenvectors, FockSpace(chosencouriermodes))
    wanniercourier = localwannierization(localisocourier, courierseeds)

    visualize(wanniercourier'*localcorrelations*wanniercourier|>eigspech)

    frozenseeds = idmap(secondhexagonalregionfock, secondhexagonalregionfock)[:,FockSpace(mode for mode in frozenrsmodes)]
    localisoforzen = columns(localspectrum |> geteigenvectors, FockSpace(chosenfrozenmodes))
    wannierfrozen = localwannierization(localisoforzen, frozenseeds)

    visualize(wannierfrozen'*localcorrelations*wannierfrozen|>eigspech)

    localwannierresults =  Dict(:wannierfrozen => wannierfrozen, :wanniercourier => wanniercourier,
                    :frozenseeds => frozenseeds,:courierseeds => courierseeds)

    wannierinfos =  Dict(secondhexagonalregionfock=>localwannierresults)

    secondshiftedhexagonalregionfocklist = [(shiftedsecondcenter1,shiftedsecondhexagonalregion1fock), (shiftedsecondcenter2,shiftedsecondhexagonalregion2fock)]


    for (shifted,hexagonalregionfock) in secondshiftedhexagonalregionfocklist
        localcorrelations = regioncorrelations(firstrgedcorrelations,hexagonalregionfock)
        localspectrum = localcorrelations|>eigspech
        filleddict,emptydict = separatefilledandemptymodes(localspectrum)
        
        sortedfilled = sortgroupdictwifvaluesec(filleddict,false)
        sortedempty = sortgroupdictwifvaluesec(emptydict,true)
        rsmodedict = Dict(mode=>norm(euclidean((mode|>getattr(:b))+(mode|>getattr(:r))-(secondcenter+shifted))) for mode in localspectrum |> geteigenvectors|>getoutspace|>orderedmodes)
        sortedrsmode = sortgroupdictwifvalue(rsmodedict,false)
        frozenrsmodes = sum([pair[2] for pair in sortedrsmode[1:4]])
        frozenrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in frozenrsmodes)
        visualize(frozenrsmodesoffset)
        courierrsmodes = sum([pair[2] for pair in sortedrsmode[5:end]])
        courierrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in courierrsmodes)
        visualize(courierrsmodesoffset)

        sortedwifevalfrozen = [(pair[1][1],pair[2][1],pair[1][2]+pair[2][2]) for pair in zip(sortedfilled,sortedempty)]
        frozenmodeswifoverlap = [calculateavgoverlapofmodes(frozen,localspectrum,frozenrsmodes) for frozen in sortedwifevalfrozen]
        sortedfrozenmodeswifoverlap = sort(frozenmodeswifoverlap,by=x->x[2],rev=true)
        sortedfrozenmodeswifoverlap[1:8]
        # filledmodeswifoverlap = [calculateavgoverlapofmodes(filled,localspectrum,frozenrsmodes) for filled in sortedfilled]
        # emptymodeswifoverlap = [calculateavgoverlapofmodes(empty,localspectrum,frozenrsmodes) for empty in sortedempty]
        # sortedfilledmodeswifoverlap = sort(filledmodeswifoverlap,by=x->x[2],rev=true)
        # sortedfilledmodeswifoverlap
        # sortedemptymodeswifoverlap = sort(emptymodeswifoverlap,by=x->x[2],rev=true)
        # sortedemptymodeswifoverlap
        # rankedfilledandemptymodes = [pair[1][1]+pair[2][1] for pair in zip(sortedfilledmodeswifoverlap,sortedemptymodeswifoverlap) if length(pair[1][1]) == length(pair[2][1])]
        chosenfrozenmodes = sum([data[1] for data in sortedfrozenmodeswifoverlap[1:10]])
        chosencouriermodes = (localspectrum |> geteigenmodes)-chosenfrozenmodes

        courierseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in courierrsmodes)]
        localisocourier = columns(localspectrum |> geteigenvectors, FockSpace(chosencouriermodes))
        wanniercourier = localwannierization(localisocourier, courierseeds)

        display(visualize(wanniercourier'*localcorrelations*wanniercourier|>eigspech))

        frozenseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in frozenrsmodes)]
        localisoforzen = columns(localspectrum |> geteigenvectors, FockSpace(chosenfrozenmodes))
        wannierfrozen = localwannierization(localisoforzen, frozenseeds)

        display(visualize(wannierfrozen'*localcorrelations*wannierfrozen|>eigspech))

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

    # shiftedthirdcenter1 = [0,1] ∈ secondrgedspace
    # shiftedthirdhexagonalregion1 = thirdhexagonalregion.+shiftedthirdcenter1
    # shiftedthirdhexagonalregion1fock = quantize(shiftedthirdhexagonalregion1,1)

    # shiftedthirdcenter2 = [1,0] ∈ secondrgedspace
    # shiftedthirdhexagonalregion2 = thirdhexagonalregion.+shiftedthirdcenter2
    # shiftedthirdhexagonalregion2fock = quantize(shiftedthirdhexagonalregion2,1)

    allregion3 = thirdhexagonalregion
    visualize(allregion3)
    visualize(secondrgedcrystal|>getunitcell)
    visualize(intersect(allregion3,firstrgedcrystal|>getunitcell))
    if intersect(allregion3,firstrgedcrystal|>getunitcell) == secondrgedcrystal|>getunitcell
        @info("cover the whole unitcell")
    else
        @error("the allregion cannot cover the whole unitcell ")
    end

    localcorrelations = regioncorrelations(secondrgedcorrelations,thirdhexagonalregionfock)
    localspectrum = localcorrelations|>eigspech
    display(visualize(localspectrum))

    filleddict,emptydict = separatefilledandemptymodes(localspectrum)

    
    sortedfilled = sortgroupdictwifvaluethird(filleddict,false)
    sortedempty = sortgroupdictwifvaluethird(emptydict,true)
    for data in sortedempty
        println(data)
    end
    for data in sortedfilled
        println(data)
    end
    rsmodedict = Dict(mode=>norm(euclidean((mode|>getattr(:b))+(mode|>getattr(:r))-thirdcenter)) for mode in localspectrum |> geteigenvectors|>getoutspace|>orderedmodes)
    sortedrsmode = sortgroupdictwifvalue(rsmodedict,false)
    sortedrsmode
    frozenrsmodes = sum([pair[2] for pair in sortedrsmode[1:3]])
    frozenrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in frozenrsmodes)
    visualize(frozenrsmodesoffset)
    courierrsmodes = sum([pair[2] for pair in sortedrsmode[4:end]])
    courierrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in courierrsmodes)
    visualize(courierrsmodesoffset)

    sortedwifevalfrozen = [(pair[1][1],pair[2][1],pair[1][2]+pair[2][2]) for pair in zip(sortedfilled,sortedempty)]
    frozenmodeswifoverlap = [calculateavgoverlapofmodes(frozen,localspectrum,frozenrsmodes) for frozen in sortedwifevalfrozen]
    sortedfrozenmodeswifoverlap = sort(frozenmodeswifoverlap,by=x->x[2],rev=true)
    sortedfrozenmodeswifoverlap[1:3]
    # filledmodeswifoverlap = [calculateavgoverlapofmodes(filled,localspectrum,frozenrsmodes) for filled in sortedfilled]
    # emptymodeswifoverlap = [calculateavgoverlapofmodes(empty,localspectrum,frozenrsmodes) for empty in sortedempty]
    # sortedfilledmodeswifoverlap = sort(filledmodeswifoverlap,by=x->x[2],rev=true)
    # sortedfilledmodeswifoverlap
    # sortedemptymodeswifoverlap = sort(emptymodeswifoverlap,by=x->x[2],rev=true)
    # sortedemptymodeswifoverlap
    # rankedfilledandemptymodes = [pair[1][1]+pair[2][1] for pair in zip(sortedfilledmodeswifoverlap,sortedemptymodeswifoverlap) if length(pair[1][1]) == length(pair[2][1])]
    chosenfrozenmodes = sum([data[1] for data in sortedfrozenmodeswifoverlap[1:5]])
    chosencouriermodes = (localspectrum |> geteigenmodes)-chosenfrozenmodes

    courierseeds = idmap(thirdhexagonalregionfock, thirdhexagonalregionfock)[:,FockSpace(mode for mode in courierrsmodes)]
    localisocourier = columns(localspectrum |> geteigenvectors, FockSpace(chosencouriermodes))
    wanniercourier = localwannierization(localisocourier, courierseeds)

    visualize(wanniercourier'*localcorrelations*wanniercourier|>eigspech)

    frozenseeds = idmap(thirdhexagonalregionfock, thirdhexagonalregionfock)[:,FockSpace(mode for mode in frozenrsmodes)]
    localisoforzen = columns(localspectrum |> geteigenvectors, FockSpace(chosenfrozenmodes))
    wannierfrozen = localwannierization(localisoforzen, frozenseeds)

    visualize(wannierfrozen'*localcorrelations*wannierfrozen|>eigspech)

    localwannierresults =  Dict(:wannierfrozen => wannierfrozen, :wanniercourier => wanniercourier,
                    :frozenseeds => frozenseeds,:courierseeds => courierseeds)

    wannierinfos =  Dict(thirdhexagonalregionfock=>localwannierresults)

            
    # origin = [0, 0] ∈ firstrgedspace
    # refunictcellfockfrozen = FockSpace(Subset(mode for mode in wannierinfos[thirdhexagonalregionfock][:wannierfrozen] |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    # refunictcellfockcourier = FockSpace(Subset(mode for mode in wannierinfos[thirdhexagonalregionfock][:wanniercourier] |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))


    wannierforzenisometry = globalwannierfunction(secondrgedcorrelations,wannierinfos[thirdhexagonalregionfock][:wannierfrozen])
    wanniercourierisometry = globalwannierfunction(secondrgedcorrelations,wannierinfos[thirdhexagonalregionfock][:wanniercourier])

    frozencorrelations = wannierforzenisometry' * secondrgedcorrelations * wannierforzenisometry
    frozencorrelations|>eigspech|>visualize

    couriercorrelations = wanniercourierisometry' * secondrgedcorrelations * wanniercourierisometry
    couriercorrelations|>eigspech|>visualize

    visualize(couriercorrelations|>getinspace|>getcrystal|>getunitcell)

