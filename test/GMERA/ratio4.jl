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

# function findseedsbyweightanddist(sortedrsmodeslist)
#     refweight = sortedrsmodeslist[1][1]
#     refdist = sortedrsmodeslist[1][2]
#     modelist = Subset(data[3] for data in sortedrsmodeslist if (data[1]==refweight) & (data[2]==refdist))
#     return tuple(refdist,modelist)
# end

# function generatesortedlistfromspectrum(localspectrum::EigenSpectrum)
#     modelist = [mode for mode in localspectrum|> geteigenmodes]
#     realspacemode = localspectrum |> geteigenvectors|>getoutspace|>orderedmodes
#     modewifsortedrsmodeslist = sort([(findseedsbyweightanddist(sort([(round(norm((localspectrum |> geteigenvectors)[rm,mode])^2,digits=3),round(norm(euclidean((rm|>getattr(:b))+(rm|>getattr(:r)))),digits=3),rm) for rm in realspacemode],by=first,rev=true))[1],mode,(localspectrum|>geteigenvalues)[mode],findseedsbyweightanddist(sort([(round(norm((localspectrum |> geteigenvectors)[rm,mode])^2,digits=3),round(norm(euclidean((rm|>getattr(:b))+(rm|>getattr(:r)))),digits=3),rm) for rm in realspacemode],by=first,rev=true))[2]) for mode in modelist],by=first)
#     return modewifsortedrsmodeslist
# end

function sortgroupdictwifvalue(dict::Dict{Mode, Float64},rev::Bool)
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
    overlaps = [(norm(rows(columns(spectrum |> geteigenvectors, FockSpace(mode)),FockSpace(rsmodes))))^2 for mode in evalwifmodes[2]]
    avgoverlap = sum(overlaps)/length(overlaps)
    return tuple(evalwifmodes[2],avgoverlap,evalwifmodes[1])
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
    (regioncorrelations(blockedcorrelations,firsthexagonalregionfock)|> rep |> Matrix |> Hermitian |> eigvals)[1:32]
    localspectrum = localcorrelations|>eigspech
    filleddict,emptydict = separatefilledandemptymodes(localspectrum)
    
    sortedfilled = sortgroupdictwifvalue(filleddict,false)
    sortedempty = sortgroupdictwifvalue(emptydict,true)
    rsmodedict = Dict(mode=>norm(euclidean((mode|>getattr(:b))+(mode|>getattr(:r)))) for mode in localspectrum |> geteigenvectors|>getoutspace|>orderedmodes)
    sortedrsmode = sortedgroupdictwifvalue(rsmodedict,false)
    frozenrsmodes = sum([pair[2] for pair in sortedrsmode[1:4]])
    frozenrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in frozenmodes)
    visualize(frozenrsmodesoffset)
    courierrsmodes = sum([pair[2] for pair in sortedrsmode[5:end]])
    courierrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in courierrsmodes)
    visualize(couriermodesoffset)

    filledmodeswifoverlap = [calculateavgoverlapofmodes(filled,localspectrum,frozenmodes) for filled in sortedfilled]
    emptymodeswifoverlap = [calculateavgoverlapofmodes(empty,localspectrum,frozenmodes) for empty in sortedempty]
    sortedfilledmodeswifoverlap = sort(filledmodeswifoverlap,by=x->x[2],rev=true)
    sortedfilledmodeswifoverlap
    sortedemptymodeswifoverlap = sort(emptymodeswifoverlap,by=x->x[2],rev=true)
    sortedemptymodeswifoverlap
    rankedfilledandemptymodes = [pair[1][1]+pair[2][1] for pair in zip(sortedfilledmodeswifoverlap,sortedemptymodeswifoverlap) if length(pair[1][1]) == length(pair[2][1])]
    # for pair in zip(sortedfilledmodeswifoverlap,sortedemptymodeswifoverlap)
    #     if length(pair[1][1]) == length(pair[2][1])
    #         append!(rankedfilledandemptymodes,tuple(pair[1][1],pair[2][1]))
    #     else
    #         @error("number of filled and empty modes not matched")
    #     end
    # end
    chosenfrozenmodes = sum(rankedfilledandemptymodes[1:8])
    chosencouriermodes = (localspectrum |> geteigenmodes)-chosenfrozenmodes
    

    
    courierseeds = idmap(firsthexagonalregionfock, firsthexagonalregionfock)[:,FockSpace(mode for mode in courierrsmodes)]
    localisocourier = columns(localspectrum |> geteigenvectors, FockSpace(chosencouriermodes))
    wanniercourier = localwannierization(localisocourier, courierseeds)

    frozenseeds = idmap(firsthexagonalregionfock, firsthexagonalregionfock)[:,FockSpace(mode for mode in frozenrsmodes)]
    localisoforzen = columns(localspectrum |> geteigenvectors, FockSpace(chosenfrozenmodes))
    wannierfrozen = localwannierization(localisoforzen, frozenseeds)

    
    display(visualize(localspectrum))
    localiso = localcorrelations|>modeselectionbycount(18)

    localwannierseeds = localwannierseedslists(localiso)

    emptyseeds = idmap(firsthexagonalregionfock, firsthexagonalregionfock)[:,FockSpace(mode for mode in localwannierseeds[:empty])] 
    filledseeds = idmap(firsthexagonalregionfock, firsthexagonalregionfock)[:,FockSpace(mode for mode in localwannierseeds[:filled])] 
    courierseeds = idmap(firsthexagonalregionfock, firsthexagonalregionfock)[:,FockSpace(mode for mode in localwannierseeds[:courier])]

    wannierempty = localwannierization(localiso[:empty], emptyseeds)
    wannierfilled = localwannierization(localiso[:filled], filledseeds)
    wanniercourier = localwannierization(localiso[:courier], courierseeds)

    localwannierresults =  Dict(:wannierempty => wannierempty, :wannierfilled => wannierfilled, :wanniercourier => wanniercourier,
                    :emptyseeds => emptyseeds, :filledseeds => filledseeds, :courierseeds => courierseeds,
                    :bempty=>localwannierseeds[:bempty],:bfilled=>localwannierseeds[:bfilled],:bcourier=>localwannierseeds[:bcourier])

    wannierinfos =  Dict(firsthexagonalregionfock=>localwannierresults)

    firstshiftedhexagonalregionfocklist = [shiftedfirsthexagonalregion1fock, shiftedfirsthexagonalregion2fock, shiftedfirsthexagonalregion3fock]

    for hexagonalregionfock in firstshiftedhexagonalregionfocklist
        localcorrelations = regioncorrelations(blockedcorrelations,hexagonalregionfock)
        localspectrum = localcorrelations|>eigspech
        localspectrum|>visualize
        localiso = localcorrelations|>modeselectionbycount(18)

        localwannierresults = localwannierseedslists(localiso)

        emptyseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in localwannierresults[:empty])] 
        filledseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in localwannierresults[:filled])] 
        courierseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in localwannierresults[:courier])]

        wannierempty = localwannierization(localiso[:empty], emptyseeds)
        wannierfilled = localwannierization(localiso[:filled], filledseeds)
        wanniercourier = localwannierization(localiso[:courier], courierseeds)

        shiftedlocalwannierresults =  Dict(:wannierempty => wannierempty, :wannierfilled => wannierfilled, :wanniercourier => wanniercourier,
                    :emptyseeds => emptyseeds, :filledseeds => filledseeds, :courierseeds => courierseeds,
                    :bempty=>localwannierresults[:bempty],:bfilled=>localwannierresults[:bfilled],:bcourier=>localwannierresults[:bcourier])

        wannierinfos[hexagonalregionfock] = shiftedlocalwannierresults
    end 

    firsthexagonalregionfocklist = [firsthexagonalregionfock, shiftedfirsthexagonalregion1fock, shiftedfirsthexagonalregion2fock, shiftedfirsthexagonalregion3fock]

    extendedwannierizedempty =  sum(wannierinfos[regionfock][:wannierempty] for regionfock in firsthexagonalregionfocklist)
    extendedwannierizedfilled =  sum(wannierinfos[regionfock][:wannierfilled] for regionfock in firsthexagonalregionfocklist)
    extendedwannierizedcourier =  sum(wannierinfos[regionfock][:wanniercourier] for regionfock in firsthexagonalregionfocklist)
            
    origin = [0, 0] ∈ blockedspace
    refunictcellfockempty = FockSpace(Subset(mode for mode in extendedwannierizedempty |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    refunictcellfockfilled = FockSpace(Subset(mode for mode in extendedwannierizedfilled |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    refunictcellfockcourier = FockSpace(Subset(mode for mode in extendedwannierizedcourier |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))


    wannieremptyisometry = globalwannierfunction(blockedcorrelations,extendedwannierizedempty[:,refunictcellfockempty])
    wannierfilledisometry = globalwannierfunction(blockedcorrelations,extendedwannierizedfilled[:,refunictcellfockfilled])
    wanniercourierisometry = globalwannierfunction(blockedcorrelations,extendedwannierizedcourier[:,refunictcellfockcourier])

    couriercorrelations = wanniercourierisometry' * blockedcorrelations * wanniercourierisometry
    couriercorrelationspectrum = couriercorrelations |> crystalspectrum
    purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
    purifiedcouriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap

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

    localcorrelations = regioncorrelations(purifiedcouriercorrelations,secondhexagonalregionfock)
    localspectrum = localcorrelations|>eigspech
    display(visualize(localspectrum))
    localiso = localcorrelations|>modeselectionbycount(10)

    localwannierseeds = localwannierseedslists(localiso)

    emptyseeds = idmap(secondhexagonalregionfock, secondhexagonalregionfock)[:,FockSpace(mode for mode in localwannierseeds[:empty])] 
    filledseeds = idmap(secondhexagonalregionfock, secondhexagonalregionfock)[:,FockSpace(mode for mode in localwannierseeds[:filled])] 
    courierseeds = idmap(secondhexagonalregionfock, secondhexagonalregionfock)[:,FockSpace(mode for mode in localwannierseeds[:courier])]

    wannierempty = localwannierization(localiso[:empty], emptyseeds)
    wannierfilled = localwannierization(localiso[:filled], filledseeds)
    wanniercourier = localwannierization(localiso[:courier], courierseeds)

    localwannierresults =  Dict(:wannierempty => wannierempty, :wannierfilled => wannierfilled, :wanniercourier => wanniercourier,
                    :emptyseeds => emptyseeds, :filledseeds => filledseeds, :courierseeds => courierseeds,
                    :bempty=>localwannierseeds[:bempty],:bfilled=>localwannierseeds[:bfilled],:bcourier=>localwannierseeds[:bcourier])

    wannierinfos =  Dict(secondhexagonalregionfock =>localwannierresults)

    secondshiftedhexagonalregionfocklist = [shiftedsecondhexagonalregion1fock, shiftedsecondhexagonalregion2fock]

    for hexagonalregionfock in secondshiftedhexagonalregionfocklist
        localcorrelations = regioncorrelations(purifiedcouriercorrelations,hexagonalregionfock)
        localspectrum = localcorrelations|>eigspech
        localspectrum|>visualize
        localiso = localcorrelations|>modeselectionbycount(10)

        localwannierresults = localwannierseedslists(localiso)

        emptyseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in localwannierresults[:empty])] 
        filledseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in localwannierresults[:filled])] 
        courierseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in localwannierresults[:courier])]

        wannierempty = localwannierization(localiso[:empty], emptyseeds)
        wannierfilled = localwannierization(localiso[:filled], filledseeds)
        wanniercourier = localwannierization(localiso[:courier], courierseeds)

        shiftedlocalwannierresults =  Dict(:wannierempty => wannierempty, :wannierfilled => wannierfilled, :wanniercourier => wanniercourier,
                    :emptyseeds => emptyseeds, :filledseeds => filledseeds, :courierseeds => courierseeds,
                    :bempty=>localwannierseeds[:bempty],:bfilled=>localwannierseeds[:bfilled],:bcourier=>localwannierseeds[:bcourier])

        wannierinfos[hexagonalregionfock] = shiftedlocalwannierresults
    end 

    secondhexagonalregionfocklist = [secondhexagonalregionfock, shiftedsecondhexagonalregion1fock, shiftedsecondhexagonalregion2fock]

    extendedwannierizedempty =  sum(wannierinfos[regionfock][:wannierempty] for regionfock in secondhexagonalregionfocklist)
    extendedwannierizedfilled =  sum(wannierinfos[regionfock][:wannierfilled] for regionfock in secondhexagonalregionfocklist)
    extendedwannierizedcourier =  sum(wannierinfos[regionfock][:wanniercourier] for regionfock in secondhexagonalregionfocklist)
            
    origin = [0, 0] ∈ firstrgedspace
    refunictcellfockempty = FockSpace(Subset(mode for mode in extendedwannierizedempty |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    refunictcellfockfilled = FockSpace(Subset(mode for mode in extendedwannierizedfilled |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    refunictcellfockcourier = FockSpace(Subset(mode for mode in extendedwannierizedcourier |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))

    wannieremptyisometry = globalwannierfunction(firstrgedcorrelations,extendedwannierizedempty[:,refunictcellfockempty])
    wannierfilledisometry = globalwannierfunction(firstrgedcorrelations,extendedwannierizedfilled[:,refunictcellfockfilled])
    wanniercourierisometry = globalwannierfunction(firstrgedcorrelations,extendedwannierizedcourier[:,refunictcellfockcourier])

    couriercorrelations2 = wanniercourierisometry' * firstrgedcorrelations * wanniercourierisometry
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

    shiftedthirdcenter1 = [0,1] ∈ secondrgedspace
    shiftedthirdhexagonalregion1 = thirdhexagonalregion.+shiftedthirdcenter1
    shiftedthirdhexagonalregion1fock = quantize(shiftedthirdhexagonalregion1,1)

    shiftedthirdcenter2 = [1,0] ∈ secondrgedspace
    shiftedthirdhexagonalregion2 = thirdhexagonalregion.+shiftedthirdcenter2
    shiftedthirdhexagonalregion2fock = quantize(shiftedthirdhexagonalregion2,1)

    allregion3 = thirdhexagonalregion+shiftedthirdhexagonalregion1+shiftedthirdhexagonalregion2
    visualize(allregion3)
    visualize(secondrgedcrystal|>getunitcell)
    visualize(intersect(allregion3,firstrgedcrystal|>getunitcell))
    if intersect(allregion3,firstrgedcrystal|>getunitcell) == secondrgedcrystal|>getunitcell
        @info("cover the whole unitcell")
    else
        @error("the allregion cannot cover the whole unitcell ")
    end

    localcorrelations = regioncorrelations(firstrgedcorrelations,thirdhexagonalregionfock)
    localspectrum = localcorrelations|>eigspech
    display(visualize(localspectrum))
    localiso = localcorrelations|>modeselectionbycount(12)

    localwannierseeds = localwannierseedslists(localiso)

    emptyseeds = idmap(thirdhexagonalregionfock, thirdhexagonalregionfock)[:,FockSpace(mode for mode in localwannierseeds[:empty])] 
    filledseeds = idmap(thirdhexagonalregionfock, thirdhexagonalregionfock)[:,FockSpace(mode for mode in localwannierseeds[:filled])] 
    courierseeds = idmap(thirdhexagonalregionfock, thirdhexagonalregionfock)[:,FockSpace(mode for mode in localwannierseeds[:courier])]

    wannierempty = localwannierization(localiso[:empty], emptyseeds)
    wannierfilled = localwannierization(localiso[:filled], filledseeds)
    wanniercourier = localwannierization(localiso[:courier], courierseeds)

    localwannierresults =  Dict(:wannierempty => wannierempty, :wannierfilled => wannierfilled, :wanniercourier => wanniercourier,
                    :emptyseeds => emptyseeds, :filledseeds => filledseeds, :courierseeds => courierseeds,
                    :bempty=>localwannierresults[:bempty],:bfilled=>localwannierresults[:bfilled],:bcourier=>localwannierresults[:bcourier])

    wannierinfos =  Dict(thirdhexagonalregionfock =>localwannierresults)

    thirdshiftedhexagonalregionfocklist = [shiftedthirdhexagonalregion1fock, shiftedthirdhexagonalregion2fock]

    for hexagonalregionfock in thirdshiftedhexagonalregionfocklist
        localcorrelations = regioncorrelations(purifiedcouriercorrelations,hexagonalregionfock)
        localspectrum = localcorrelations|>eigspech
        localspectrum|>visualize
        localiso = localcorrelations|>modeselectionbycount(12)

        localwannierresults = localwannierseedslists(localiso)

        emptyseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in localwannierresults[:empty])] 
        filledseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in localwannierresults[:filled])] 
        courierseeds = idmap(hexagonalregionfock, hexagonalregionfock)[:,FockSpace(mode for mode in localwannierresults[:courier])]

        wannierempty = localwannierization(localiso[:empty], emptyseeds)
        wannierfilled = localwannierization(localiso[:filled], filledseeds)
        wanniercourier = localwannierization(localiso[:courier], courierseeds)

        shiftedlocalwannierresults =  Dict(:wannierempty => wannierempty, :wannierfilled => wannierfilled, :wanniercourier => wanniercourier,
                    :emptyseeds => emptyseeds, :filledseeds => filledseeds, :courierseeds => courierseeds,
                    :bempty=>localwannierseeds[:bempty],:bfilled=>localwannierseeds[:bfilled],:bcourier=>localwannierseeds[:bcourier])

        wannierinfos[hexagonalregionfock] = shiftedlocalwannierresults
    end 

    thirdhexagonalregionfocklist = [thirdhexagonalregionfock, shiftedthirdhexagonalregion1fock, shiftedthirdhexagonalregion2fock]

    extendedwannierizedempty =  sum(wannierinfos[regionfock][:wannierempty] for regionfock in thirdhexagonalregionfocklist)
    extendedwannierizedfilled =  sum(wannierinfos[regionfock][:wannierfilled] for regionfock in thirdhexagonalregionfocklist)
    extendedwannierizedcourier =  sum(wannierinfos[regionfock][:wanniercourier] for regionfock in thirdhexagonalregionfocklist)
            
    origin = [0, 0] ∈ secondrgedspace
    refunictcellfockempty = FockSpace(Subset(mode for mode in extendedwannierizedempty |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    refunictcellfockfilled = FockSpace(Subset(mode for mode in extendedwannierizedfilled |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    refunictcellfockcourier = FockSpace(Subset(mode for mode in extendedwannierizedcourier |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))

    wannieremptyisometry = globalwannierfunction(secondrgedcorrelations,extendedwannierizedempty[:,refunictcellfockempty])
    wannierfilledisometry = globalwannierfunction(secondrgedcorrelations,extendedwannierizedfilled[:,refunictcellfockfilled])
    wanniercourierisometry = globalwannierfunction(secondrgedcorrelations,extendedwannierizedcourier[:,refunictcellfockcourier])

    couriercorrelations3 = wanniercourierisometry' * secondrgedcorrelations * wanniercourierisometry
    couriercorrelationspectrum3 = couriercorrelations3 |> crystalspectrum
    purifiedcorrelationspectrum3 = couriercorrelationspectrum3 |> roundingpurification
    purifiedcouriercorrelations3 = purifiedcorrelationspectrum3 |> CrystalFockMap

    thirdrgedcorrelations = purifiedcouriercorrelations3
    thirdrgedcrystalfock = purifiedcouriercorrelations3|>getoutspace
    thirdrgedcrystal::Crystal = thirdrgedcrystalfock|>getcrystal
    thirdrgedspace::RealSpace = thirdrgedcrystal|>getspace


    thirdrgedcrystal|>getunitcell|>visualize

rg1 = @time gmera(correlations)
rg2 = @time gmera(rg1)