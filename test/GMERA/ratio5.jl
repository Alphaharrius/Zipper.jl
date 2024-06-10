using Zipper
using LinearAlgebra

function calculateavgoverlapofmodes(modes,spectrum::EigenSpectrum,rsmodes::Subset{Mode})
    overlaps = [(norm(rows(columns(spectrum |> geteigenvectors, FockSpace(mode)),FockSpace(rsmodes))))^2 for mode in modes]
    avgoverlap = sum(overlaps)/length(overlaps)
    return tuple(avgoverlap,modes)
end

function calculateoverlapofmode(evalwifmode,spectrum::EigenSpectrum,rsmodes::Subset{Mode})
    overlap = (norm(rows(columns(spectrum |> geteigenvectors, FockSpace(evalwifmode[2])),FockSpace(rsmodes))))^2
    return tuple(evalwifmode[2],evalwifmode[1],overlap)
end

function sortdictwifvaluefirst(dict::Dict{Mode, Float64},rev::Bool)
    sortedtupledata = sort([(value,key) for (key,value) in dict],rev=rev,by=first)
    return sortedtupledata
end

function sortgroupdictwifvalue(dict::Dict{Mode, Float64},rev::Bool)
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

# function sortgroupdictwifvalue(dict::Dict{Mode, Float64},rev::Bool)
#     sortedtupledata = sort([(value,key) for (key,value) in dict],rev=rev,by=first)
#     refvalue = sortedtupledata[1][1]
#     result = []
#     subresult = Subset(sortedtupledata[1][2])
#     for pair in sortedtupledata
#         if pair[1]==refvalue
#             subresult = subresult+Subset(pair[2])
#         else
#             append!(result,[subresult])
#             refvalue = pair[1]
#             subresult = Subset(pair[2])
#         end
#     end
#     append!(result,[subresult])

#     return result
# end

function sortgroupdictwifdist(dict::Dict{Mode, Float64},rev::Bool)
    sortedtupledata = sort([(value,key) for (key,value) in dict],rev=rev,by=first)
    refvalue = round(sortedtupledata[1][1],digits=8)
    result = []
    subresult = Subset(sortedtupledata[1][2])
    for pair in sortedtupledata
        if round(pair[1],digits=8)==refvalue
            subresult = subresult+Subset(pair[2])
        else
            append!(result,tuple([refvalue,subresult]))
            refvalue = round(pair[1],digits=8)
            subresult = Subset(pair[2])
        end
    end
    append!(result,tuple([refvalue,subresult]))

    return result
end

function sortgroupdictwifvaluefirst(dict::Dict{Mode, Float64},rev::Bool)
    sortedtupledata = sort([(value,key) for (key,value) in dict],rev=rev,by=first)
    # return sortedtupledata
    refvalue = round(sortedtupledata[1][1],digits=8)
    result = []
    subresult = Subset(sortedtupledata[1][2])+Subset(sortedtupledata[2][2])
    for pair in sortedtupledata
        if round(pair[1],digits=8)==refvalue
            subresult = subresult+Subset(pair[2])
        else
            append!(result,tuple([refvalue,subresult]))
            refvalue = round(pair[1],digits=8)
            subresult = Subset(pair[2])
        end
    end
    append!(result,tuple([refvalue,subresult]))

    return result
end

function sortgroupdictwifvaluesec(dict::Dict{Mode, Float64},rev::Bool)
    sortedtupledata = sort([(value,key) for (key,value) in dict],rev=rev,by=first)
    # return sortedtupledata
    refvalue = sortedtupledata[1][1]
    result = []
    subresult = Subset(sortedtupledata[1][2])+Subset(sortedtupledata[2][2])
    for pair in sortedtupledata
        if pair[1]==refvalue
            subresult = subresult+Subset(pair[2])
        else
            append!(result,tuple([refvalue,subresult]))
            refvalue = pair[1]
            subresult = Subset(pair[2])
        end
    end
    append!(result,tuple([refvalue,subresult]))

    return result
end

function findmodeswiftotalno(sortedgroup,totalno::Number)
    ref = totalno 
    result = [] 
    for data in sortedgroup
        condition = ref-length(data[2])
        if condition>0
            ref = condition
            append!(result,tuple([data[1],data[2]]))
        elseif condition<0
            ref = ref
            @info("more than the number of modes we need")
        else
            ref = condition
            append!(result,tuple([data[1],data[2]]))
            return result
        end
    end
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
crystal = Crystal(unitcell, [20, 20])
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
localspectrum = localcorrelations|>eigspech
display(visualize(localspectrum))
filleddict,emptydict = separatefilledandemptymodes(localspectrum)

sortedgroupfilledwifeval = sortgroupdictwifvalue(filleddict,false)
sortedgroupemptywifeval = sortgroupdictwifvalue(emptydict,true)


rsmodedict = Dict(mode=>norm(euclidean((mode|>getattr(:b))+(mode|>getattr(:r))-firstcenter)) for mode in localspectrum |> geteigenvectors|>getoutspace|>orderedmodes)
sortedrsmode = sortgroupdictwifdist(rsmodedict,false)
frozenrsmodes = sum([pair[2] for pair in sortedrsmode[1:6]])
frozenrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in frozenrsmodes)
visualize(frozenrsmodesoffset)
courierrsmodes = sum([pair[2] for pair in sortedrsmode[7:end]])
courierrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in courierrsmodes)
visualize(courierrsmodesoffset)

sortedgroupfilledwifoverlap = sort([calculateavgoverlapofmodes(evalwifmodes,localspectrum,frozenrsmodes) for evalwifmodes in sortedgroupfilledwifeval],rev=true,by=x->x[1])
sortedgroupemptywifoverlap = sort([calculateavgoverlapofmodes(evalwifmodes,localspectrum,frozenrsmodes) for evalwifmodes in sortedgroupemptywifeval],rev=true,by=x->x[1])

chosensortedgroupfilledwifoverlap = findmodeswiftotalno(sortedgroupfilledwifoverlap,div(length(frozenrsmodes),2))
chosensortedgroupemptywifoverlap = findmodeswiftotalno(sortedgroupemptywifoverlap,div(length(frozenrsmodes),2))
chosensortedgroupfilledwifoverlap
chosensortedgroupemptywifoverlap

chosenfilledmodes = sum([data[2] for data in chosensortedgroupfilledwifoverlap])
chosenemptymodes = sum([data[2] for data in chosensortedgroupemptywifoverlap])
chosenfrozenmodes = chosenfilledmodes+chosenemptymodes
chosencouriermodes = (localspectrum |> geteigenmodes)-chosenfrozenmodes

courierseeds = idmap(firsthexagonalregionfock, firsthexagonalregionfock)[:,FockSpace(mode for mode in courierrsmodes)]
localisocourier = columns(localspectrum |> geteigenvectors, FockSpace(chosencouriermodes))
wanniercourier = localwannierization(localisocourier, courierseeds)

frozenseeds = idmap(firsthexagonalregionfock, firsthexagonalregionfock)[:,FockSpace(mode for mode in frozenrsmodes)]
localisoforzen = columns(localspectrum |> geteigenvectors, FockSpace(chosenfrozenmodes))
wannierfrozen = localwannierization(localisoforzen, frozenseeds)

(wanniercourier'*localcorrelations*wanniercourier)|>eigspech|>visualize
(wannierfrozen'*localcorrelations*wannierfrozen)|>eigspech|>visualize


localwannierresults =  Dict(:wannierfrozen => wannierfrozen, :wanniercourier => wanniercourier,
                    :frozenseeds => frozenseeds,:courierseeds => courierseeds)

wannierinfos =  Dict(firsthexagonalregionfock=>localwannierresults)

firstshiftedhexagonalregionfocklist = [(shiftedfirstcenter1,shiftedfirsthexagonalregion1fock), (shiftedfirstcenter2,shiftedfirsthexagonalregion2fock), (shiftedfirstcenter3,shiftedfirsthexagonalregion3fock)]

for (shifted,hexagonalregionfock) in firstshiftedhexagonalregionfocklist
    localcorrelations = regioncorrelations(blockedcorrelations,hexagonalregionfock)
    localspectrum = localcorrelations|>eigspech
    filleddict,emptydict = separatefilledandemptymodes(localspectrum)
    
    sortedgroupfilledwifeval = sortgroupdictwifvalue(filleddict,false)
    sortedgroupemptywifeval = sortgroupdictwifvalue(emptydict,true)

    rsmodedict = Dict(mode=>norm(euclidean((mode|>getattr(:b))+(mode|>getattr(:r))-(firstcenter+shifted))) for mode in localspectrum |> geteigenvectors|>getoutspace|>orderedmodes)
    sortedrsmode = sortgroupdictwifdist(rsmodedict,false)
    frozenrsmodes = sum([pair[2] for pair in sortedrsmode[1:6]])
    frozenrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in frozenrsmodes)
    courierrsmodes = sum([pair[2] for pair in sortedrsmode[7:end]])
    courierrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in courierrsmodes)

    sortedgroupfilledwifoverlap = sort([calculateavgoverlapofmodes(evalwifmodes,localspectrum,frozenrsmodes) for evalwifmodes in sortedgroupfilledwifeval],rev=true,by=x->x[1])
    sortedgroupemptywifoverlap = sort([calculateavgoverlapofmodes(evalwifmodes,localspectrum,frozenrsmodes) for evalwifmodes in sortedgroupemptywifeval],rev=true,by=x->x[1])

    chosensortedgroupfilledwifoverlap = findmodeswiftotalno(sortedgroupfilledwifoverlap,div(length(frozenrsmodes),2))
    chosensortedgroupemptywifoverlap = findmodeswiftotalno(sortedgroupemptywifoverlap,div(length(frozenrsmodes),2))

    chosenfilledmodes = sum([data[2] for data in chosensortedgroupfilledwifoverlap])
    chosenemptymodes = sum([data[2] for data in chosensortedgroupemptywifoverlap])
    chosenfrozenmodes = chosenfilledmodes+chosenemptymodes
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
# wannierisometry = globalwannierfunction(blockedcorrelations,extendedwannierizedfrozen[:,refunictcellfockfrozen]+extendedwannierizedcourier[:,refunictcellfockcourier])

couriercorrelations = wanniercourierisometry' * blockedcorrelations * wanniercourierisometry
couriercorrelations|>eigspech|>visualize
couriercorrelationspectrum = couriercorrelations |> crystalspectrum
purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
purifiedcouriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap

# refcorrelations = wannierisometry' * blockedcorrelations * wannierisometry

frozencorrelations = wannierforzenisometry' * blockedcorrelations * wannierforzenisometry
frozencorrelations|>eigspech|>visualize

firstrgedcorrelations = purifiedcouriercorrelations
firstrgedcrystalfock = purifiedcouriercorrelations|>getoutspace
firstrgedcrystal::Crystal = firstrgedcrystalfock|>getcrystal
firstrgedspace::RealSpace = firstrgedcrystal|>getspace

refrgedcrystalfock = refcorrelations|>getoutspace
refrgedcrystal::Crystal = refrgedcrystalfock|>getcrystal
refrgedspace::RealSpace = refrgedcrystal|>getspace

@info("Computing local correlations...")
secondcenter = [2/3,-1/3] ∈ firstrgedspace
secondhexagonalregion = gethexagonalregion(crystal=firstrgedcrystal, center=secondcenter, metricspace=firstrgedspace)
visualize(secondhexagonalregion)
secondhexagonalregionfock = quantize(secondhexagonalregion,1)

# secondcenter = [2/3,-1/3] ∈ firstrgedspace
# secondhexagonalregion = gethexagonalregion(crystal=refrgedcrystal, center=secondcenter, metricspace=refrgedspace)
# visualize(secondhexagonalregion)
# secondhexagonalregionfock = quantize(secondhexagonalregion,1)

# localcorrelations = regioncorrelations(refcorrelations,secondhexagonalregionfock)
# localspectrum = localcorrelations|>eigspech|>visualize

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

    sortedgroupfilledwifeval = sortgroupdictwifvalue(filleddict,false)

    sortedgroupemptywifeval = sortgroupdictwifvalue(emptydict,true)

    # for modes in sortedgroupfilledwifeval
    #     for (ind,m) in enumerate(modes)
    #         println((ind,(localspectrum|>geteigenvalues)[m]))
    #     end
    # end

    # for modes in sortedgroupemptywifeval
    #     for (ind,m) in enumerate(modes)
    #         println((ind,(localspectrum|>geteigenvalues)[m]))
    #     end
    # end
    
    rsmodedict = Dict(mode=>norm(euclidean((mode|>getattr(:b))+(mode|>getattr(:r))-secondcenter)) for mode in localspectrum |> geteigenvectors|>getoutspace|>orderedmodes)
    sortedrsmode = sortgroupdictwifdist(rsmodedict,false)
    frozenrsmodes = sum([pair[2] for pair in sortedrsmode[1:8]])
    frozenrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in frozenrsmodes)
    visualize(frozenrsmodesoffset)
    courierrsmodes = sum([pair[2] for pair in sortedrsmode[9:end]])
    courierrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in courierrsmodes)
    visualize(courierrsmodesoffset)

    sortedgroupfilledwifoverlap = sort([calculateavgoverlapofmodes(evalwifmodes,localspectrum,frozenrsmodes) for evalwifmodes in sortedgroupfilledwifeval],rev=true,by=x->x[1])
    sortedgroupemptywifoverlap = sort([calculateavgoverlapofmodes(evalwifmodes,localspectrum,frozenrsmodes) for evalwifmodes in sortedgroupemptywifeval],rev=true,by=x->x[1])

    # chosensortedgroupfilledwifoverlap = findmodeswiftotalno(sortedgroupfilledwifoverlap,div(length(frozenrsmodes),2))
    # chosensortedgroupemptywifoverlap = findmodeswiftotalno(sortedgroupemptywifoverlap,div(length(frozenrsmodes),2))

    # chosenfilledmodes = sum([data[2] for data in chosensortedgroupfilledwifoverlap])
    # chosenemptymodes = sum([data[2] for data in chosensortedgroupemptywifoverlap])
    # chosenfrozenmodes = chosenfilledmodes+chosenemptymodes
    # chosencouriermodes = (localspectrum |> geteigenmodes)-chosenfrozenmodes

    chosenfilledmodes = sum([data[2] for data in sortedgroupfilledwifoverlap[1:21]])
    chosenemptymodes = sum([data[2] for data in sortedgroupemptywifoverlap[1:20]])
    chosenfrozenmodes = chosenfilledmodes+chosenemptymodes
    chosencouriermodes = (localspectrum |> geteigenmodes)-chosenfrozenmodes
    

    sortedgroupfilledwifoverlap
    for pair in sortedgroupfilledwifoverlap
        for (ind,m) in enumerate(pair[2])
            println((ind,(localspectrum|>geteigenvalues)[m]))
        end
    end

    for pair in sortedgroupemptywifoverlap
        for (ind,m) in enumerate(pair[2])
            println((ind,(localspectrum|>geteigenvalues)[m]))
        end
    end

    for pair in sortedgroupemptywifoverlap
        println((pair))
    end

    for m in sortedgroupemptywifoverlap[3][2]
        println((localspectrum|>geteigenvalues)[m])
    end
    
    courierseeds = idmap(secondhexagonalregionfock, secondhexagonalregionfock)[:,FockSpace(mode for mode in courierrsmodes)]
    # localisocourier = columns(localspectrum |> geteigenvectors, FockSpace(chosencouriermodes))
    # wanniercourier = localwannierization(localisocourier, courierseeds)

    frozenseeds = idmap(secondhexagonalregionfock, secondhexagonalregionfock)[:,FockSpace(mode for mode in frozenrsmodes)]
    # localisoforzen = columns(localspectrum |> geteigenvectors, FockSpace(chosenfrozenmodes))
    # wannierfrozen = localwannierization(localisoforzen, frozenseeds)

    wanniercourierref = localwannierization(localspectrum |> geteigenvectors, courierseeds)
    wannierfrozenref = localwannierization(localspectrum |> geteigenvectors, frozenseeds)

    sortedfrozen = sort([(norm(rows((localspectrum |> geteigenvectors)'*wannierfrozenref,FockSpace(mode)))^2,mode,(localspectrum |> geteigenvalues)[mode]) for mode in ((localspectrum |> geteigenvectors)'*wannierfrozenref)|>getoutspace],rev=true,by=first)
    sortedcourier =  sort([(norm(rows((localspectrum |> geteigenvectors)'*wanniercourierref,FockSpace(mode)))^2,mode,(localspectrum |> geteigenvalues)[mode]) for mode in ((localspectrum |> geteigenvectors)'*wanniercourierref)|>getoutspace],rev=true,by=first)

    sortedfilled = [data for data in sortedfrozen if data[3]<0.5]
    sortedempty = [data for data in sortedfrozen if data[3]>0.5]

    for data in sortedfilled
        println(data)
    end

    for data in sortedempty
        println(data)
    end

    chosenfrozenmodes = sum(Subset(pair[2]) for pair in sortedfilled[1:20])+sum(Subset(pair[2]) for pair in sortedempty[1:22])
    chosencouriermodes = (localspectrum |> geteigenmodes)-chosenfrozenmodes

    

    localisoforzen = columns(localspectrum |> geteigenvectors, FockSpace(chosenfrozenmodes))
    localisocourier = columns(localspectrum |> geteigenvectors, FockSpace(chosencouriermodes))
    wannierfrozen = localwannierization(localisoforzen, frozenseeds)
    wanniercourier = localwannierization(localisocourier, courierseeds)


    visualize(wanniercourier'*localcorrelations*wanniercourier|>eigspech)
    visualize(wannierfrozen'*localcorrelations*wannierfrozen|>eigspech)

    localwannierresults =  Dict(:wannierfrozen => wannierfrozen, :wanniercourier => wanniercourier,
                    :frozenseeds => frozenseeds,:courierseeds => courierseeds)

    wannierinfos =  Dict(secondhexagonalregionfock=>localwannierresults)

    secondshiftedhexagonalregionfocklist = [(shiftedsecondcenter1,shiftedsecondhexagonalregion1fock), (shiftedsecondcenter2,shiftedsecondhexagonalregion2fock)]

    for (shifted,hexagonalregionfock) in secondshiftedhexagonalregionfocklist
        localcorrelations = regioncorrelations(firstrgedcorrelations,hexagonalregionfock)
        localspectrum = localcorrelations|>eigspech
        filleddict,emptydict = separatefilledandemptymodes(localspectrum)

        sortedgroupfilledwifeval = sortgroupdictwifvalue(filleddict,false)
        sortedgroupemptywifeval = sortgroupdictwifvalue(emptydict,true)
        
        rsmodedict = Dict(mode=>norm(euclidean((mode|>getattr(:b))+(mode|>getattr(:r))-(secondcenter+shifted))) for mode in localspectrum |> geteigenvectors|>getoutspace|>orderedmodes)
        sortedrsmode = sortgroupdictwifdist(rsmodedict,false)
        frozenrsmodes = sum([pair[2] for pair in sortedrsmode[1:8]])
        frozenrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in frozenrsmodes)
        visualize(frozenrsmodesoffset)
        courierrsmodes = sum([pair[2] for pair in sortedrsmode[9:end]])
        courierrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in courierrsmodes)
        visualize(courierrsmodesoffset)

        sortedgroupfilledwifoverlap = sort([calculateavgoverlapofmodes(evalwifmodes,localspectrum,frozenrsmodes) for evalwifmodes in sortedgroupfilledwifeval],rev=true,by=x->x[1])
        sortedgroupemptywifoverlap = sort([calculateavgoverlapofmodes(evalwifmodes,localspectrum,frozenrsmodes) for evalwifmodes in sortedgroupemptywifeval],rev=true,by=x->x[1])


        chosensortedgroupfilledwifoverlap = findmodeswiftotalno(sortedgroupfilledwifoverlap,div(length(frozenrsmodes),2))
        chosensortedgroupemptywifoverlap = findmodeswiftotalno(sortedgroupemptywifoverlap,div(length(frozenrsmodes),2))

        chosenfilledmodes = sum([data[2] for data in chosensortedgroupfilledwifoverlap])
        chosenemptymodes = sum([data[2] for data in chosensortedgroupemptywifoverlap])
        chosenfrozenmodes = chosenfilledmodes+chosenemptymodes
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

    sortedgroupfilledwifeval = sortgroupdictwifvalue(filleddict,false)
    sortedgroupemptywifeval = sortgroupdictwifvalue(emptydict,true)

    rsmodedict = Dict(mode=>norm(euclidean((mode|>getattr(:b))+(mode|>getattr(:r))-thirdcenter)) for mode in localspectrum |> geteigenvectors|>getoutspace|>orderedmodes)
    sortedrsmode = sortgroupdictwifdist(rsmodedict,false)
    frozenrsmodes = sum([pair[2] for pair in sortedrsmode[1:7]])
    frozenrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in frozenrsmodes)
    visualize(frozenrsmodesoffset)
    courierrsmodes = sum([pair[2] for pair in sortedrsmode[8:end]])
    courierrsmodesoffset = Subset((mode|>getattr(:b))+(mode|>getattr(:r)) for mode in courierrsmodes)
    visualize(courierrsmodesoffset)

    sortedgroupfilledwifoverlap = sort([calculateavgoverlapofmodes(evalwifmodes,localspectrum,frozenrsmodes) for evalwifmodes in sortedgroupfilledwifeval],rev=true,by=x->x[1])
    sortedgroupemptywifoverlap = sort([calculateavgoverlapofmodes(evalwifmodes,localspectrum,frozenrsmodes) for evalwifmodes in sortedgroupemptywifeval],rev=true,by=x->x[1])

    chosensortedgroupfilledwifoverlap = findmodeswiftotalno(sortedgroupfilledwifoverlap,div(length(frozenrsmodes),2))
    chosensortedgroupemptywifoverlap = findmodeswiftotalno(sortedgroupemptywifoverlap,div(length(frozenrsmodes),2))

    chosenfilledmodes = sum([data[2] for data in chosensortedgroupfilledwifoverlap])
    chosenemptymodes = sum([data[2] for data in chosensortedgroupemptywifoverlap])
    chosenfrozenmodes = chosenfilledmodes+chosenemptymodes
    chosencouriermodes = (localspectrum |> geteigenmodes)-chosenfrozenmodes

    courierseeds = idmap(thirdhexagonalregionfock, thirdhexagonalregionfock)[:,FockSpace(mode for mode in courierrsmodes)]
    localisocourier = columns(localspectrum |> geteigenvectors, FockSpace(chosencouriermodes))
    wanniercourier = localwannierization(localisocourier, courierseeds)

    frozenseeds = idmap(thirdhexagonalregionfock, thirdhexagonalregionfock)[:,FockSpace(mode for mode in frozenrsmodes)]
    localisoforzen = columns(localspectrum |> geteigenvectors, FockSpace(chosenfrozenmodes))
    wannierfrozen = localwannierization(localisoforzen, frozenseeds)

    visualize(wanniercourier'*localcorrelations*wanniercourier|>eigspech)

    frozenseeds = idmap(thirdhexagonalregionfock, thirdhexagonalregionfock)[:,FockSpace(mode for mode in frozenrsmodes)]
    localisoforzen = columns(localspectrum |> geteigenvectors, FockSpace(chosenfrozenmodes))
    wannierfrozen = localwannierization(localisoforzen, frozenseeds)

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

    visualize(couriercorrelations|>getinspace|>getcrystal|>getunitcell)

