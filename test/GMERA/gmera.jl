using Zipper
using LinearAlgebra
using Combinatorics

# function findmaxofoverlap(truncatedeigenvec, chosenmodes,modes)
#     return maximum(abs.((columns(truncatedeigenvec,FockSpace(chosenmodes))'*columns(truncatedeigenvec,FockSpace(modes))).rep)) 
# end 


# function choosemodesforwannierization(rankedandgroupedmodes,truncatedeigenvec,threshold::Float64,noofmodes::Number)
#     chosenmodes = rankedandgroupedmodes[1]
#     for modes in rankedandgroupedmodes[2:end]
#         if length(chosenmodes)<noofmodes
#             if length(chosenmodes)+length(modes)>noofmodes
#                 @info("more than need")
#             else
#                 max = findmaxofoverlap(truncatedeigenvec, chosenmodes,modes)
#                 if max<threshold
#                     chosenmodes = chosenmodes+modes
#                     chosencolumns = columns(truncatedeigenvec, FockSpace(chosenmodes))
#                     U, Σ, Vt = svd(chosencolumns)
#                     minsvdvalue::Number = minimum(v for (_, v) in Σ)
#                     @info("min svdvalue during checking after", minsvdvalue)
#                 else
#                     @info("reject this candidate due to small max", max)
#                 end
#             end
#         else
#             return chosenmodes
#         end
#     end
#     @error("threshold too large")
# end

# function minsvdvalue(trunvatedevec,modes::Subset{Mode})
#     u,s,vt = svd(columns(trunvatedevec,FockSpace(modes)))
#     return minimum(v for (_, v) in s)
# end

# correlations,H = generatesystem(-0.3 + 0im,0.3 + 0im,-1 + 0im,-(1/(3*sqrt(3)))im,32)
correlations,H = generatesystem( 0im,0im,-1 + 0im,0im,32)

crystalfock = correlations|>getoutspace

scale = Scale([4 0; 0 4], crystalfock|>getcrystal|>getspace)
@info("Performing rgblocking...",scale)
@info("Generating rgblocking transformation...")
blocker = @time scale * crystalfock
@info("Performing rgblocking on correlations...")
blockedcorrelations = @time blocker * correlations * blocker'
blockedcrystal = blockedcorrelations|>getoutspace|>getcrystal

Asublatticeunitcellmodes,Bsublatticeunitcellmodes = generaterefAandBsublatticeunitcellmodes(blockedcrystal)
Asublatticeoffsets = Subset((mode|>getattr(:b)) for mode in Asublatticeunitcellmodes)
Bsublatticeoffsets = Subset((mode|>getattr(:b)) for mode in Bsublatticeunitcellmodes)
visualize(Asublatticeoffsets)
visualize(Bsublatticeoffsets)

# crystalfock = blockedcorrelations|>getoutspace
#     crystal::Crystal = crystalfock|>getcrystal
#     space::RealSpace = crystal|>getspace

#     @info("Computing local correlations...")
#     firstcenter = [0,0] ∈ space
#     firsthexagonalregion = gethexagonalregion(crystal=crystal, center=firstcenter, metricspace=space)
#     firsthexagonalregionfock = quantize(firsthexagonalregion,1)
    
#     shiftedfirstcenter1 = [1,0] ∈ space
#     shiftedfirsthexagonalregion1 = firsthexagonalregion.+shiftedfirstcenter1
#     shiftedfirsthexagonalregion1fock = quantize(shiftedfirsthexagonalregion1,1)
    
#     shiftedfirstcenter2 = [0,1] ∈ space
#     shiftedfirsthexagonalregion2 = firsthexagonalregion.+shiftedfirstcenter2
#     shiftedfirsthexagonalregion2fock = quantize(shiftedfirsthexagonalregion2,1)
    
#     shiftedfirstcenter3 = [1,1] ∈ space
#     shiftedfirsthexagonalregion3 = firsthexagonalregion.+shiftedfirstcenter3
#     shiftedfirsthexagonalregion3fock = quantize(shiftedfirsthexagonalregion3,1)
    
#     allregion = firsthexagonalregion+shiftedfirsthexagonalregion1+shiftedfirsthexagonalregion2+shiftedfirsthexagonalregion3
#     if intersect(allregion,crystal|>getunitcell) == crystal|>getunitcell
#         @info("cover the whole unitcell")
#     else
#         @error("the allregion cannot cover the whole unitcell ")
#     end

#     localcorrelations = regioncorrelations(blockedcorrelations,firsthexagonalregionfock)
#     localspectrum = localcorrelations|>eigspech
#     display(localspectrum|>visualize)
#     filleddict,emptydict = separatefilledandemptymodes(localspectrum)
#     rankedandgroupedfilledmodes = sortgroupdictwifvaluefilled(filleddict,0.05)
#     rankedandgroupedemptymodes = sortgroupdictwifvalueempty(emptydict,0.05)

#     rsmodedict = Dict(mode=>norm(euclidean((mode|>getattr(:b))+(mode|>getattr(:r))-firstcenter)) for mode in localspectrum |> geteigenvectors|>getoutspace|>orderedmodes)
#     sortedrsmode = sortgroupdictwifdist(rsmodedict,false)
#     frozenrsmodes, courierrsmodes = findfrozenandcourierrsmode(sortedrsmode,36)
#     filledrsmodes,emptyrsmodes = distinguishABsublatticemodesforhexagon(frozenrsmodes,Asublatticeunitcellmodes,Bsublatticeunitcellmodes)

#     filledseeds = idmap(firsthexagonalregionfock, firsthexagonalregionfock)[:,FockSpace(mode for mode in filledrsmodes)]
#     emptyseeds = idmap(firsthexagonalregionfock, firsthexagonalregionfock)[:,FockSpace(mode for mode in emptyrsmodes)]
#     frozenseeds = idmap(firsthexagonalregionfock, firsthexagonalregionfock)[:,FockSpace(mode for mode in frozenrsmodes)]
#     courierseeds = idmap(firsthexagonalregionfock, firsthexagonalregionfock)[:,FockSpace(mode for mode in courierrsmodes)]

#     truncatedeigenvectofilled = rows((localspectrum|>geteigenvectors),FockSpace(filledrsmodes))
#     truncatedeigenvectoempty = rows((localspectrum|>geteigenvectors),FockSpace(emptyrsmodes))

#     svd(columns(truncatedeigenvectofilled,FockSpace(sum(rankedandgroupedcouriermodeswifoverlaptofilled[1:9]))))
#     svd(columns(truncatedeigenvectofilled,FockSpace(sum(rankedandgroupedfilledmodes[1:3]))))

#     u1,s1,vt1 = svd(columns(truncatedeigenvectofilled,FockSpace(rankedandgroupedcouriermodeswifoverlaptofilled[1])))
#     u2,s2,vt2 = svd(columns(truncatedeigenvectofilled,FockSpace(rankedandgroupedcouriermodeswifoverlaptofilled[2])))
#     maximum(abs.((columns(truncatedeigenvectofilled,FockSpace(sum(rankedandgroupedcouriermodeswifoverlaptofilled[1:8])))'*columns(truncatedeigenvectofilled,FockSpace(rankedandgroupedcouriermodeswifoverlaptofilled[11]))).rep))
    
#     for r in range(1,20)
#         println((columns(truncatedeigenvectofilled,FockSpace(sum(rankedandgroupedfilledmodes[1:2])))'*columns(truncatedeigenvectofilled,FockSpace(rankedandgroupedfilledmodes[r]))).rep)
#     end

#     # chosencouriermodes = chosemodesforwannierization(rankedandgroupedcouriermodeswifoverlaptofilled,truncatedeigenvectocourier,0.5,div(60,2)) + chosemodesforwannierization(rankedandgroupedcouriermodeswifoverlaptoempty,truncatedeigenvectocourier,0.5,div(60,2))
#     rankedandgroupedfilledmodeswifoverlap = [pair[2] for pair in sort([calculateavgoverlapofmodes(truncatedeigenvectofilled,modes) for modes in rankedandgroupedfilledmodes],by=first,rev=true)]
#     rankedandgroupedemptymodeswifoverlap = [pair[2] for pair in sort([calculateavgoverlapofmodes(truncatedeigenvectoempty,modes) for modes in rankedandgroupedemptymodes],by=first,rev=true)]

#     chosenfilledmodes = chosemodesforwannierization(rankedandgroupedfilledmodeswifoverlap,truncatedeigenvectofilled,0.5,div(36,2))
#     chosenemptymodes = chosemodesforwannierization(rankedandgroupedemptymodeswifoverlap,truncatedeigenvectoempty,0.1,div(36,2))
#     chosenfrozenmodes = chosenfilledmodes+chosenemptymodes
#     chosencouriermodes = (localspectrum |> geteigenmodes)-chosenfrozenmodes

#     localisofilled = columns(localspectrum |> geteigenvectors, FockSpace(chosenfilledmodes))
#     localisoempty = columns(localspectrum |> geteigenvectors, FockSpace(chosenemptymodes))
#     localisofrozen = columns(localspectrum |> geteigenvectors, FockSpace(chosenfrozenmodes))
#     localisocourier = columns(localspectrum |> geteigenvectors, FockSpace(chosencouriermodes))

#     wannierfilled = localwannierization(localisofilled, filledseeds)
#     wannierempty = localwannierization(localisoempty, emptyseeds)
#     wannierfrozen = localwannierization(localisofrozen, frozenseeds)
#     wanniercourier = localwannierization(localisocourier, courierseeds)

#     (wannierfrozen'*localcorrelations*wannierfrozen)|>eigspec|>visualize|>display

firstrgedcorrelations = firstgmerastep(blockedcorrelations,0.9,36,Asublatticeunitcellmodes,Bsublatticeunitcellmodes)

secondrgedcorrelations = secondgmerastep(firstrgedcorrelations,0.9,30,Asublatticeunitcellmodes,Bsublatticeunitcellmodes)

thirdrgedcorrelations = thirdgmerastep(secondrgedcorrelations,0.4,24,Asublatticeunitcellmodes,Bsublatticeunitcellmodes)