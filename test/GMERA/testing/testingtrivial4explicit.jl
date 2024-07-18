using Zipper
using LinearAlgebra
using Plots

setmaxthreads(8)

systemsize=32

onsite = -0
    correlations,H = generatesystem(onsite+0im, -onsite+0im,-1 + 0im,0im,systemsize)
    crystalfock = correlations|>getoutspace

    scale = Scale([2 0; 0 2], crystalfock|>getcrystal|>getspace)
    @info("Performing intermediateblocking...",scale)
    @info("Generating blocking transformation...")
    intermediateblocker = @time scale * crystalfock
    @info("Performing intermediateblocking on correlations...")
    intermediateblockedcorrelations = @time intermediateblocker * correlations * intermediateblocker'

    intermediateblockedcrystal = intermediateblockedcorrelations|>getoutspace|>getcrystal
    intermediateblockedcrystalfock = intermediateblockedcorrelations|>getoutspace
    intermediateblockedspace::RealSpace = intermediateblockedcrystal|>getspace

    refrot = inv([2/3 -1/3; -1/3 2/3]')
    refcenter = [0,0] ∈ intermediateblockedspace

    intermediatehexagonalregion = gethexagonalregion(rot = refrot,crystal=intermediateblockedcrystal, center=refcenter, metricspace=intermediateblockedspace)

    scale = Scale([2 0; 0 2], intermediateblockedcrystalfock|>getcrystal|>getspace)
    @info("Performing gmerablocking...",scale)
    @info("Generating gmera blocking transformation...")
    blocker = @time scale * intermediateblockedcrystalfock
    @info("Performing gmera blocking on correlations...")
    blockedcorrelations = @time blocker * intermediateblockedcorrelations * blocker'

    blockedcrystal = blockedcorrelations|>getoutspace|>getcrystal
    blockedcrystalfock = blockedcorrelations|>getoutspace
    blockedspace::RealSpace = blockedcrystal|>getspace


    refrot = inv([2/3 -1/3; -1/3 2/3]')
    c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
    c3 = c6^2

    firstcenter = [0,0] ∈ blockedspace
    firsthexagonalregion = gethexagonalregion(rot = refrot,crystal=blockedcrystal, center=firstcenter, metricspace=blockedspace)
    firsthexagonalregionfock = quantize(firsthexagonalregion,1)
    rgshiftedcenter1 = [2/3,-1/3] ∈ blockedspace
    firstrgshiftedhexagonalregion1 = firsthexagonalregion.+rgshiftedcenter1
    rgshiftedcenter2 = [1/3,1/3] ∈ blockedspace
    firstrgshiftedhexagonalregion2 = firsthexagonalregion.+rgshiftedcenter2

    innerhexagonalregion = checkintersectbyeculidean(intermediatehexagonalregion,firsthexagonalregion)
    innerhexagonalregionfock = quantize(innerhexagonalregion,1)

    #siteAregion1
    siteAregion1 = intersect(intersect(firsthexagonalregion,firstrgshiftedhexagonalregion1 ),firstrgshiftedhexagonalregion2)
    siteAregion1closetriangle1 =  intersect(siteAregion1,innerhexagonalregion)

    siteAregion1bdytriangle1 = siteAregion1closetriangle1.+(1/2)*rgshiftedcenter1
    siteAregion1bdytriangle2 = siteAregion1closetriangle1.+(1/2)*rgshiftedcenter2
    siteAregion1bdytriangle1fock = quantize(siteAregion1bdytriangle1,1)
    siteAregion1bdytriangle2fock = quantize(siteAregion1bdytriangle2,1)

    siteAregion1centertriangle = siteAregion1 - siteAregion1closetriangle1 - siteAregion1bdytriangle1 - siteAregion1bdytriangle2
    siteAregion1centertrianglefock = quantize(siteAregion1centertriangle,1)

    #siteAregion2
    siteAregion2 = c3*siteAregion1
    siteAregion2closetriangle1 =  c3*siteAregion1closetriangle1

    siteAregion2bdytriangle1 = c3*siteAregion1bdytriangle1
    siteAregion2bdytriangle2 = c3*siteAregion1bdytriangle2
    siteAregion2bdytriangle1fock = quantize(siteAregion2bdytriangle1,1)
    siteAregion2bdytriangle2fock = quantize(siteAregion2bdytriangle2,1)

    siteAregion2centertriangle = c3*siteAregion1centertriangle
    siteAregion2centertrianglefock = quantize(siteAregion2centertriangle,1)

    #siteAregion3
    siteAregion3 = c3*siteAregion2
    siteAregion3closetriangle1 =  c3*siteAregion2closetriangle1

    siteAregion3bdytriangle1 = c3*siteAregion2bdytriangle1
    siteAregion3bdytriangle2 = c3*siteAregion2bdytriangle2
    siteAregion3bdytriangle1fock = quantize(siteAregion3bdytriangle1,1)
    siteAregion3bdytriangle2fock = quantize(siteAregion3bdytriangle2,1)

    siteAregion3centertriangle = c3*siteAregion2centertriangle
    siteAregion3centertrianglefock = quantize(siteAregion3centertriangle,1)

    #siteBregion1
    siteBregion1 = c6*siteAregion1
    siteBregion1closetriangle1 =  c6*siteAregion1closetriangle1

    siteBregion1bdytriangle1 = c6*siteAregion1bdytriangle1
    siteBregion1bdytriangle2 = c6*siteAregion1bdytriangle2
    siteBregion1bdytriangle1fock = quantize(siteBregion1bdytriangle1,1)
    siteBregion1bdytriangle2fock = quantize(siteBregion1bdytriangle2,1)

    siteBregion1centertriangle = c6*siteAregion1centertriangle
    siteBregion1centertrianglefock = quantize(siteBregion1centertriangle,1)

    #siteBregion2
    siteBregion2 = c3*siteBregion1
    siteBregion2closetriangle1 =  c3*siteBregion1closetriangle1

    siteBregion2bdytriangle1 = c3*siteBregion1bdytriangle1
    siteBregion2bdytriangle2 = c3*siteBregion1bdytriangle2
    siteBregion2bdytriangle1fock = quantize(siteBregion2bdytriangle1,1)
    siteBregion2bdytriangle2fock = quantize(siteBregion2bdytriangle2,1)

    siteBregion2centertriangle = c3*siteBregion1centertriangle
    siteBregion2centertrianglefock = quantize(siteBregion2centertriangle,1)

    #siteBregion3
    siteBregion3 = (c3)*siteBregion2
    siteBregion3closetriangle1 =  c3*siteBregion2closetriangle1
    siteBregion3closetriangle1fock = quantize(siteBregion3closetriangle1,1)
    visualize(siteBregion3closetriangle1+siteBregion3bdytriangle1+siteBregion3bdytriangle2+siteBregion3centertriangle)

    siteBregion3bdytriangle1 = c3*siteBregion2bdytriangle1
    siteBregion3bdytriangle2 = c3*siteBregion2bdytriangle2
    siteBregion3bdytriangle1fock = quantize(siteBregion3bdytriangle1,1)
    siteBregion3bdytriangle2fock = quantize(siteBregion3bdytriangle2,1)

    siteBregion3centertriangle = c3*siteBregion2centertriangle
    siteBregion3centertrianglefock = quantize(siteBregion3centertriangle,1)

    localcorrelations = regioncorrelations(blockedcorrelations,firsthexagonalregionfock)
    localspectrum = localcorrelations|>eigspec
    evec = localspectrum|>geteigenvectors
    modeevalspair = localspectrum|>geteigenvalues
    sum([-norm(eval)*log(norm(eval))-(1-norm(eval))*log(1-norm(eval)) for (mode,eval) in modeevalspair if ((norm(eval)>0.00000001) & (norm(eval)<1-0.00000001))])

    result = []
    for rsmode in firsthexagonalregionfock|>orderedmodes
        save = []
        for (mode,eval) in modeevalspair
            if ((norm(eval)>0.00000001) & (norm(eval)<1-0.00000001))
                append!(save,norm(evec[rsmode,mode])^2*(-norm(eval)*log(norm(eval))-(1-norm(eval))*log(1-norm(eval))))
            end
        end
        append!(result,tuple([rsmode,sum(save)]))
    end
    result
    ans = [(rsmode|>getattr(:b))+(rsmode|>getattr(:r)) for (rsmode,val) in sort([(rsmode,val) for (rsmode,val) in result],by=x->x[2])]
    sort([(rsmode,val) for (rsmode,val) in result],by=x->x[2])[42:54]
    test = Subset(ans[1:54])
    visualize(test)
    # localspectrum|>visualize
    courierseedsoffset = firsthexagonalregion-test
    courierseeds = quantize(courierseedsoffset,1)
    
    localstates = getregionstates(localcorrelations=localcorrelations, grouping=[27,42,27])
    localcourierisometry = localstates[2]|>FockMap
    localfilledisometry = localstates[1]|>FockMap
    localemptyisometry = localstates[3]|>FockMap
    
    localcourierseeds = idmap(firsthexagonalregionfock, firsthexagonalregionfock)[:,FockSpace(mode for mode in courierseeds)]
    localwannierization(localcourierisometry, localcourierseeds)

    ((localfilledisometry*localfilledisometry')-(localemptyisometry*localemptyisometry'))[innerhexagonalregionfock,innerhexagonalregionfock]|>eigspech|>visualize
    ((localfilledisometry*localfilledisometry')-(localemptyisometry*localemptyisometry'))[siteAregion3centertrianglefock,siteAregion3centertrianglefock]|>eigspech|>visualize

    (localcourierisometry*localcourierisometry')[siteAregion3centertrianglefock,siteAregion3centertrianglefock]|>eigspech|>visualize

    # localfilledcorrelations = idmap(localcorrelations|>getoutspace)-(localfilledisometry*localfilledisometry')
    # localemptycorrelations = idmap(localcorrelations|>getoutspace)-(localemptyisometry*localemptyisometry')
    # localcouriercorrelations = idmap(localcorrelations|>getoutspace)-(localcourierisometry*localcourierisometry')

    localfilledcorrelations[innerhexagonalregionfock,innerhexagonalregionfock]|>eigspech|>visualize

    innerhexagonalfilledseeds = getregionstates(localcorrelations=localfilledcorrelations[innerhexagonalregionfock,innerhexagonalregionfock],grouping=[12])[1]
    centerfilledseedsB1 =  getregionstates(localcorrelations=localfilledcorrelations[siteBregion1centertrianglefock,siteBregion1centertrianglefock],grouping=[2])[1]
    centerfilledseedsB2 =  getregionstates(localcorrelations=localfilledcorrelations[siteBregion2centertrianglefock,siteBregion2centertrianglefock],grouping=[2])[1]
    centerfilledseedsB3 =  getregionstates(localcorrelations=localfilledcorrelations[siteBregion3centertrianglefock,siteBregion3centertrianglefock],grouping=[2])[1]
    filledseeds = innerhexagonalfilledseeds+centerfilledseedsB1+centerfilledseedsB2+centerfilledseedsB3
    filledseeds = filledseeds|>FockMap
    filledextension = extensionmap(filledseeds|>getoutspace,localfilledisometry|>getoutspace)
    extendedfilledseeds = filledextension'*filledseeds

    localwannierfilled = localwannierization(localfilledisometry, extendedfilledseeds)

    innerhexagonalemptyseeds = getregionstates(localcorrelations=localemptycorrelations[innerhexagonalregionfock,innerhexagonalregionfock],grouping=[12])[1]
    centeremptyseedsA1 =  getregionstates(localcorrelations=localemptycorrelations[siteAregion1centertrianglefock,siteAregion1centertrianglefock],grouping=[2])[1]
    centeremptyseedsA2 =  getregionstates(localcorrelations=localemptycorrelations[siteAregion2centertrianglefock,siteAregion2centertrianglefock],grouping=[2])[1]
    centeremptyseedsA3 =  getregionstates(localcorrelations=localemptycorrelations[siteAregion3centertrianglefock,siteAregion3centertrianglefock],grouping=[2])[1]
    emptyseeds = innerhexagonalemptyseeds+centeremptyseedsA1+centeremptyseedsA2+centeremptyseedsA3
    emptyseeds  = emptyseeds |>FockMap
    emptyextension = extensionmap(emptyseeds|>getoutspace,localemptyisometry|>getoutspace)
    extendedemptyseeds = emptyextension'*emptyseeds

    localwannierempty = localwannierization(localemptyisometry, extendedemptyseeds)


    localcouriercorrelations[siteAregion1centertrianglefock,siteAregion1centertrianglefock]|>eigspech|>visualize
    localcouriercorrelations[siteAregion1bdytriangle1fock,siteAregion1bdytriangle1fock]|>eigspech|>visualize

    centertrianglecourierseedsA1 = getregionstates(localcorrelations=(localcouriercorrelations[siteAregion1centertrianglefock,siteAregion1centertrianglefock]), grouping=[2])[1]
    centertrianglecourierseedsA2 = getregionstates(localcorrelations=(localcouriercorrelations[siteAregion2centertrianglefock,siteAregion2centertrianglefock]), grouping=[2])[1]
    centertrianglecourierseedsA3 = getregionstates(localcorrelations=(localcouriercorrelations[siteAregion3centertrianglefock,siteAregion3centertrianglefock]), grouping=[2])[1]
    centertrianglecourierseedsB1 = getregionstates(localcorrelations=(localcouriercorrelations[siteBregion1centertrianglefock,siteBregion1centertrianglefock]), grouping=[2])[1]
    centertrianglecourierseedsB2 = getregionstates(localcorrelations=(localcouriercorrelations[siteBregion2centertrianglefock,siteBregion2centertrianglefock]), grouping=[2])[1]
    centertrianglecourierseedsB3 = getregionstates(localcorrelations=(localcouriercorrelations[siteBregion3centertrianglefock,siteBregion3centertrianglefock]), grouping=[2])[1]

    bdytriangle1courierseedsA1 = getregionstates(localcorrelations=(localcouriercorrelations[siteAregion1bdytriangle1fock,siteAregion1bdytriangle1fock]), grouping=[4])[1]
    bdytriangle2courierseedsA1 = getregionstates(localcorrelations=(localcouriercorrelations[siteAregion1bdytriangle2fock,siteAregion1bdytriangle2fock]), grouping=[4])[1]
    bdytriangle1courierseedsA2 = getregionstates(localcorrelations=(localcouriercorrelations[siteAregion2bdytriangle1fock,siteAregion2bdytriangle1fock]), grouping=[4])[1]
    bdytriangle2courierseedsA2 = getregionstates(localcorrelations=(localcouriercorrelations[siteAregion2bdytriangle2fock,siteAregion2bdytriangle2fock]), grouping=[4])[1]
    bdytriangle1courierseedsA3 = getregionstates(localcorrelations=(localcouriercorrelations[siteAregion3bdytriangle1fock,siteAregion3bdytriangle1fock]), grouping=[4])[1]
    bdytriangle2courierseedsA3 = getregionstates(localcorrelations=(localcouriercorrelations[siteAregion3bdytriangle2fock,siteAregion3bdytriangle2fock]), grouping=[4])[1]
    bdytriangle1courierseedsB1 = getregionstates(localcorrelations=(localcouriercorrelations[siteBregion1bdytriangle1fock,siteBregion1bdytriangle1fock]), grouping=[4])[1]
    bdytriangle2courierseedsB1 = getregionstates(localcorrelations=(localcouriercorrelations[siteBregion1bdytriangle2fock,siteBregion1bdytriangle2fock]), grouping=[4])[1]
    bdytriangle1courierseedsB2 = getregionstates(localcorrelations=(localcouriercorrelations[siteBregion2bdytriangle1fock,siteBregion2bdytriangle1fock]), grouping=[4])[1]
    bdytriangle2courierseedsB2 = getregionstates(localcorrelations=(localcouriercorrelations[siteBregion2bdytriangle2fock,siteBregion2bdytriangle2fock]), grouping=[4])[1]
    bdytriangle1courierseedsB3 = getregionstates(localcorrelations=(localcouriercorrelations[siteBregion3bdytriangle1fock,siteBregion3bdytriangle1fock]), grouping=[4])[1]
    bdytriangle2courierseedsB3 = getregionstates(localcorrelations=(localcouriercorrelations[siteBregion3bdytriangle2fock,siteBregion3bdytriangle2fock]), grouping=[4])[1]

    courierseeds = (centertrianglecourierseedsA1+centertrianglecourierseedsA2+centertrianglecourierseedsA3+centertrianglecourierseedsB1+centertrianglecourierseedsB2+centertrianglecourierseedsB3
                +bdytriangle1courierseedsA1+bdytriangle2courierseedsA1+bdytriangle1courierseedsA2+bdytriangle2courierseedsA2+bdytriangle1courierseedsA3+bdytriangle2courierseedsA3
                +bdytriangle1courierseedsB1+bdytriangle2courierseedsB1+bdytriangle1courierseedsB2+bdytriangle2courierseedsB2+bdytriangle1courierseedsB3+bdytriangle2courierseedsB3)|>FockMap
    courierextension = extensionmap(courierseeds|>getoutspace,localcourierisometry|>getoutspace)
    extendedcourierseeds = courierextension'*courierseeds
    localwanniercourier = localwannierization(localcourierisometry, extendedcourierseeds)

    localwannierresults =  Dict(:localwanniercourier => localwanniercourier, :localwannierfilled => localwannierfilled,
                                :localwannierempty => localwannierempty, :courierseeds => courierseeds)

    wannierinfos =  Dict(firsthexagonalregionfock=>localwannierresults)

    shiftedfirstcenter1 = [1,0] ∈ blockedspace
    shiftedfirstcenter2 = [0,1] ∈ blockedspace
    shiftedfirstcenter3 = [1,1] ∈ blockedspace

    firstshiftedhexagonalcenters = [shiftedfirstcenter1, shiftedfirstcenter2, shiftedfirstcenter3]