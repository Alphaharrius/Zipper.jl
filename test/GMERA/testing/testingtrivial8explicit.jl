using Zipper
using LinearAlgebra
using Plots

setmaxthreads(8)

systemsize=32

onsite = -0.1
    correlations,H = generatesystem(onsite+0im, -onsite+0im,-1 + 0im,0im,systemsize)
    crystalfock = correlations|>getoutspace

    scale = Scale([4 0; 0 4], crystalfock|>getcrystal|>getspace)
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
    siteAregion1closetriangle =  intersect(siteAregion1,innerhexagonalregion)
    siteAregion1closetrianglefock = quantize(siteAregion1closetriangle,1)

    siteAregion1bdytriangle1 = siteAregion1closetriangle.+(1/2)*rgshiftedcenter1
    siteAregion1bdytriangle2 = siteAregion1closetriangle.+(1/2)*rgshiftedcenter2
    siteAregion1bdytriangle1fock = quantize(siteAregion1bdytriangle1,1)
    siteAregion1bdytriangle2fock = quantize(siteAregion1bdytriangle2,1)

    siteAregion1centertriangle = siteAregion1 - siteAregion1closetriangle - siteAregion1bdytriangle1 - siteAregion1bdytriangle2
    siteAregion1centertrianglefock = quantize(siteAregion1centertriangle,1)

    #siteAregion2
    siteAregion2 = c3*siteAregion1
    siteAregion2closetriangle =  c3*siteAregion1closetriangle
    siteAregion2closetrianglefock = quantize(siteAregion2closetriangle,1)

    siteAregion2bdytriangle1 = c3*siteAregion1bdytriangle1
    siteAregion2bdytriangle2 = c3*siteAregion1bdytriangle2
    siteAregion2bdytriangle1fock = quantize(siteAregion2bdytriangle1,1)
    siteAregion2bdytriangle2fock = quantize(siteAregion2bdytriangle2,1)

    siteAregion2centertriangle = c3*siteAregion1centertriangle
    siteAregion2centertrianglefock = quantize(siteAregion2centertriangle,1)

    #siteAregion3
    siteAregion3 = c3*siteAregion2
    siteAregion3closetriangle =  c3*siteAregion2closetriangle
    siteAregion3closetrianglefock = quantize(siteAregion3closetriangle,1)

    siteAregion3bdytriangle1 = c3*siteAregion2bdytriangle1
    siteAregion3bdytriangle2 = c3*siteAregion2bdytriangle2
    siteAregion3bdytriangle1fock = quantize(siteAregion3bdytriangle1,1)
    siteAregion3bdytriangle2fock = quantize(siteAregion3bdytriangle2,1)

    siteAregion3centertriangle = c3*siteAregion2centertriangle
    siteAregion3centertrianglefock = quantize(siteAregion3centertriangle,1)

    #siteBregion1
    siteBregion1 = c6*siteAregion1
    siteBregion1closetriangle =  c6*siteAregion1closetriangle
    siteBregion1closetrianglefock = quantize(siteBregion1closetriangle,1)

    siteBregion1bdytriangle1 = c6*siteAregion1bdytriangle1
    siteBregion1bdytriangle2 = c6*siteAregion1bdytriangle2
    siteBregion1bdytriangle1fock = quantize(siteBregion1bdytriangle1,1)
    siteBregion1bdytriangle2fock = quantize(siteBregion1bdytriangle2,1)

    siteBregion1centertriangle = c6*siteAregion1centertriangle
    siteBregion1centertrianglefock = quantize(siteBregion1centertriangle,1)

    (siteAregion1closetriangle+siteAregion2closetriangle+siteAregion3closetriangle+siteBregion1closetriangle+siteBregion2closetriangle+siteBregion3closetriangle)|>visualize

    #siteBregion2
    siteBregion2 = c3*siteBregion1
    siteBregion2closetriangle =  c3*siteBregion1closetriangle
    siteBregion2closetrianglefock = quantize(siteBregion2closetriangle,1)

    siteBregion2bdytriangle1 = c3*siteBregion1bdytriangle1
    siteBregion2bdytriangle2 = c3*siteBregion1bdytriangle2
    siteBregion2bdytriangle1fock = quantize(siteBregion2bdytriangle1,1)
    siteBregion2bdytriangle2fock = quantize(siteBregion2bdytriangle2,1)

    siteBregion2centertriangle = c3*siteBregion1centertriangle
    siteBregion2centertrianglefock = quantize(siteBregion2centertriangle,1)

    #siteBregion3
    siteBregion3 = (c3)*siteBregion2
    siteBregion3closetriangle =  c3*siteBregion2closetriangle
    siteBregion3closetrianglefock = quantize(siteBregion3closetriangle,1)

    siteBregion3bdytriangle1 = c3*siteBregion2bdytriangle1
    siteBregion3bdytriangle2 = c3*siteBregion2bdytriangle2
    siteBregion3bdytriangle1fock = quantize(siteBregion3bdytriangle1,1)
    siteBregion3bdytriangle2fock = quantize(siteBregion3bdytriangle2,1)

    siteBregion3centertriangle = c3*siteBregion2centertriangle
    siteBregion3centertrianglefock = quantize(siteBregion3centertriangle,1)

    localcorrelations = regioncorrelations(blockedcorrelations,firsthexagonalregionfock)
    localspectrum = localcorrelations|>eigspec
    localspectrum|>visualize

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
    sort([(rsmode,val) for (rsmode,val) in result],by=x->x[2])[174:181]
    test = Subset(ans[1:120])
    visualize(test)

    # visualize(firsthexagonalregion)
    innerhexagonalregion = rankedandgroupoffsets(firsthexagonalregion,11)

    courierseedsoffset = firsthexagonalregion-test
    visualize(courierseedsoffset)
    visualize(c6*courierseedsoffset)
    visualize(c3*courierseedsoffset)
    courierseeds = quantize(courierseedsoffset,1)
    
    localstates = getregionstates(localcorrelations=localcorrelations, grouping=[60,264,60])
    localcourierisometry = localstates[2]|>FockMap
    # localfilledisometry = localstates[1]|>FockMap
    # localemptyisometry = localstates[3]|>FockMap

    localcourierseeds = idmap(firsthexagonalregionfock, firsthexagonalregionfock)[:,FockSpace(mode for mode in courierseeds)]
    localwannierization(localcourierisometry, localcourierseeds)

    localfilledcorrelations = idmap(localcorrelations|>getoutspace)-(localfilledisometry*localfilledisometry')
    localemptycorrelations = idmap(localcorrelations|>getoutspace)-(localemptyisometry*localemptyisometry')
    localcouriercorrelations = idmap(localcorrelations|>getoutspace)-(localcourierisometry*localcourierisometry')

    # localfilledcorrelations[siteBregion3closetrianglefock,siteBregion3closetrianglefock]|>eigspech|>visualize
    # localcouriercorrelations[siteAregion1bdytriangle1fock,siteAregion1bdytriangle1fock]|>eigspech|>visualize

    # subregion = rankedandgroupoffsets(firsthexagonalregion,10)
    # subregionfock = quantize(subregion,1)
    # # subregion|>visualize
    # localfilledcorrelations[subregionfock,subregionfock]|>eigspech|>visualize

    # localfilledcorrelations[siteAregion1centertrianglefock,siteAregion1centertrianglefock]|>eigspech|>visualize
    # localfilledcorrelations[innerhexagonalregionfock,innerhexagonalregionfock]|>eigspech|>visualize
    # localfilledcorrelations[siteAregion1bdytriangle1fock,siteAregion1bdytriangle1fock]|>eigspech|>visualize
    # localfilledcorrelations[siteAregion1closetrianglefock,siteAregion1closetrianglefock]|>eigspech|>visualize

    # innerhexagonalfilledseeds = getregionstates(localcorrelations=localfilledcorrelations[innerhexagonalregionfock,innerhexagonalregionfock],grouping=[36])[1]
    # innerhexagonalfilledseeds = getregionstates(localcorrelations=localfilledcorrelations[subregionfock,subregionfock],grouping=[48])[1]

    # siteAregion1closetrianglefilledseeds = getregionstates(localcorrelations=localfilledcorrelations[siteAregion1closetrianglefock,siteAregion1closetrianglefock],grouping=[6])[1]
    # siteAregion2closetrianglefilledseeds = getregionstates(localcorrelations=localfilledcorrelations[siteAregion2closetrianglefock,siteAregion2closetrianglefock],grouping=[6])[1]
    # siteAregion3closetrianglefilledseeds = getregionstates(localcorrelations=localfilledcorrelations[siteAregion3closetrianglefock,siteAregion3closetrianglefock],grouping=[6])[1]
    # siteBregion1closetrianglefilledseeds = getregionstates(localcorrelations=localfilledcorrelations[siteBregion1closetrianglefock,siteBregion1closetrianglefock],grouping=[6])[1]
    # siteBregion2closetrianglefilledseeds = getregionstates(localcorrelations=localfilledcorrelations[siteBregion2closetrianglefock,siteBregion2closetrianglefock],grouping=[6])[1]
    # siteBregion3closetrianglefilledseeds = getregionstates(localcorrelations=localfilledcorrelations[siteBregion3closetrianglefock,siteBregion3closetrianglefock],grouping=[6])[1]

    # centerfilledseedsA1 =  getregionstates(localcorrelations=localfilledcorrelations[siteAregion1centertrianglefock,siteAregion1centertrianglefock],grouping=[4])[1]
    # centerfilledseedsA2 =  getregionstates(localcorrelations=localfilledcorrelations[siteAregion2centertrianglefock,siteAregion2centertrianglefock],grouping=[4])[1]
    # centerfilledseedsA3 =  getregionstates(localcorrelations=localfilledcorrelations[siteAregion3centertrianglefock,siteAregion3centertrianglefock],grouping=[4])[1]
    # centerfilledseedsB1 =  getregionstates(localcorrelations=localfilledcorrelations[siteBregion1centertrianglefock,siteBregion1centertrianglefock],grouping=[4])[1]
    # centerfilledseedsB2 =  getregionstates(localcorrelations=localfilledcorrelations[siteBregion2centertrianglefock,siteBregion2centertrianglefock],grouping=[4])[1]
    # centerfilledseedsB3 =  getregionstates(localcorrelations=localfilledcorrelations[siteBregion3centertrianglefock,siteBregion3centertrianglefock],grouping=[4])[1]

    # bdytriangle1filledseedsA1 = getregionstates(localcorrelations=(localfilledcorrelations[siteAregion1bdytriangle1fock,siteAregion1bdytriangle1fock]), grouping=[1])[1]
    # bdytriangle2filledseedsA1 = getregionstates(localcorrelations=(localfilledcorrelations[siteAregion1bdytriangle2fock,siteAregion1bdytriangle2fock]), grouping=[1])[1]
    # bdytriangle1filledseedsA2 = getregionstates(localcorrelations=(localfilledcorrelations[siteAregion2bdytriangle1fock,siteAregion2bdytriangle1fock]), grouping=[1])[1]
    # bdytriangle2filledseedsA2 = getregionstates(localcorrelations=(localfilledcorrelations[siteAregion2bdytriangle2fock,siteAregion2bdytriangle2fock]), grouping=[1])[1]
    # bdytriangle1filledseedsA3 = getregionstates(localcorrelations=(localfilledcorrelations[siteAregion3bdytriangle1fock,siteAregion3bdytriangle1fock]), grouping=[1])[1]
    # bdytriangle2filledseedsA3 = getregionstates(localcorrelations=(localfilledcorrelations[siteAregion3bdytriangle2fock,siteAregion3bdytriangle2fock]), grouping=[1])[1]
    # bdytriangle1filledseedsB1 = getregionstates(localcorrelations=(localfilledcorrelations[siteBregion1bdytriangle1fock,siteBregion1bdytriangle1fock]), grouping=[1])[1]
    # bdytriangle2filledseedsB1 = getregionstates(localcorrelations=(localfilledcorrelations[siteBregion1bdytriangle2fock,siteBregion1bdytriangle2fock]), grouping=[1])[1]
    # bdytriangle1filledseedsB2 = getregionstates(localcorrelations=(localfilledcorrelations[siteBregion2bdytriangle1fock,siteBregion2bdytriangle1fock]), grouping=[1])[1]
    # bdytriangle2filledseedsB2 = getregionstates(localcorrelations=(localfilledcorrelations[siteBregion2bdytriangle2fock,siteBregion2bdytriangle2fock]), grouping=[1])[1]
    # bdytriangle1filledseedsB3 = getregionstates(localcorrelations=(localfilledcorrelations[siteBregion3bdytriangle1fock,siteBregion3bdytriangle1fock]), grouping=[1])[1]
    # bdytriangle2filledseedsB3 = getregionstates(localcorrelations=(localfilledcorrelations[siteBregion3bdytriangle2fock,siteBregion3bdytriangle2fock]), grouping=[1])[1]

    # filledseeds = (innerhexagonalfilledseeds
    #             +centerfilledseedsA1+centerfilledseedsA2+centerfilledseedsA3+centerfilledseedsB1+centerfilledseedsB2+centerfilledseedsB3
    #             +bdytriangle1filledseedsA1+bdytriangle2filledseedsA1+bdytriangle1filledseedsA2+bdytriangle2filledseedsA2+bdytriangle1filledseedsA3+bdytriangle2filledseedsA3
    #             +bdytriangle1filledseedsB1+bdytriangle2filledseedsB1+bdytriangle1filledseedsB2+bdytriangle2filledseedsB2+bdytriangle1filledseedsB3+bdytriangle2filledseedsB3)|>FockMap
    # # filledseeds = (bdytriangle1filledseedsA1+bdytriangle2filledseedsA1+bdytriangle1filledseedsA2+bdytriangle2filledseedsA2+bdytriangle1filledseedsA3+bdytriangle2filledseedsA3
    # # +bdytriangle1filledseedsB1+bdytriangle2filledseedsB1+bdytriangle1filledseedsB2+bdytriangle2filledseedsB2+bdytriangle1filledseedsB3+bdytriangle2filledseedsB3)|>FockMap
    # # filledseeds = (siteAregion1closetrianglefilledseeds+siteAregion2closetrianglefilledseeds+siteAregion3closetrianglefilledseeds
    # # +siteBregion1closetrianglefilledseeds+siteBregion2closetrianglefilledseeds+siteBregion3closetrianglefilledseeds)|>FockMap
    # filledextension = extensionmap(filledseeds|>getoutspace,localfilledisometry|>getoutspace)
    # extendedfilledseeds = filledextension'*filledseeds

    # localwannierfilled = localwannierization(localfilledisometry, extendedfilledseeds)

    # innerhexagonalemptyseeds = getregionstates(localcorrelations=localemptycorrelations[innerhexagonalregionfock,innerhexagonalregionfock],grouping=[12])[1]
    # centeremptyseedsA1 =  getregionstates(localcorrelations=localemptycorrelations[siteAregion1centertrianglefock,siteAregion1centertrianglefock],grouping=[2])[1]
    # centeremptyseedsA2 =  getregionstates(localcorrelations=localemptycorrelations[siteAregion2centertrianglefock,siteAregion2centertrianglefock],grouping=[2])[1]
    # centeremptyseedsA3 =  getregionstates(localcorrelations=localemptycorrelations[siteAregion3centertrianglefock,siteAregion3centertrianglefock],grouping=[2])[1]
    # emptyseeds = innerhexagonalemptyseeds+centeremptyseedsA1+centeremptyseedsA2+centeremptyseedsA3
    # emptyseeds  = emptyseeds |>FockMap
    # emptyextension = extensionmap(emptyseeds|>getoutspace,localemptyisometry|>getoutspace)
    # extendedemptyseeds = emptyextension'*emptyseeds

    # localwannierempty = localwannierization(localemptyisometry, extendedemptyseeds)

    # wholeA1 = (siteAregion1bdytriangle1fock+siteAregion1bdytriangle2fock+siteAregion1closetrianglefock+siteAregion1centertrianglefock)|>RegionFock
    # wholeA2 = (siteAregion2bdytriangle1fock+siteAregion2bdytriangle2fock+siteAregion2closetrianglefock+siteAregion2centertrianglefock)|>RegionFock
    # wholeA3 = (siteAregion3bdytriangle1fock+siteAregion3bdytriangle2fock+siteAregion3closetrianglefock+siteAregion3centertrianglefock)|>RegionFock

    # wholeB1 = (siteBregion1bdytriangle1fock+siteBregion1bdytriangle2fock+siteBregion1closetrianglefock+siteBregion1centertrianglefock)|>RegionFock
    # wholeB2 = (siteBregion2bdytriangle1fock+siteBregion2bdytriangle2fock+siteBregion2closetrianglefock+siteBregion2centertrianglefock)|>RegionFock
    # wholeB3 = (siteBregion3bdytriangle1fock+siteBregion3bdytriangle2fock+siteBregion3closetrianglefock+siteBregion3centertrianglefock)|>RegionFock

    # wholeA1seeds = getregionstates(localcorrelations=localcouriercorrelations[wholeA1,wholeA1],grouping=[40])[1]
    # wholeA2seeds = getregionstates(localcorrelations=localcouriercorrelations[wholeA2,wholeA2],grouping=[40])[1]
    # wholeA3seeds = getregionstates(localcorrelations=localcouriercorrelations[wholeA3,wholeA3],grouping=[40])[1]

    # wholeB1seeds = getregionstates(localcorrelations=localcouriercorrelations[wholeB1,wholeB1],grouping=[40])[1]
    # wholeB2seeds = getregionstates(localcorrelations=localcouriercorrelations[wholeB2,wholeB2],grouping=[40])[1]
    # wholeB3seeds = getregionstates(localcorrelations=localcouriercorrelations[wholeB3,wholeB3],grouping=[40])[1]

    # localcouriercorrelations[wholeA1,wholeA1]|>eigspech|>visualize
    # localcouriercorrelations[siteAregion2bdytriangle1fock,siteAregion2bdytriangle1fock]|>eigspech|>visualize
    # localcouriercorrelations[siteAregion2centertrianglefock,siteAregion2centertrianglefock]|>eigspech|>visualize
    # localcouriercorrelations[innerhexagonalregionfock,innerhexagonalregionfock]|>eigspech|>visualize


    siteAregion1closetrianglecourierseeds = getregionstates(localcorrelations=localcouriercorrelations[siteAregion1closetrianglefock,siteAregion1closetrianglefock],grouping=[4])[1]
    siteAregion2closetrianglecourierseeds = getregionstates(localcorrelations=localcouriercorrelations[siteAregion2closetrianglefock,siteAregion2closetrianglefock],grouping=[4])[1]
    siteAregion3closetrianglecourierseeds = getregionstates(localcorrelations=localcouriercorrelations[siteAregion3closetrianglefock,siteAregion3closetrianglefock],grouping=[4])[1]
    siteBregion1closetrianglecourierseeds = getregionstates(localcorrelations=localcouriercorrelations[siteBregion1closetrianglefock,siteBregion1closetrianglefock],grouping=[4])[1]
    siteBregion2closetrianglecourierseeds = getregionstates(localcorrelations=localcouriercorrelations[siteBregion2closetrianglefock,siteBregion2closetrianglefock],grouping=[4])[1]
    siteBregion3closetrianglecourierseeds = getregionstates(localcorrelations=localcouriercorrelations[siteBregion3closetrianglefock,siteBregion3closetrianglefock],grouping=[4])[1]

    centertrianglecourierseedsA1 = getregionstates(localcorrelations=(localcouriercorrelations[siteAregion1centertrianglefock,siteAregion1centertrianglefock]), grouping=[8])[1]
    centertrianglecourierseedsA2 = getregionstates(localcorrelations=(localcouriercorrelations[siteAregion2centertrianglefock,siteAregion2centertrianglefock]), grouping=[8])[1]
    centertrianglecourierseedsA3 = getregionstates(localcorrelations=(localcouriercorrelations[siteAregion3centertrianglefock,siteAregion3centertrianglefock]), grouping=[8])[1]
    centertrianglecourierseedsB1 = getregionstates(localcorrelations=(localcouriercorrelations[siteBregion1centertrianglefock,siteBregion1centertrianglefock]), grouping=[8])[1]
    centertrianglecourierseedsB2 = getregionstates(localcorrelations=(localcouriercorrelations[siteBregion2centertrianglefock,siteBregion2centertrianglefock]), grouping=[8])[1]
    centertrianglecourierseedsB3 = getregionstates(localcorrelations=(localcouriercorrelations[siteBregion3centertrianglefock,siteBregion3centertrianglefock]), grouping=[8])[1]

    bdytriangle1courierseedsA1 = getregionstates(localcorrelations=(localcouriercorrelations[siteAregion1bdytriangle1fock,siteAregion1bdytriangle1fock]), grouping=[14])[1]
    bdytriangle2courierseedsA1 = getregionstates(localcorrelations=(localcouriercorrelations[siteAregion1bdytriangle2fock,siteAregion1bdytriangle2fock]), grouping=[14])[1]
    bdytriangle1courierseedsA2 = getregionstates(localcorrelations=(localcouriercorrelations[siteAregion2bdytriangle1fock,siteAregion2bdytriangle1fock]), grouping=[14])[1]
    bdytriangle2courierseedsA2 = getregionstates(localcorrelations=(localcouriercorrelations[siteAregion2bdytriangle2fock,siteAregion2bdytriangle2fock]), grouping=[14])[1]
    bdytriangle1courierseedsA3 = getregionstates(localcorrelations=(localcouriercorrelations[siteAregion3bdytriangle1fock,siteAregion3bdytriangle1fock]), grouping=[14])[1]
    bdytriangle2courierseedsA3 = getregionstates(localcorrelations=(localcouriercorrelations[siteAregion3bdytriangle2fock,siteAregion3bdytriangle2fock]), grouping=[14])[1]
    bdytriangle1courierseedsB1 = getregionstates(localcorrelations=(localcouriercorrelations[siteBregion1bdytriangle1fock,siteBregion1bdytriangle1fock]), grouping=[14])[1]
    bdytriangle2courierseedsB1 = getregionstates(localcorrelations=(localcouriercorrelations[siteBregion1bdytriangle2fock,siteBregion1bdytriangle2fock]), grouping=[14])[1]
    bdytriangle1courierseedsB2 = getregionstates(localcorrelations=(localcouriercorrelations[siteBregion2bdytriangle1fock,siteBregion2bdytriangle1fock]), grouping=[14])[1]
    bdytriangle2courierseedsB2 = getregionstates(localcorrelations=(localcouriercorrelations[siteBregion2bdytriangle2fock,siteBregion2bdytriangle2fock]), grouping=[14])[1]
    bdytriangle1courierseedsB3 = getregionstates(localcorrelations=(localcouriercorrelations[siteBregion3bdytriangle1fock,siteBregion3bdytriangle1fock]), grouping=[14])[1]
    bdytriangle2courierseedsB3 = getregionstates(localcorrelations=(localcouriercorrelations[siteBregion3bdytriangle2fock,siteBregion3bdytriangle2fock]), grouping=[14])[1]

    # courierseeds = (centertrianglecourierseedsA1+centertrianglecourierseedsA2+centertrianglecourierseedsA3+centertrianglecourierseedsB1+centertrianglecourierseedsB2+centertrianglecourierseedsB3
    #     +bdytriangle1courierseedsA1+bdytriangle2courierseedsA1+bdytriangle1courierseedsA2+bdytriangle2courierseedsA2+bdytriangle1courierseedsA3+bdytriangle2courierseedsA3
    #     +bdytriangle1courierseedsB1+bdytriangle2courierseedsB1+bdytriangle1courierseedsB2+bdytriangle2courierseedsB2+bdytriangle1courierseedsB3+bdytriangle2courierseedsB3)|>FockMap
    courierseeds = (siteAregion1closetrianglecourierseeds+siteAregion2closetrianglecourierseeds+siteAregion3closetrianglecourierseeds+siteBregion1closetrianglecourierseeds+siteBregion2closetrianglecourierseeds+siteBregion3closetrianglecourierseeds
        +centertrianglecourierseedsA1+centertrianglecourierseedsA2+centertrianglecourierseedsA3+centertrianglecourierseedsB1+centertrianglecourierseedsB2+centertrianglecourierseedsB3
        +bdytriangle1courierseedsA1+bdytriangle2courierseedsA1+bdytriangle1courierseedsA2+bdytriangle2courierseedsA2+bdytriangle1courierseedsA3+bdytriangle2courierseedsA3
        +bdytriangle1courierseedsB1+bdytriangle2courierseedsB1+bdytriangle1courierseedsB2+bdytriangle2courierseedsB2+bdytriangle1courierseedsB3+bdytriangle2courierseedsB3)|>FockMap
    # courierextension = extensionmap(courierseeds|>getoutspace,localcourierisometry|>getoutspace)
    # extendedcourierseeds = courierextension'*courierseeds
    # localwanniercourier = localwannierization(localcourierisometry, extendedcourierseeds)
    # courierseeds = (wholeA1seeds+wholeA2seeds+wholeA3seeds+wholeB1seeds+wholeB2seeds+wholeB3seeds)|>FockMap
    localwanniercourier,svd = localwannierization(localcourierisometry, courierseeds)

    localwannierresults =  Dict(:localwanniercourier => localwanniercourier, :localwannierfilled => localfilledisometry,
                                :localwannierempty => localemptyisometry, :courierseeds => courierseeds)

    wannierinfos =  Dict(firsthexagonalregionfock=>localwannierresults)

    shiftedfirstcenter1 = [1,0] ∈ blockedspace
    shiftedfirstcenter2 = [0,1] ∈ blockedspace
    shiftedfirstcenter3 = [1,1] ∈ blockedspace

    firstshiftedhexagonalcenters = [shiftedfirstcenter1, shiftedfirstcenter2, shiftedfirstcenter3]

    for hexagonalcenter in firstshiftedhexagonalcenters
        shiftedfirsthexagonalregion = firsthexagonalregion.+hexagonalcenter
        shiftedfirsthexagonalregionfock = quantize(shiftedfirsthexagonalregion,1)
        shiftedinnerhexagonalregion = innerhexagonalregion.+hexagonalcenter
        shiftedinnerhexagonalregionfock = quantize(shiftedinnerhexagonalregion,1)

        c6recenter = recenter(c6,hexagonalcenter)
        c3recenter = recenter(c3,hexagonalcenter)

        #siteAregion1
        shiftedsiteAregion1 = siteAregion1.+hexagonalcenter
        shiftedsiteAregion1closetriangle =  intersect(shiftedsiteAregion1,shiftedinnerhexagonalregion)
        shiftedsiteAregion1closetrianglefock = quantize(shiftedsiteAregion1closetriangle,1)

        shiftedsiteAregion1bdytriangle1 = shiftedsiteAregion1closetriangle.+(1/2)*rgshiftedcenter1
        shiftedsiteAregion1bdytriangle2 = shiftedsiteAregion1closetriangle.+(1/2)*rgshiftedcenter2
        shiftedsiteAregion1bdytriangle1fock = quantize(shiftedsiteAregion1bdytriangle1,1)
        shiftedsiteAregion1bdytriangle2fock = quantize(shiftedsiteAregion1bdytriangle2,1)

        shiftedsiteAregion1centertriangle = shiftedsiteAregion1 - shiftedsiteAregion1closetriangle - shiftedsiteAregion1bdytriangle1 - shiftedsiteAregion1bdytriangle2
        shiftedsiteAregion1centertrianglefock = quantize(shiftedsiteAregion1centertriangle,1)

        #siteAregion2
        shiftedsiteAregion2 = c3recenter*shiftedsiteAregion1
        shiftedsiteAregion2closetriangle =  c3recenter*shiftedsiteAregion1closetriangle
        shiftedsiteAregion2closetrianglefock = quantize(shiftedsiteAregion2closetriangle,1)

        shiftedsiteAregion2bdytriangle1 = c3recenter*shiftedsiteAregion1bdytriangle1
        shiftedsiteAregion2bdytriangle2 = c3recenter*shiftedsiteAregion1bdytriangle2
        shiftedsiteAregion2bdytriangle1fock = quantize(shiftedsiteAregion2bdytriangle1,1)
        shiftedsiteAregion2bdytriangle2fock = quantize(shiftedsiteAregion2bdytriangle2,1)

        shiftedsiteAregion2centertriangle = c3recenter*shiftedsiteAregion1centertriangle
        shiftedsiteAregion2centertrianglefock = quantize(shiftedsiteAregion2centertriangle,1)

        #siteAregion3
        shiftedsiteAregion3 = c3recenter*shiftedsiteAregion2
        shiftedsiteAregion3closetriangle =  c3recenter*shiftedsiteAregion2closetriangle
        shiftedsiteAregion3closetrianglefock = quantize(shiftedsiteAregion3closetriangle,1)

        shiftedsiteAregion3bdytriangle1 = c3recenter*shiftedsiteAregion2bdytriangle1
        shiftedsiteAregion3bdytriangle2 = c3recenter*shiftedsiteAregion2bdytriangle2
        shiftedsiteAregion3bdytriangle1fock = quantize(shiftedsiteAregion3bdytriangle1,1)
        shiftedsiteAregion3bdytriangle2fock = quantize(shiftedsiteAregion3bdytriangle2,1)

        shiftedsiteAregion3centertriangle = c3recenter*shiftedsiteAregion2centertriangle
        shiftedsiteAregion3centertrianglefock = quantize(shiftedsiteAregion3centertriangle,1)

        #siteBregion1
        shiftedsiteBregion1 = c6recenter*shiftedsiteAregion1
        shiftedsiteBregion1closetriangle =  c6recenter*shiftedsiteAregion1closetriangle
        shiftedsiteBregion1closetrianglefock = quantize(shiftedsiteBregion1closetriangle,1)

        shiftedsiteBregion1bdytriangle1 = c6recenter*shiftedsiteAregion1bdytriangle1
        shiftedsiteBregion1bdytriangle2 = c6recenter*shiftedsiteAregion1bdytriangle2
        shiftedsiteBregion1bdytriangle1fock = quantize(shiftedsiteBregion1bdytriangle1,1)
        shiftedsiteBregion1bdytriangle2fock = quantize(shiftedsiteBregion1bdytriangle2,1)

        shiftedsiteBregion1centertriangle = c6recenter*shiftedsiteAregion1centertriangle
        shiftedsiteBregion1centertrianglefock = quantize(shiftedsiteBregion1centertriangle,1)

        #siteBregion2
        shiftedsiteBregion2 = c3recenter*shiftedsiteBregion1
        shiftedsiteBregion2closetriangle =  c3recenter*shiftedsiteBregion1closetriangle
        shiftedsiteBregion2closetrianglefock = quantize(shiftedsiteBregion2closetriangle,1)

        shiftedsiteBregion2bdytriangle1 = c3recenter*shiftedsiteBregion1bdytriangle1
        shiftedsiteBregion2bdytriangle2 = c3recenter*shiftedsiteBregion1bdytriangle2
        shiftedsiteBregion2bdytriangle1fock = quantize(shiftedsiteBregion2bdytriangle1,1)
        shiftedsiteBregion2bdytriangle2fock = quantize(shiftedsiteBregion2bdytriangle2,1)

        shiftedsiteBregion2centertriangle = c3recenter*shiftedsiteBregion1centertriangle
        shiftedsiteBregion2centertrianglefock = quantize(shiftedsiteBregion2centertriangle,1)

        #siteBregion3
        shiftedsiteBregion3 = (c3recenter)*shiftedsiteBregion2
        shiftedsiteBregion3closetriangle =  c3recenter*shiftedsiteBregion2closetriangle
        shiftedsiteBregion3closetrianglefock = quantize(shiftedsiteBregion3closetriangle,1)

        shiftedsiteBregion3bdytriangle1 = c3recenter*shiftedsiteBregion2bdytriangle1
        shiftedsiteBregion3bdytriangle2 = c3recenter*shiftedsiteBregion2bdytriangle2
        shiftedsiteBregion3bdytriangle1fock = quantize(shiftedsiteBregion3bdytriangle1,1)
        shiftedsiteBregion3bdytriangle2fock = quantize(shiftedsiteBregion3bdytriangle2,1)

        shiftedsiteBregion3centertriangle = c3recenter*shiftedsiteBregion2centertriangle
        shiftedsiteBregion3centertrianglefock = quantize(shiftedsiteBregion3centertriangle,1)
    
        localcorrelations = regioncorrelations(blockedcorrelations,shiftedfirsthexagonalregionfock)
        localspectrum = localcorrelations|>eigspec
        # localspectrum|>visualize
    
        localstates = getregionstates(localcorrelations=localcorrelations, grouping=[72,240,72])
        localcourierisometry = localstates[2]|>FockMap
        localfilledisometry = localstates[1]|>FockMap
        localemptyisometry = localstates[3]|>FockMap

        localfilledcorrelations = idmap(localcorrelations|>getoutspace)-(localfilledisometry*localfilledisometry')
        localemptycorrelations = idmap(localcorrelations|>getoutspace)-(localemptyisometry*localemptyisometry')
        localcouriercorrelations = idmap(localcorrelations|>getoutspace)-(localcourierisometry*localcourierisometry')
    
    
        siteAregion1closetrianglecourierseeds = getregionstates(localcorrelations=localcouriercorrelations[shiftedsiteAregion1closetrianglefock,shiftedsiteAregion1closetrianglefock],grouping=[4])[1]
        siteAregion2closetrianglecourierseeds = getregionstates(localcorrelations=localcouriercorrelations[shiftedsiteAregion2closetrianglefock,shiftedsiteAregion2closetrianglefock],grouping=[4])[1]
        siteAregion3closetrianglecourierseeds = getregionstates(localcorrelations=localcouriercorrelations[shiftedsiteAregion3closetrianglefock,shiftedsiteAregion3closetrianglefock],grouping=[4])[1]
        siteBregion1closetrianglecourierseeds = getregionstates(localcorrelations=localcouriercorrelations[shiftedsiteBregion1closetrianglefock,shiftedsiteBregion1closetrianglefock],grouping=[4])[1]
        siteBregion2closetrianglecourierseeds = getregionstates(localcorrelations=localcouriercorrelations[shiftedsiteBregion2closetrianglefock,shiftedsiteBregion2closetrianglefock],grouping=[4])[1]
        siteBregion3closetrianglecourierseeds = getregionstates(localcorrelations=localcouriercorrelations[shiftedsiteBregion3closetrianglefock,shiftedsiteBregion3closetrianglefock],grouping=[4])[1]

        centertrianglecourierseedsA1 = getregionstates(localcorrelations=(localcouriercorrelations[shiftedsiteAregion1centertrianglefock,shiftedsiteAregion1centertrianglefock]), grouping=[8])[1]
        centertrianglecourierseedsA2 = getregionstates(localcorrelations=(localcouriercorrelations[shiftedsiteAregion2centertrianglefock,shiftedsiteAregion2centertrianglefock]), grouping=[8])[1]
        centertrianglecourierseedsA3 = getregionstates(localcorrelations=(localcouriercorrelations[shiftedsiteAregion3centertrianglefock,shiftedsiteAregion3centertrianglefock]), grouping=[8])[1]
        centertrianglecourierseedsB1 = getregionstates(localcorrelations=(localcouriercorrelations[shiftedsiteBregion1centertrianglefock,shiftedsiteBregion1centertrianglefock]), grouping=[8])[1]
        centertrianglecourierseedsB2 = getregionstates(localcorrelations=(localcouriercorrelations[shiftedsiteBregion2centertrianglefock,shiftedsiteBregion2centertrianglefock]), grouping=[8])[1]
        centertrianglecourierseedsB3 = getregionstates(localcorrelations=(localcouriercorrelations[shiftedsiteBregion3centertrianglefock,shiftedsiteBregion3centertrianglefock]), grouping=[8])[1]

        bdytriangle1courierseedsA1 = getregionstates(localcorrelations=(localcouriercorrelations[shiftedsiteAregion1bdytriangle1fock,shiftedsiteAregion1bdytriangle1fock]), grouping=[14])[1]
        bdytriangle2courierseedsA1 = getregionstates(localcorrelations=(localcouriercorrelations[shiftedsiteAregion1bdytriangle2fock,shiftedsiteAregion1bdytriangle2fock]), grouping=[14])[1]
        bdytriangle1courierseedsA2 = getregionstates(localcorrelations=(localcouriercorrelations[shiftedsiteAregion2bdytriangle1fock,shiftedsiteAregion2bdytriangle1fock]), grouping=[14])[1]
        bdytriangle2courierseedsA2 = getregionstates(localcorrelations=(localcouriercorrelations[shiftedsiteAregion2bdytriangle2fock,shiftedsiteAregion2bdytriangle2fock]), grouping=[14])[1]
        bdytriangle1courierseedsA3 = getregionstates(localcorrelations=(localcouriercorrelations[shiftedsiteAregion3bdytriangle1fock,shiftedsiteAregion3bdytriangle1fock]), grouping=[14])[1]
        bdytriangle2courierseedsA3 = getregionstates(localcorrelations=(localcouriercorrelations[shiftedsiteAregion3bdytriangle2fock,shiftedsiteAregion3bdytriangle2fock]), grouping=[14])[1]
        bdytriangle1courierseedsB1 = getregionstates(localcorrelations=(localcouriercorrelations[shiftedsiteBregion1bdytriangle1fock,shiftedsiteBregion1bdytriangle1fock]), grouping=[14])[1]
        bdytriangle2courierseedsB1 = getregionstates(localcorrelations=(localcouriercorrelations[shiftedsiteBregion1bdytriangle2fock,shiftedsiteBregion1bdytriangle2fock]), grouping=[14])[1]
        bdytriangle1courierseedsB2 = getregionstates(localcorrelations=(localcouriercorrelations[shiftedsiteBregion2bdytriangle1fock,shiftedsiteBregion2bdytriangle1fock]), grouping=[14])[1]
        bdytriangle2courierseedsB2 = getregionstates(localcorrelations=(localcouriercorrelations[shiftedsiteBregion2bdytriangle2fock,shiftedsiteBregion2bdytriangle2fock]), grouping=[14])[1]
        bdytriangle1courierseedsB3 = getregionstates(localcorrelations=(localcouriercorrelations[shiftedsiteBregion3bdytriangle1fock,shiftedsiteBregion3bdytriangle1fock]), grouping=[14])[1]
        bdytriangle2courierseedsB3 = getregionstates(localcorrelations=(localcouriercorrelations[shiftedsiteBregion3bdytriangle2fock,shiftedsiteBregion3bdytriangle2fock]), grouping=[14])[1]
    
        courierseeds = (siteAregion1closetrianglecourierseeds+siteAregion2closetrianglecourierseeds+siteAregion3closetrianglecourierseeds+siteBregion1closetrianglecourierseeds+siteBregion2closetrianglecourierseeds+siteBregion3closetrianglecourierseeds
            +centertrianglecourierseedsA1+centertrianglecourierseedsA2+centertrianglecourierseedsA3+centertrianglecourierseedsB1+centertrianglecourierseedsB2+centertrianglecourierseedsB3
            +bdytriangle1courierseedsA1+bdytriangle2courierseedsA1+bdytriangle1courierseedsA2+bdytriangle2courierseedsA2+bdytriangle1courierseedsA3+bdytriangle2courierseedsA3
            +bdytriangle1courierseedsB1+bdytriangle2courierseedsB1+bdytriangle1courierseedsB2+bdytriangle2courierseedsB2+bdytriangle1courierseedsB3+bdytriangle2courierseedsB3)|>FockMap
        localwanniercourier,svd = localwannierization(localcourierisometry, courierseeds)
    
    
        shiftedlocalwannierresults =  Dict(:localwanniercourier => localwanniercourier, :localwannierfilled => localfilledisometry,
                                            :localwannierempty => localemptyisometry, :courierseeds => courierseeds)
        wannierinfos[shiftedfirsthexagonalregionfock] =  shiftedlocalwannierresults
    end

    ref = [quantize(firsthexagonalregion.+hexagonalcenter,1) for hexagonalcenter in firstshiftedhexagonalcenters]
    firsthexagonalregionfocklist = [firsthexagonalregionfock,ref...]
    
    extendedwannierizedcourier =  sum(wannierinfos[regionfock][:localwanniercourier] for regionfock in firsthexagonalregionfocklist)
    extendedwannierizedfilled =  sum(wannierinfos[regionfock][:localwannierfilled] for regionfock in firsthexagonalregionfocklist)
    extendedwannierizedempty =  sum(wannierinfos[regionfock][:localwannierempty] for regionfock in firsthexagonalregionfocklist)
    
                    
    origin = [0, 0] ∈ blockedspace
    refunictcellfockcourier = FockSpace(Subset(mode for mode in extendedwannierizedcourier |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    refunictcellfockfilled = FockSpace(Subset(mode for mode in extendedwannierizedfilled |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    refunictcellfockempty = FockSpace(Subset(mode for mode in extendedwannierizedempty |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))

    wanniercourierisometry1 = globalwannierfunction(blockedcorrelations,extendedwannierizedcourier[:,refunictcellfockcourier])
    wannierfilledisometry1 = globalwannierfunction(blockedcorrelations,extendedwannierizedfilled[:,refunictcellfockfilled])
    wannieremptyisometry1 = globalwannierfunction(blockedcorrelations,extendedwannierizedempty[:,refunictcellfockempty])

    @info "Computing local courier states..."
    leftrestrict = fourier(wanniercourierisometry1|>getoutspace, firsthexagonalregionfock) / (blockedcrystal|>vol|>sqrt)
    rightrestrict = fourier(wanniercourierisometry1|>getinspace, wanniercourierisometry1|>getinspace|>unitcellfock|>RegionFock)
    wanniercourierstate = leftrestrict' * wanniercourierisometry1 * rightrestrict

    wanniercourierrsstate = wanniercourierstate|>RegionState
    wanniercourierrsstate[28]|>visualize

    couriercorrelations = wanniercourierisometry1' * blockedcorrelations * wanniercourierisometry1
    display(couriercorrelations|>eigspech|>visualize)

    filledcorrelations = wannierfilledisometry1' * blockedcorrelations * wannierfilledisometry1

    emptycorrelations = wannieremptyisometry1' * blockedcorrelations * wannieremptyisometry1

    couriercorrelationspectrum = couriercorrelations |> crystalspectrum

    purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
    purifiedcouriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap

    firstgmeracorrelations = purifiedcouriercorrelations

    firstgmeracrystalfock = firstgmeracorrelations|>getoutspace
    firstgmeracrystal::Crystal = firstgmeracrystalfock|>getcrystal
    firstgmeraspace::RealSpace = firstgmeracrystal|>getspace
    
    refrot = inv([2/3 -1/3; -1/3 2/3]')
    c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
    c3 = c6^2

    secondcenter = [2/3,-1/3] ∈ firstgmeraspace
    secondhexagonalregion = gethexagonalregion(rot=refrot,crystal=firstgmeracrystal, center=secondcenter, metricspace=firstgmeraspace)
    # secondhexagonalregionfock = quantize(secondhexagonalregion,noofflavourpermode)
    rgshiftedcenter1 = [2/3,-1/3] ∈ firstgmeraspace
    secondrgshiftedhexagonalregion1 = secondhexagonalregion.+rgshiftedcenter1
    rgshiftedcenter2 = [1/3,1/3] ∈ firstgmeraspace
    secondrgshiftedhexagonalregion2 = secondhexagonalregion.+rgshiftedcenter2

    c6recenter = recenter(c6,secondcenter)
    c3recenter = recenter(c3,secondcenter)

    innerhexagonalregion = rankedandgroupoffsets(secondhexagonalregion,1)
    visualize(innerhexagonalregion)

    #siteAregion1
    siteAregion1 = intersect(intersect(secondhexagonalregion,secondrgshiftedhexagonalregion1),secondrgshiftedhexagonalregion2)
    siteAregion1closetriangle =  intersect(siteAregion1,innerhexagonalregion)
    siteAregion1closetrianglefock = quantize(siteAregion1closetriangle,14)

    siteAregion1bdytriangle1 = siteAregion1closetriangle.+(1/2)*rgshiftedcenter1
    siteAregion1bdytriangle2 = siteAregion1closetriangle.+(1/2)*rgshiftedcenter2
    siteAregion1bdytriangle1fock = quantize(siteAregion1bdytriangle1,14)
    siteAregion1bdytriangle2fock = quantize(siteAregion1bdytriangle2,4)

    siteAregion1centertriangle = siteAregion1 - siteAregion1closetriangle - siteAregion1bdytriangle1 - siteAregion1bdytriangle2
    siteAregion1centertrianglefock = quantize(siteAregion1centertriangle,8)

    #siteAregion2
    siteAregion2 = c3recenter*siteAregion1
    siteAregion2closetriangle =  c3recenter*siteAregion1closetriangle
    siteAregion2closetrianglefock = quantize(siteAregion2closetriangle,14)

    siteAregion2bdytriangle1 = c3recenter*siteAregion1bdytriangle1
    siteAregion2bdytriangle2 = c3recenter*siteAregion1bdytriangle2
    siteAregion2bdytriangle1fock = quantize(siteAregion2bdytriangle1,14)
    siteAregion2bdytriangle2fock = quantize(siteAregion2bdytriangle2,4)

    siteAregion2centertriangle = c3recenter*siteAregion1centertriangle
    siteAregion2centertrianglefock = quantize(siteAregion2centertriangle,8)


(siteAregion1closetriangle+siteAregion1bdytriangle1+siteAregion1bdytriangle2+siteAregion1centertriangle+siteAregion2closetriangle+siteAregion2bdytriangle1+siteAregion2bdytriangle2+siteAregion2centertriangle)|>visualize
    #siteAregion3
    siteAregion3 = c3recenter*siteAregion2
    siteAregion3closetriangle =  c3recenter*siteAregion2closetriangle
    siteAregion3closetrianglefock = quantize(siteAregion3closetriangle,14)

    siteAregion3bdytriangle1 = c3recenter*siteAregion2bdytriangle1
    siteAregion3bdytriangle2 = c3recenter*siteAregion2bdytriangle2
    siteAregion3bdytriangle1fock = quantize(siteAregion3bdytriangle1,14)
    siteAregion3bdytriangle2fock = quantize(siteAregion3bdytriangle2,4)

    siteAregion3centertriangle = c3recenter*siteAregion2centertriangle
    siteAregion3centertrianglefock = quantize(siteAregion3centertriangle,8)

    #siteBregion1
    siteBregion1 = c6recenter*siteAregion1
    siteBregion1closetriangle =  c6recenter*siteAregion1closetriangle
    siteBregion1closetrianglefock = quantize(siteBregion1closetriangle,14)

    siteBregion1bdytriangle1 = c6recenter*siteAregion1bdytriangle1
    siteBregion1bdytriangle2 = c6recenter*siteAregion1bdytriangle2
    siteBregion1bdytriangle1fock = quantize(siteBregion1bdytriangle1,4)
    siteBregion1bdytriangle2fock = quantize(siteBregion1bdytriangle2,14)

    siteBregion1centertriangle = c6recenter*siteAregion1centertriangle
    siteBregion1centertrianglefock = quantize(siteBregion1centertriangle,8)

    #siteBregion2
    siteBregion2 = c3recenter*siteBregion1
    siteBregion2closetriangle =  c3recenter*siteBregion1closetriangle
    siteBregion2closetrianglefock = quantize(siteBregion2closetriangle,14)

    siteBregion2bdytriangle1 = c3recenter*siteBregion1bdytriangle1
    siteBregion2bdytriangle2 = c3recenter*siteBregion1bdytriangle2
    siteBregion2bdytriangle1fock = quantize(siteBregion2bdytriangle1,4)
    siteBregion2bdytriangle2fock = quantize(siteBregion2bdytriangle2,14)

    siteBregion2centertriangle = c3recenter*siteBregion1centertriangle
    siteBregion2centertrianglefock = quantize(siteBregion2centertriangle,8)

    #siteBregion3
    siteBregion3 = (c3recenter)*siteBregion2
    siteBregion3closetriangle =  c3recenter*siteBregion2closetriangle
    siteBregion3closetrianglefock = quantize(siteBregion3closetriangle,14)

    siteBregion3bdytriangle1 = c3recenter*siteBregion2bdytriangle1
    siteBregion3bdytriangle2 = c3recenter*siteBregion2bdytriangle2
    siteBregion3bdytriangle1fock = quantize(siteBregion3bdytriangle1,4)
    siteBregion3bdytriangle2fock = quantize(siteBregion3bdytriangle2,14)

    siteBregion3centertriangle = c3recenter*siteBregion2centertriangle
    siteBregion3centertrianglefock = quantize(siteBregion3centertriangle,8)

    secondhexagonalregionfock = (siteAregion1closetrianglefock+siteAregion2closetrianglefock+siteAregion3closetrianglefock
                            + siteBregion1closetrianglefock+siteBregion2closetrianglefock+siteBregion3closetrianglefock
                            + siteAregion1centertrianglefock + siteAregion2centertrianglefock + siteAregion3centertrianglefock
                            + siteBregion1centertrianglefock + siteBregion2centertrianglefock + siteBregion3centertrianglefock
                            + siteAregion1bdytriangle1fock + siteAregion1bdytriangle2fock 
                            + siteAregion2bdytriangle1fock + siteAregion2bdytriangle2fock 
                            + siteAregion3bdytriangle1fock + siteAregion3bdytriangle2fock 
                            + siteBregion1bdytriangle1fock + siteBregion1bdytriangle2fock 
                            + siteBregion2bdytriangle1fock + siteBregion2bdytriangle2fock 
                            + siteBregion3bdytriangle1fock + siteBregion3bdytriangle2fock)|>RegionFock

    localcorrelations = regioncorrelations(firstgmeracorrelations,secondhexagonalregionfock)
    localspectrum = localcorrelations|>eigspec
    localspectrum|>visualize

    localstates = getregionstates(localcorrelations=localcorrelations, grouping=[48,144,48])
    localcourierisometry = localstates[2]|>FockMap
    localfilledisometry = localstates[1]|>FockMap
    localemptyisometry = localstates[3]|>FockMap

    localfilledcorrelations = idmap(localcorrelations|>getoutspace)-(localfilledisometry*localfilledisometry')
    localemptycorrelations = idmap(localcorrelations|>getoutspace)-(localemptyisometry*localemptyisometry')
    localcouriercorrelations = idmap(localcorrelations|>getoutspace)-(localcourierisometry*localcourierisometry')

    localfilledcorrelations[siteBregion1centertrianglefock,siteBregion1centertrianglefock]|>eigspech|>visualize
    # localfilledcorrelations[innerhexagonalregionfock,innerhexagonalregionfock]|>eigspech|>visualize
    localfilledcorrelations[siteAregion1bdytriangle1fock,siteAregion1bdytriangle1fock]|>eigspech|>visualize
    localfilledcorrelations[siteAregion1closetrianglefock,siteAregion1closetrianglefock]|>eigspech|>visualize

    localcouriercorrelations[siteAregion1closetrianglefock,siteAregion1closetrianglefock]|>eigspech|>visualize
    localcouriercorrelations[siteAregion1bdytriangle1fock,siteAregion1bdytriangle1fock]|>eigspech|>visualize
    localcouriercorrelations[siteBregion1centertrianglefock,siteBregion1centertrianglefock]|>eigspech|>visualize