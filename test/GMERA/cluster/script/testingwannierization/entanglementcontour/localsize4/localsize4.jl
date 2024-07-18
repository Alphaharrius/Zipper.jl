using Zipper
using LinearAlgebra

setmaxthreads(8)

systemsize=32
correlations,H = generatesystem(-0.1+0im,0.1+ 0im,-1 + 0im,0im,systemsize)
crystalfock = correlations|>getoutspace

scale = Scale([4 0; 0 4], crystalfock|>getcrystal|>getspace)
@info("Performing rgblocking...",scale)
@info("Generating rgblocking transformation...")
blocker = @time scale * crystalfock
@info("Performing rgblocking on correlations...")
blockedcorrelations = @time blocker * correlations * blocker'

function localwannierizationresult(blockedcorrelations)

    blockedcrystal = blockedcorrelations|>getoutspace|>getcrystal
    blockedcrystalfock = blockedcorrelations|>getoutspace
    blockedspace::RealSpace = blockedcrystal|>getspace

    refrot = inv([2/3 -1/3; -1/3 2/3]')
    c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
    c3 = c6^2

    firstcenter = [0,0] ∈ blockedspace
    firsthexagonalregion = gethexagonalregion(rot = refrot,crystal=blockedcrystal, center=firstcenter, metricspace=blockedspace)
    firsthexagonalregionfock = quantize(firsthexagonalregion,1)

    localcorrelations = regioncorrelations(blockedcorrelations,firsthexagonalregionfock)
    localspectrum = localcorrelations|>eigspec

    evec = localspectrum|>geteigenvectors
    modeevalspair = localspectrum|>geteigenvalues

    entanglementcontour = []
    for rsmode in firsthexagonalregionfock|>orderedmodes
        save = []
        for (mode,eval) in modeevalspair
            if ((norm(eval)>0.00000001) & (norm(eval)<1-0.00000001))
                append!(save,norm(evec[rsmode,mode])^2*(-norm(eval)*log(norm(eval))-(1-norm(eval))*log(1-norm(eval))))
            end
        end
        rsoffset = (rsmode|>getattr(:b))+(rsmode|>getattr(:r))
        append!(entanglementcontour,tuple([rsoffset,sum(save)]))
    end

    sortedentangelementcontour = sort([(rsoffset,round(val,digits=4)) for (rsoffset,val) in entanglementcontour],by=x->x[2])
    ref = sortedentangelementcontour[1][2]
    samecontrib = []
    groupedoffsetbyentanglementcontour = []
    for pair in sortedentangelementcontour
        if pair[2] == ref
            push!(samecontrib,pair[1])
        else
            ref = pair[2]
            push!(groupedoffsetbyentanglementcontour,samecontrib)
            samecontrib = []
            push!(samecontrib,pair[1])
        end
    end
    for chosen in range(1,length(groupedoffsetbyentanglementcontour))
        frozenseedsoffset = Subset([offset for offsets in groupedoffsetbyentanglementcontour[1:chosen] for offset in offsets])
        courierseedsoffset = firsthexagonalregion-frozenseedsoffset
        courierseedsfock = quantize(courierseedsoffset,1)

        nooffilledmodes = div(length(frozenseedsoffset),2)
        noofemptymodes = nooffilledmodes
        noofcouriermodes = length(courierseedsoffset)
        println(nooffilledmodes+noofemptymodes+noofcouriermodes)

        localstates = getregionstates(localcorrelations=localcorrelations, grouping=[nooffilledmodes,noofcouriermodes,noofemptymodes])
        localcourierisometry = localstates[2]|>FockMap


        courierseeds = idmap(firsthexagonalregionfock, firsthexagonalregionfock)[:,FockSpace(mode for mode in courierseedsfock)]
        localwanniercourier,minsvd,traceofsv=localwannierization(localcourierisometry, courierseeds)

        fiodir("/Users/slwongag/Desktop/data/testingwannierization/entanglementcontour/systemsize$systemsize/localsize4/chosen_$chosen")
        fiosave(minsvd, name="minsvd")
        fiosave(traceofsv, name="traceofsv")
        fiosave(frozenseedsoffset, name="frozenseedsoffset")
        fiosave(courierseedsoffset, name="courierseedsoffset")
    end
end

localwannierizationresult(blockedcorrelations)

