using Zipper
using LinearAlgebra
using Plots

# function save(
#     states::RegionState{2};
#     visualspace::RealSpace = states|>getoutspace|>getregion|>getspace|>euclidean,
#     markersize::Integer = 1, logscale::Real = 1)

#     function savestate(order::Integer, statemap::SparseFockMap)
#         r = [visualspace*getpos(mode)|>vec for mode in statemap|>getoutspace]
#         x = [v|>first for v in r]
#         y = [v|>last for v in r]
#         f = ((statemap[mode, :]|>rep)[1, 1] for mode in statemap|>getoutspace)
#         color = [convert(RGB{Float64}, HSV(angle(v) / 2π * 360, 1, 1)) for v in f]
#         size = [max(markersize*(abs(v)^logscale), 1e-5) for v in f]
#         hover = [trunc(v|>real, digits=3)+trunc(v|>imag, digits=3)*im|>string for v in f]
#         scatter(
#             x, y, markersize=size, markercolor=color, hover=hover,
#             markerstrokealpha=0.01,
#             aspect_ratio=1.1, xlims=(minimum(x)-1, maximum(x)+1), ylims=(minimum(y)-1, maximum(y)+1),
#             title="mode $order")
#     end

#     source = [(n, mode, state) for (n, (mode, state)) in states|>enumerate]
#     for scatter in [savestate(n, state) for (n, _, state) in source]
#         display(plot(scatter))
#     end
#     # p = plot((plotstate(n, state) for (n, _, state) in source)..., legend=false)
#     # @info "Showing corresponding modes..."
#     # showmodes(mode for (_, mode, _) in source)
#     # return p
# end

systemsize=32
correlations,H = generatesystem(-0.3 + 0im,0.3 + 0im,-1 + 0im,0im,systemsize)
crystalfock = correlations|>getoutspace

scale = Scale([8 0; 0 8], crystalfock|>getcrystal|>getspace)
@info("Performing rgblocking...",scale)
@info("Generating rgblocking transformation...")
blocker = @time scale * crystalfock
@info("Performing rgblocking on correlations...")
blockedcorrelations = @time blocker * correlations * blocker'
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

siteAregion1 = intersect(intersect(firsthexagonalregion,firstrgshiftedhexagonalregion1 ),firstrgshiftedhexagonalregion2)
siteAregionfock1 = quantize(siteAregion1,1)
siteAregion2 = c3*siteAregion1
siteAregionfock2 = quantize(siteAregion2,1)
siteAregion3 = (c3)*siteAregion2
siteAregionfock3 = quantize(siteAregion3,1)

siteBregion1 = (c6)*siteAregion1
siteBregionfock1 = quantize(siteBregion1,1)
siteBregion2 = (c3)*siteBregion1
siteBregionfock2 = quantize(siteBregion2,1)
siteBregion3 = (c3)*siteBregion2
siteBregionfock3 = quantize(siteBregion3,1)

# allregion = siteAregion1+siteAregion2+siteAregion3+siteBregion1+siteBregion2+siteBregion3
# allregion|>visualize

localcorrelations = regioncorrelations(blockedcorrelations,firsthexagonalregionfock)
localcorrelations|>eigspec|>visualize
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[69, 246, 69])
# save(localstates[2])
localcourierisometry = localstates[2]|>FockMap
courierproj = localcourierisometry*localcourierisometry'
columns(rows(courierproj,siteAregionfock1),siteAregionfock1)|>eigspech|>visualize

localAcorrelations1 = columns(rows(courierproj,siteAregionfock1),siteAregionfock1)
localAcorrelations1|>eigspec|>visualize
localAstates1 = getregionstates(localcorrelations=localAcorrelations1, grouping=[23, 41])
localAseeds1 = localAstates1[2]

localAcorrelations2 = columns(rows(courierproj,siteAregionfock2),siteAregionfock2)
localAcorrelations2|>eigspec|>visualize
localAstates2 = getregionstates(localcorrelations=localAcorrelations2, grouping=[23, 41])
localAseeds2 = localAstates2[2]


localAcorrelations3 = columns(rows(courierproj,siteAregionfock3),siteAregionfock3)
localAcorrelations3|>eigspec|>visualize
localAstates3 = getregionstates(localcorrelations=localAcorrelations3, grouping=[23, 41])
localAseeds3 = localAstates3[2]

allAseedsstate = localAseeds1+localAseeds2+localAseeds3

localBcorrelations1 = columns(rows(courierproj,siteBregionfock1),siteBregionfock1)
localBcorrelations1|>eigspec|>visualize
localBstates1 = getregionstates(localcorrelations=localBcorrelations1, grouping=[23, 41])
localBseeds1 = localBstates1[2]

localBcorrelations2 = columns(rows(courierproj,siteBregionfock2),siteBregionfock2)
localBcorrelations2|>eigspec|>visualize
localBstates2 = getregionstates(localcorrelations=localBcorrelations2, grouping=[23, 41])
localBseeds2 = localBstates2[2]

localBcorrelations3 = columns(rows(courierproj,siteBregionfock3),siteBregionfock3)
localBcorrelations3|>eigspec|>visualize
localBstates3 = getregionstates(localcorrelations=localBcorrelations3, grouping=[23, 41])
localBseeds3 = localBstates3[2]

allBseedsstate = localBseeds1+localBseeds2+localBseeds3

allseedsstates = (allAseedsstate+allBseedsstate)
courierseeds = allseedsstates|>FockMap
localcourierisometry

wanniercourier = localwannierization(localcourierisometry, courierseeds)
(wanniercourier'*localcorrelations*wanniercourier)|>eigspec|>visualize

localwannierresults =  Dict(:wanniercourier => wanniercourier, :courierseeds => courierseeds)

wannierinfos =  Dict(firsthexagonalregionfock=>localwannierresults)

shiftedfirstcenter1 = [1,0] ∈ blockedspace
shiftedfirstcenter2 = [0,1] ∈ blockedspace
shiftedfirstcenter3 = [1,1] ∈ blockedspace

firstshiftedhexagonalcenterlist = [shiftedfirstcenter1, shiftedfirstcenter2, shiftedfirstcenter3]

for hexagonalcenter in firstshiftedhexagonalcenterlist
    shiftedfirsthexagonalregion = firsthexagonalregion.+hexagonalcenter
    shiftedfirsthexagonalregionfock = quantize(shiftedfirsthexagonalregion,1)
    c6recenter = recenter(c6,hexagonalcenter)
    c3recenter = recenter(c3,hexagonalcenter)

    shiftedsiteAregion1 = siteAregion1.+hexagonalcenter
    shiftedsiteAregionfock1 = quantize(shiftedsiteAregion1,1)
    shiftedsiteAregion2 = (c3recenter)*shiftedsiteAregion1
    shiftedsiteAregionfock2 = quantize(shiftedsiteAregion2,1)
    shiftedsiteAregion3 = (c3recenter)*shiftedsiteAregion2
    shiftedsiteAregionfock3 = quantize(shiftedsiteAregion3,1)

    shiftedsiteBregion1 = (c6recenter)*shiftedsiteAregion1
    shiftedsiteBregionfock1 = quantize(shiftedsiteBregion1,1)
    shiftedsiteBregion2 = (c3recenter)*shiftedsiteBregion1
    shiftedsiteBregionfock2 = quantize(shiftedsiteBregion2,1)
    shiftedsiteBregion3 = (c3recenter)*shiftedsiteBregion2
    shiftedsiteBregionfock3 = quantize(shiftedsiteBregion3,1)

    shiftedlocalcorrelations = regioncorrelations(blockedcorrelations,shiftedfirsthexagonalregionfock)
    shiftedlocalstates = getregionstates(localcorrelations=shiftedlocalcorrelations, grouping=[69, 246, 69])
    shiftedlocalcourierisometry = shiftedlocalstates[2]|>FockMap
    shiftedcourierproj = shiftedlocalcourierisometry*shiftedlocalcourierisometry'

    shiftedlocalAcorrelations1 = columns(rows(shiftedcourierproj,shiftedsiteAregionfock1),shiftedsiteAregionfock1)
    shiftedlocalAstates1 = getregionstates(localcorrelations=shiftedlocalAcorrelations1, grouping=[23, 41])
    shiftedlocalAseeds1 = shiftedlocalAstates1[2]

    shiftedlocalAcorrelations2 = columns(rows(shiftedcourierproj,shiftedsiteAregionfock2),shiftedsiteAregionfock2)
    shiftedlocalAstates2 = getregionstates(localcorrelations=shiftedlocalAcorrelations2, grouping=[23, 41])
    shiftedlocalAseeds2 = shiftedlocalAstates2[2]

    shiftedlocalAcorrelations3 = columns(rows(shiftedcourierproj,shiftedsiteAregionfock3),shiftedsiteAregionfock3)
    shiftedlocalAstates3 = getregionstates(localcorrelations=shiftedlocalAcorrelations3, grouping=[23, 41])
    shiftedlocalAseeds3 = shiftedlocalAstates3[2]

    shiftedallAseedsstate = shiftedlocalAseeds1+shiftedlocalAseeds2+shiftedlocalAseeds3

    shiftedlocalBcorrelations1 = columns(rows(shiftedcourierproj,shiftedsiteBregionfock1),shiftedsiteBregionfock1)
    shiftedlocalBstates1 = getregionstates(localcorrelations=shiftedlocalBcorrelations1, grouping=[23, 41])
    shiftedlocalBseeds1 = shiftedlocalBstates1[2]

    shiftedlocalBcorrelations2 = columns(rows(shiftedcourierproj,shiftedsiteBregionfock2),shiftedsiteBregionfock2)
    shiftedlocalBstates2 = getregionstates(localcorrelations=shiftedlocalBcorrelations2, grouping=[23, 41])
    shiftedlocalBseeds2 = shiftedlocalBstates2[2]

    shiftedlocalBcorrelations3 = columns(rows(shiftedcourierproj,shiftedsiteBregionfock3),shiftedsiteBregionfock3)
    shiftedlocalBstates3 = getregionstates(localcorrelations=shiftedlocalBcorrelations3, grouping=[23, 41])
    shiftedlocalBseeds3 = shiftedlocalBstates3[2]

    shiftedallBseedsstate = shiftedlocalBseeds1+shiftedlocalBseeds2+shiftedlocalBseeds3

    shiftedallseedsstates = (shiftedallAseedsstate+shiftedallBseedsstate)
    shiftedcourierseeds = shiftedallseedsstates|>FockMap

    shiftedwanniercourier = localwannierization(shiftedlocalcourierisometry, shiftedcourierseeds)

    shiftedlocalwannierresults =  Dict(:wanniercourier => shiftedwanniercourier, :courierseeds => shiftedcourierseeds)
    wannierinfos[shiftedfirsthexagonalregionfock] =  shiftedlocalwannierresults
end

ref = [quantize(firsthexagonalregion.+hexagonalcenter,1) for hexagonalcenter in firstshiftedhexagonalcenterlist]
firsthexagonalregionfocklist = [firsthexagonalregionfock,ref...]

extendedwannierizedcourier =  sum(wannierinfos[regionfock][:wanniercourier] for regionfock in firsthexagonalregionfocklist)
                
origin = [0, 0] ∈ blockedspace
refunictcellfockcourier = FockSpace(Subset(mode for mode in extendedwannierizedcourier |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
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
secondhexagonalregionfock = quantize(secondhexagonalregion,41)
rgshiftedcenter1 = [2/3,-1/3] ∈ firstrgedspace
secondrgshiftedhexagonalregion1 = secondhexagonalregion.+rgshiftedcenter1
rgshiftedcenter2 = [1/3,1/3] ∈ firstrgedspace
secondrgshiftedhexagonalregion2 = secondhexagonalregion.+rgshiftedcenter2

c6recenter = recenter(c6,secondcenter)
c3recenter = recenter(c3,secondcenter)

siteAregion1 = intersect(intersect(secondhexagonalregion,secondrgshiftedhexagonalregion1 ),secondrgshiftedhexagonalregion2)
siteAregionfock1 = quantize(siteAregion1,41)
siteAregion2 = (c3recenter)*siteAregion1
siteAregionfock2 = quantize(siteAregion2,41)
siteAregion3 = (c3recenter)*siteAregion2
siteAregionfock3 = quantize(siteAregion3,41)

siteBregion1 = (c6recenter)*siteAregion1
siteBregionfock1 = quantize(siteBregion1,41)
siteBregion2 = (c3recenter)*siteBregion1
siteBregionfock2 = quantize(siteBregion2,41)
siteBregion3 = (c3recenter)*siteBregion2
siteBregionfock3 = quantize(siteBregion3,41)


localcorrelations = regioncorrelations(firstrgedcorrelations,secondhexagonalregionfock)
localspectrum = localcorrelations|>eigspech
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[63, 120, 63])
# save(localstates[2])
localcourierisometry = localstates[2]|>FockMap
courierproj = localcourierisometry*localcourierisometry'
columns(rows(courierproj,siteAregionfock1),siteAregionfock1)|>eigspech|>visualize

localAcorrelations1 = columns(rows(courierproj,siteAregionfock1),siteAregionfock1)
localAcorrelations1|>eigspec|>visualize
localAstates1 = getregionstates(localcorrelations=localAcorrelations1, grouping=[21, 20])
localAseeds1 = localAstates1[2]

localAcorrelations2 = columns(rows(courierproj,siteAregionfock2),siteAregionfock2)
localAcorrelations2|>eigspec|>visualize
localAstates2 = getregionstates(localcorrelations=localAcorrelations2, grouping=[21, 20])
localAseeds2 = localAstates2[2]


localAcorrelations3 = columns(rows(courierproj,siteAregionfock3),siteAregionfock3)
localAcorrelations3|>eigspec|>visualize
localAstates3 = getregionstates(localcorrelations=localAcorrelations3, grouping=[21, 20])
localAseeds3 = localAstates3[2]

allAseedsstate = localAseeds1+localAseeds2+localAseeds3

localBcorrelations1 = columns(rows(courierproj,siteBregionfock1),siteBregionfock1)
localBcorrelations1|>eigspec|>visualize
localBstates1 = getregionstates(localcorrelations=localBcorrelations1, grouping=[21, 20])
localBseeds1 = localBstates1[2]

localBcorrelations2 = columns(rows(courierproj,siteBregionfock2),siteBregionfock2)
localBcorrelations2|>eigspec|>visualize
localBstates2 = getregionstates(localcorrelations=localBcorrelations2, grouping=[21, 20])
localBseeds2 = localBstates2[2]

localBcorrelations3 = columns(rows(courierproj,siteBregionfock3),siteBregionfock3)
localBcorrelations3|>eigspec|>visualize
localBstates3 = getregionstates(localcorrelations=localBcorrelations3, grouping=[21, 20])
localBseeds3 = localBstates3[2]

allBseedsstate = localBseeds1+localBseeds2+localBseeds3

allseedsstates = (allAseedsstate+allBseedsstate)
courierseeds = allseedsstates|>FockMap
localcourierisometry

wanniercourier = localwannierization(localcourierisometry, courierseeds)
(wanniercourier'*localcorrelations*wanniercourier)|>eigspec|>visualize

localwannierresults =  Dict(:wanniercourier => wanniercourier, :courierseeds => courierseeds)

wannierinfos =  Dict(secondhexagonalregionfock=>localwannierresults)

shiftedsecondcenter1 = [0,1] ∈ firstrgedspace
shiftedsecondcenter2 = [-1,1] ∈ firstrgedspace

secondshiftedhexagonalcenterlist = [shiftedsecondcenter1, shiftedsecondcenter2]

for hexagonalcenter in secondshiftedhexagonalcenterlist
    shiftedsecondhexagonalregion = secondhexagonalregion.+hexagonalcenter
    shiftedsecondhexagonalregionfock = quantize(shiftedsecondhexagonalregion,41)
    c6recenter = recenter(c6,secondcenter+hexagonalcenter)
    c3recenter = recenter(c3,secondcenter+hexagonalcenter)

    shiftedsiteAregion1 = siteAregion1.+hexagonalcenter
    shiftedsiteAregionfock1 = quantize(shiftedsiteAregion1,41)
    shiftedsiteAregion2 = (c3recenter)*shiftedsiteAregion1
    shiftedsiteAregionfock2 = quantize(shiftedsiteAregion2,41)
    shiftedsiteAregion3 = (c3recenter)*shiftedsiteAregion2
    shiftedsiteAregionfock3 = quantize(shiftedsiteAregion3,41)

    shiftedsiteBregion1 = (c6recenter)*shiftedsiteAregion1
    shiftedsiteBregionfock1 = quantize(shiftedsiteBregion1,41)
    shiftedsiteBregion2 = (c3recenter)*shiftedsiteBregion1
    shiftedsiteBregionfock2 = quantize(shiftedsiteBregion2,41)
    shiftedsiteBregion3 = (c3recenter)*shiftedsiteBregion2
    shiftedsiteBregionfock3 = quantize(shiftedsiteBregion3,41)

    shiftedlocalcorrelations = regioncorrelations(firstrgedcorrelations,shiftedsecondhexagonalregionfock)
    shiftedlocalstates = getregionstates(localcorrelations=shiftedlocalcorrelations, grouping=[63, 120, 63])
    shiftedlocalcourierisometry = shiftedlocalstates[2]|>FockMap
    shiftedcourierproj = shiftedlocalcourierisometry*shiftedlocalcourierisometry'

    shiftedlocalAcorrelations1 = columns(rows(shiftedcourierproj,shiftedsiteAregionfock1),shiftedsiteAregionfock1)
    shiftedlocalAstates1 = getregionstates(localcorrelations=shiftedlocalAcorrelations1, grouping=[21, 20])
    shiftedlocalAseeds1 = shiftedlocalAstates1[2]

    shiftedlocalAcorrelations2 = columns(rows(shiftedcourierproj,shiftedsiteAregionfock2),shiftedsiteAregionfock2)
    shiftedlocalAstates2 = getregionstates(localcorrelations=shiftedlocalAcorrelations2, grouping=[21, 20])
    shiftedlocalAseeds2 = shiftedlocalAstates2[2]

    shiftedlocalAcorrelations3 = columns(rows(shiftedcourierproj,shiftedsiteAregionfock3),shiftedsiteAregionfock3)
    shiftedlocalAstates3 = getregionstates(localcorrelations=shiftedlocalAcorrelations3, grouping=[21, 20])
    shiftedlocalAseeds3 = shiftedlocalAstates3[2]

    shiftedallAseedsstate = shiftedlocalAseeds1+shiftedlocalAseeds2+shiftedlocalAseeds3

    shiftedlocalBcorrelations1 = columns(rows(shiftedcourierproj,shiftedsiteBregionfock1),shiftedsiteBregionfock1)
    shiftedlocalBstates1 = getregionstates(localcorrelations=shiftedlocalBcorrelations1, grouping=[21, 20])
    shiftedlocalBseeds1 = shiftedlocalBstates1[2]

    shiftedlocalBcorrelations2 = columns(rows(shiftedcourierproj,shiftedsiteBregionfock2),shiftedsiteBregionfock2)
    shiftedlocalBstates2 = getregionstates(localcorrelations=shiftedlocalBcorrelations2, grouping=[21, 20])
    shiftedlocalBseeds2 = shiftedlocalBstates2[2]

    shiftedlocalBcorrelations3 = columns(rows(shiftedcourierproj,shiftedsiteBregionfock3),shiftedsiteBregionfock3)
    shiftedlocalBstates3 = getregionstates(localcorrelations=shiftedlocalBcorrelations3, grouping=[21, 20])
    shiftedlocalBseeds3 = shiftedlocalBstates3[2]

    shiftedallBseedsstate = shiftedlocalBseeds1+shiftedlocalBseeds2+shiftedlocalBseeds3

    shiftedallseedsstates = (shiftedallAseedsstate+shiftedallBseedsstate)
    shiftedcourierseeds = shiftedallseedsstates|>FockMap

    shiftedwanniercourier = localwannierization(shiftedlocalcourierisometry, shiftedcourierseeds)

    shiftedlocalwannierresults =  Dict(:wanniercourier => shiftedwanniercourier, :courierseeds => shiftedcourierseeds)
    wannierinfos[shiftedsecondhexagonalregionfock] =  shiftedlocalwannierresults
end

ref = [quantize(secondhexagonalregion.+hexagonalcenter,41) for hexagonalcenter in secondshiftedhexagonalcenterlist]
secondhexagonalregionfocklist = [secondhexagonalregionfock,ref...]

extendedwannierizedcourier =  sum(wannierinfos[regionfock][:wanniercourier] for regionfock in secondhexagonalregionfocklist)
                                
origin = [0, 0] ∈ firstrgedspace
refunictcellfockcourier = FockSpace(Subset(mode for mode in extendedwannierizedcourier |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
                    
wanniercourierisometry = globalwannierfunction(firstrgedcorrelations,extendedwannierizedcourier[:,refunictcellfockcourier])
                    
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
thirdcenter = [0,0] ∈ secondrgedspace
thirdhexagonalregion = gethexagonalregion(rot = refrot,crystal=secondrgedcrystal, center=thirdcenter, metricspace=secondrgedspace)
thirdhexagonalregionfock = quantize(thirdhexagonalregion,20)
localcorrelations = regioncorrelations(secondrgedcorrelations,thirdhexagonalregionfock)
localspectrum = localcorrelations|>eigspech

display(localspectrum|>visualize)


# Subset([(m|>getattr(:b))+(m|>getattr(:r)) for m in wanniercourier|>getinspace|>orderedmodes])|>visualize
# [(m|>getattr(:flavor)) for m in wanniercourier|>getinspace|>orderedmodes]

# A1 = [(norm(columns(rows(wanniercourier,siteAregionfock1),FockSpace(m)))^2,m) for m in rows(wanniercourier,siteAregionfock1)|>getinspace|>orderedmodes]
# A2 = [(norm(columns(rows(wanniercourier,siteAregionfock2),FockSpace(m)))^2,m) for m in rows(wanniercourier,siteAregionfock2)|>getinspace|>orderedmodes]
# A3 = [(norm(columns(rows(wanniercourier,siteAregionfock3),FockSpace(m)))^2,m) for m in rows(wanniercourier,siteAregionfock3)|>getinspace|>orderedmodes]
# B1 = [(norm(columns(rows(wanniercourier,siteBregionfock1),FockSpace(m)))^2,m) for m in rows(wanniercourier,siteBregionfock1)|>getinspace|>orderedmodes]
# B2 = [(norm(columns(rows(wanniercourier,siteBregionfock2),FockSpace(m)))^2,m) for m in rows(wanniercourier,siteBregionfock2)|>getinspace|>orderedmodes]
# B3 = [(norm(columns(rows(wanniercourier,siteBregionfock3),FockSpace(m)))^2,m) for m in rows(wanniercourier,siteBregionfock3)|>getinspace|>orderedmodes]

# A1ranked = sort(A1,rev=true,by=first)
# A2ranked = sort(A2,rev=true,by=first)
# A3ranked = sort(A3,rev=true,by=first)
# B1ranked = sort(B1,rev=true,by=first)
# B2ranked = sort(B2,rev=true,by=first)
# B3ranked = sort(B3,rev=true,by=first)


# A1rankedmodes = Subset([pair[2] for pair in A1ranked[1:38]])
# A2rankedmodes = Subset([pair[2] for pair in A2ranked[1:38]])
# A3rankedmodes = Subset([pair[2] for pair in A3ranked[1:38]])
# B1rankedmodes = Subset([pair[2] for pair in B1ranked[1:38]])
# B2rankedmodes = Subset([pair[2] for pair in B2ranked[1:38]])
# B3rankedmodes = Subset([pair[2] for pair in B3ranked[1:38]])

# A1ranked[1:39]
# B3ranked[1:38]

# A1rankedmodes+A2rankedmodes+A3rankedmodes+B1rankedmodes+B2rankedmodes+B3rankedmodes


# (wanniercourier'*localcorrelations*wanniercourier)|>eigspec|>visualize

# save(RegionState(wanniercourier))


@info("Computing local correlations...")

