using Plotly, SmithNormalForm, LinearAlgebra, OrderedCollections, SparseArrays, Combinatorics
using Zipper

setmaxthreads(Threads.nthreads())


triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')
kspace = convert(MomentumSpace, triangular)

pa = [1/3, 2/3] ∈ triangular
pb = [2/3, 1/3] ∈ triangular
pc = (pa + pb) / 2
spatialsnappingcalibration((pa, pb, pc))

c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])

unitcell = Subset(pa, pb)
crystal = Crystal(unitcell, [48, 48])
reciprocalhashcalibration(crystal.sizes)

modes::Subset{Mode} = quantize(unitcell, 1)|>orderedmodes
m0, m1 = members(modes)

tₙ = ComplexF64(-1)
tₕ = 0.1im

nearestneighbor = [
    (m0, m1) => tₙ,
    (m0, setattr(m1, :r => Point([-1, 0], triangular))) => tₙ,
    (m0, setattr(m1, :r => Point([0, 1], triangular))) => tₙ]

haldane = [
    (m0, setattr(m0, :r => Point([1, 1], triangular))) => tₕ,
    (m0, setattr(m0, :r => Point([-1, 0], triangular))) => tₕ,
    (m0, setattr(m0, :r => Point([0, -1], triangular))) => tₕ,
    (m1, setattr(m1, :r => Point([1, 1], triangular))) => -tₕ,
    (m1, setattr(m1, :r => Point([-1, 0], triangular))) => -tₕ,
    (m1, setattr(m1, :r => Point([0, -1], triangular))) => -tₕ]

bonds::FockMap = bondmap([nearestneighbor..., haldane...])

energyspectrum = @time computeenergyspectrum(bonds, crystal=crystal)
energyspectrum|>visualize

groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=1)
groundstates|>visualize

groundstateprojector = groundstates|>crystalprojector
correlations = idmap(groundstateprojector|>getoutspace) - groundstateprojector

H = CrystalFockMap(energyspectrum)

@info("Starting RG...")
crystalfock = correlations|>getoutspace

scale = Scale([8 0; 0 8], crystalfock|>getcrystal|>getspace)
@info("Performing unitcellblocking...")
@info("Generating unitcellblocking transformation...")
blocker48 = @time scale * crystalfock
@info("Performing unitcellblocking on correlations...")
blockedcorrelations48 = @time blocker48 * correlations * blocker48'
blockedcrystalfock48 = blockedcorrelations48|>getoutspace
blockedcrystal48::Crystal = blockedcrystalfock48|>getcrystal
blockedspace48::RealSpace = blockedcrystal48|>getspace

# localfock = RegionFock(blockedcrystalfock48|>unitcellfock)
# visualize(regioncorrelations(blockedcorrelations48,localfock)|>eigspech)
# localiso = localisometries(blockedcorrelations48,localfock,selectionstrategy=modeselectionbycount(16))

# refcorr = regioncorrelations(blockedcorrelations48,RegionFock(blockedcrystalfock48|>unitcellfock))
# Xpos = genposX(refcorr)
# Ypos = genposY(refcorr)
# localfock = RegionFock(blockedcrystalfock48|>unitcellfock)
# visualize(regioncorrelations(blockedcorrelations48,localfock)|>eigspech)
# localiso = localisometries(blockedcorrelations48,localfock,selectionstrategy=modeselectionbycount(5))
# projfilled|>getinspace|>orderedmodes
# projfilled = localiso[:filled]*localiso[:filled]'
# projcourier = localiso[:courier]*localiso[:courier]'
# localiso[:courier][:isometry]|>getinspace|>orderedmodes|>length
# chosencouriermodes = modewifdiagpairs(projector=projcourier,noofchosen=localiso[:courier]|>getinspace|>orderedmodes|>length)
# couriermodes = Subset([pair[1] for pair in chosencouriermodes])
# courierseeds = idmap(localfock, localfock)[:,FockSpace(mode for mode in couriermodes)]
# wanniercourier = localwannierization(localiso[:courier], courierseeds)

# frozenmodes = Subset(projfilled|>getinspace)-couriermodes
# chosenfrozenmodes = modewifdiagpairs(projector=projfilled,modes=frozenmodes)
# groupedfrozenmodes = groupmodesbydiag(chosenfrozenmodes)
# filledmodes = Subset([pair[2][i][2] for pair in groupedfrozenmodes for i in range(Int((pair[2]|>length)/2)+1,Int((pair[2]|>length)))])
# filledmodes = frozenmodes[3:7]

# filledseeds = idmap(localfock, localfock)[:,FockSpace(mode for mode in filledmodes)]
# wannierfilled = localwannierization(localiso[:filled][:isometry], filledseeds)

# Subset(mode for mode in projfilled|>getinspace)
# sort([norm(projfilled[mode,mode]) for mode in projfilled|>getinspace])[end-12:end]
# sort([norm(d) for d in diag(projfilled|>rep)])
# sort([norm(d) for d in diag(projcourier|>rep)])[end-22:end]
# sort([calnorm(projcourier[i,:]|>rep) for i in range(1,72)])
# sort([norm(projcourier[i,:]) for i in range(1,32)])
# sort([norm(projcourier[:,i]) for i in range(1,72)])
# sort([norm(projfilled[:,i]) for i in range(1,32)])[end-10:end]
# visualize(localiso[:filled][:isometry]'*Ypos*localiso[:filled][:isometry]|>eigspech)
# # visualize(localiso[:filled][:isometry]*localiso[:filled][:isometry]')
# uni = localiso[:filled][:isometry]'*Xpos*localiso[:filled][:isometry]|>eigspech|>geteigenvectors
# rotlocaliso = localiso[:filled][:isometry]*uni[:,5:7]
# visualize(rotlocaliso'*Ypos*rotlocaliso|>eigspech)

# [ComplexF64(0),ComplexF64(0),ComplexF64(0)]
# spdiagm([0,1,2])
# sparse(I,4,4)
# [m|>getattr(:b) for m in blockedcrystalfock48|>unitcellfock|>orderedmodes]
# RegionFock(blockedcrystalfock48|>unitcellfock)
# Complex(0)+1
# id = FockMap(RegionFock(blockedcrystalfock48|>unitcellfock),RegionFock(blockedcrystalfock48|>unitcellfock),2*sparse(I,8,8))
# (id*refcorr)|>rep
# [m|>getattr(:b)  for m in FockMap(RegionFock(blockedcrystalfock48|>unitcellfock),RegionFock(blockedcrystalfock48|>unitcellfock),sparse(I,8,8))|>getinspace|>orderedmodes]
# Point([1, 1], triangular).pos[1]

function gmera(correlations,reftransistionmap)
    @info("Starting RG...")
    crystalfock = correlations|>getoutspace

    rgscale = Scale([2 0; 0 2], crystalfock|>getcrystal|>getspace)
    @info("Performing rgblocking...")
    @info("Generating rgblocking transformation...")
    rgblocker = @time rgscale * crystalfock
    @info("Performing rgblocking on correlations...")
    rgblockedcorrelations = @time rgblocker * correlations * rgblocker'
    rgblockedcrystalfock = rgblockedcorrelations|>getoutspace
    rgblockedcrystal::Crystal = rgblockedcrystalfock|>getcrystal
    rgblockedspace::RealSpace = rgblockedcrystal|>getspace

    transistionmap = reftransistionmap*rgblocker'

    firstcenterlist = [[1/2,0] ∈ rgblockedspace,[-1/2,0] ∈ rgblockedspace]
    secondcenterlist = [[0,1/2] ∈ rgblockedspace,[0,-1/2] ∈ rgblockedspace]
    thirdcenterlist = [[1/2,1/2] ∈ rgblockedspace,[-1/2,-1/2] ∈ rgblockedspace, [-1/2,1/2] ∈ rgblockedspace,[1/2,-1/2] ∈ rgblockedspace]
    finalcenterlist = [[0,0] ∈ rgblockedspace]
    @info ("1st gmera step...")
    gmera1 = @time gmerastep(rgblockedcorrelations,rgblockedcorrelations,firstcenterlist, modeselectionbycount(3))
    # gmera1 = @time gmerastep1(rgblockedcorrelations,firstcenterlist)
    @info ("1st gmera approximation to correlations...")
    gmera1approx = transistionmap*gmera1[:emptyisometry]*gmera1[:emptyisometry]'*transistionmap'
    transistionmap = transistionmap*gmera1[:courierisometry]

    @info ("2nd gmera step...")
    gmera2 = @time gmerastep(rgblockedcorrelations,gmera1[:correlations],secondcenterlist, modeselectionbycount(3))
    @info ("2nd gmera approximation to correlations...")
    gmera2approx = transistionmap*gmera2[:emptyisometry]*gmera2[:emptyisometry]'*transistionmap'
    transistionmap = transistionmap*gmera2[:courierisometry]

    @info ("3rd gmera step...")
    gmera3 = @time gmerastep(rgblockedcorrelations,gmera2[:correlations],thirdcenterlist, modeselectionbycount(3))
    @info ("3rd gmera approximation to correlations...")
    gmera3approx = transistionmap*gmera3[:emptyisometry]*gmera3[:emptyisometry]'*transistionmap'
    transistionmap = transistionmap*gmera3[:courierisometry]

    @info ("final gmera step...")
    gmera4 = @time gmerastep(rgblockedcorrelations,gmera3[:correlations],finalcenterlist, modeselectionbycount(3))
    @info ("4th gmera approximation to correlations...")
    gmera4approx = transistionmap*gmera4[:emptyisometry]*gmera4[:emptyisometry]'*transistionmap'
    transistionmap = transistionmap*gmera4[:courierisometry]
    
    # return gmera1
    return Dict(
        :rgblockedmap => rgblocker,
        :gmera1stemptyisometry => gmera1[:emptyisometry],
        :gmera1stfilledisometry => gmera1[:filledisometry],
        :gmera1stcourierisometry => gmera1[:courierisometry],
        :gmera1stcorrelations => gmera1[:correlations],
        :gmera1stapprox => gmera1approx,
        :gmera2ndemptyisometry => gmera2[:emptyisometry],
        :gmera2ndfilledisometry => gmera2[:filledisometry],
        :gmera2ndcourierisometry => gmera2[:courierisometry],
        :gmera2ndcorrelations => gmera2[:correlations],
        :gmera2ndapprox => gmera2approx,
        :gmera3rdemptyisometry => gmera3[:emptyisometry],
        :gmera3rdfilledisometry => gmera3[:filledisometry],
        :gmera3rdcourierisometry => gmera3[:courierisometry],
        :gmera3rdcorrelations => gmera3[:correlations],
        :gmera3rdapprox => gmera3approx,
        :gmera4themptyisometry => gmera4[:emptyisometry],
        :gmera4thfilledisometry => gmera4[:filledisometry],
        :gmera4thcourierisometry => gmera4[:courierisometry],
        :gmera4thapprox => gmera4approx,
        :correlations => gmera4[:correlations],
        :transistionmap => transistionmap)
end

rg1size48 = gmera(blockedcorrelations48,idmap(blockedcorrelations48|>getinspace))
# rgblockedcorrelations = gmera(blockedcorrelations48,idmap(blockedcorrelations48|>getinspace))
rgblockedcorrelations|>getoutspace|>unitcellfock
visualize(regioncorrelations(rgblockedcorrelations,RegionFock(rgblockedcorrelations|>getoutspace|>unitcellfock))|>eigspech)
rg2size48 = gmera(rg1size48[:correlations],rg1size48[:transistionmap])
rg3size48 = gmera(rg2size48[:correlations],rg2size48[:transistionmap])

core48 = @time distillation(rg3size48[:correlations]|>crystalspectrum, :filled=> v -> v < 1e-5, :empty => v -> v > 1e-5)
coreemptyprojector48 = core48[:empty]|>crystalprojector
rg3size48[:transistionmap]*coreemptyprojector48*rg3size48[:transistionmap]'


rg1size48approx = rg1size48[:gmera1stapprox]+rg1size48[:gmera2ndapprox]+rg1size48[:gmera3rdapprox]+rg1size48[:gmera4thapprox]
rg2size48approx = rg2size48[:gmera1stapprox]+rg2size48[:gmera2ndapprox]+rg2size48[:gmera3rdapprox]+rg2size48[:gmera4thapprox]
rg3size48approx = rg3size48[:gmera1stapprox]+rg3size48[:gmera2ndapprox]+rg3size48[:gmera3rdapprox]+rg3size48[:gmera4thapprox]

sum(FockMap(abs(blockedcorrelations48-rg1size48approx-rg2size48approx-rg3size48approx-rg3size48[:transistionmap]*coreemptyprojector48*rg3size48[:transistionmap]'))|>rep)/4608

