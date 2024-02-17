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
crystal = Crystal(unitcell, [24, 24])
reciprocalhashcalibration(crystal.sizes)

modes::Subset{Mode} = quantize(unitcell, 1)|>orderedmodes
m0, m1 = members(modes)

tₙ = ComplexF64(-1)
t_a = ComplexF64(-0.1)
t_b = ComplexF64(-0.9)

bonds::FockMap = bondmap([
    (m0, m0) => t_a,
    (m1, m1) => t_b,
    (m0, m1) => tₙ,
    (m0, setattr(m1, :r => Point([-1, 0], triangular))) => tₙ,
    (m0, setattr(m1, :r => Point([0, 1], triangular))) => tₙ])

energyspectrum = @time computeenergyspectrum(bonds, crystal=crystal)
energyspectrum|>visualize

groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=1)
groundstates|>visualize

groundstateprojector = groundstates|>crystalprojector
correlations = idmap(groundstateprojector|>getoutspace) - groundstateprojector

H = CrystalFockMap(energyspectrum)

@info("Starting RG...")
crystalfock = correlations|>getoutspace

scale = Scale([2 0; 0 2], crystalfock|>getcrystal|>getspace)
@info("Performing unitcellblocking...")
@info("Generating unitcellblocking transformation...")
blocker24 = @time scale * crystalfock
@info("Performing unitcellblocking on correlations...")
blockedcorrelations24 = @time blocker24 * correlations * blocker24'
blockedcrystalfock24 = blockedcorrelations24|>getoutspace
blockedcrystal24::Crystal = blockedcrystalfock24|>getcrystal
blockedspace24::RealSpace = blockedcrystal24|>getspace

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
    gmera1 = @time gmerastep(rgblockedcorrelations,rgblockedcorrelations,firstcenterlist, modeselection1stbycountthenbythreshold(1,0.0005))
    # gmera1 = @time gmerastep1(rgblockedcorrelations,firstcenterlist)
    @info ("1st gmera approximation to correlations...")
    gmera1approx = transistionmap*gmera1[:emptyisometry]*gmera1[:emptyisometry]'*transistionmap'
    # gmera1approx = transistionmap*gmera1[:emptyisometry]
    transistionmap = transistionmap*gmera1[:courierisometry]

    @info ("2nd gmera step...")
    gmera2 = @time gmerastep(rgblockedcorrelations,gmera1[:correlations],secondcenterlist,modeselection1stbycountthenbythreshold(1,0.0005))
    @info ("2nd gmera approximation to correlations...")
    gmera2approx = transistionmap*gmera2[:emptyisometry]*gmera2[:emptyisometry]'*transistionmap'
    # gmera2approx = transistionmap*gmera2[:emptyisometry]
    transistionmap = transistionmap*gmera2[:courierisometry]

    @info ("3rd gmera step...")
    gmera3 = @time gmerastep(rgblockedcorrelations,gmera2[:correlations],thirdcenterlist, modeselection1stbycountthenbythreshold(1,0.0005))
    @info ("3rd gmera approximation to correlations...")
    gmera3approx = transistionmap*gmera3[:emptyisometry]*gmera3[:emptyisometry]'*transistionmap'
    # gmera3approx = transistionmap*gmera3[:emptyisometry]
    transistionmap = transistionmap*gmera3[:courierisometry]

    @info ("final gmera step...")
    gmera4 = @time gmerastep(rgblockedcorrelations,gmera3[:correlations],finalcenterlist, modeselection1stbycountthenbythreshold(1,0.0005))
    @info ("4th gmera approximation to correlations...")
    gmera4approx = transistionmap*gmera4[:emptyisometry]*gmera4[:emptyisometry]'*transistionmap'
    # gmera4approx = transistionmap*gmera4[:emptyisometry]
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

# visualize(rg1size24[:gmera1stcorrelations]|>getinspace|>getcrystal|>getunitcell)
rg1size24 = gmera(blockedcorrelations24,idmap(blockedcorrelations24|>getinspace))
rg2size24 = gmera(rg1size24[:correlations],rg1size24[:transistionmap])

rg1size24approx = rg1size24[:gmera1stapprox]+rg1size24[:gmera2ndapprox]+rg1size24[:gmera3rdapprox]+rg1size24[:gmera4thapprox]
rg2size24approx = rg2size24[:gmera1stapprox]+rg2size24[:gmera2ndapprox]+rg2size24[:gmera3rdapprox]+rg2size24[:gmera4thapprox]

core24 = @time distillation(rg2size24[:correlations]|>crystalspectrum, :filled=> v -> v < 1e-5, :empty => v -> v > 1e-5)
coreemptyprojector24 = core24[:empty]|>crystalprojector
coreapprox = rg2size24[:transistionmap]*coreemptyprojector24*rg2size24[:transistionmap]'

errormat = blockedcorrelations24-rg1size24approx-rg2size24approx-coreapprox
visualize((errormat)|>crystalspectrum)
tr(FockMap((errormat)*(errormat)))/1152
evals = [abs(eval) for (mode,eval) in FockMap((errormat))|>eigspech|>geteigenvalues]
sum(evals)/1152

sum(evals)/1152