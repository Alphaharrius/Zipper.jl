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

scale = Scale([3 0; 0 3], crystalfock|>getcrystal|>getspace)
@info("Performing unitcellblocking...")
@info("Generating unitcellblocking transformation...")
blockerloca3size24 = @time scale * crystalfock
@info("Performing unitcellblocking on correlations...")
blockedcorrelationslocal3size24 = @time blockerloca3size24 * correlations * blockerloca3size24'
blockedcrystalfockloca3size24 = blockedcorrelationsloca3size24|>getoutspace
blockedcrystalloca3size24::Crystal = blockedcrystalfockloca3size24|>getcrystal
blockedspaceloca3size24::RealSpace = blockedcrystalloca3size24|>getspace

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
    gmera1 = @time gmerastep(rgblockedcorrelations,rgblockedcorrelations,firstcenterlist, modeselection1stbycountthenbythreshold(1,0.0001))
    # gmera1 = @time gmerastep1(rgblockedcorrelations,firstcenterlist)
    @info ("1st gmera approximation to correlations...")
    gmera1approx = transistionmap*gmera1[:emptyisometry]*gmera1[:emptyisometry]'*transistionmap'
    transistionmap = transistionmap*gmera1[:courierisometry]

    @info ("2nd gmera step...")
    gmera2 = @time gmerastep(rgblockedcorrelations,gmera1[:correlations],secondcenterlist, modeselection1stbycountthenbythreshold(1,0.0001))
    @info ("2nd gmera approximation to correlations...")
    gmera2approx = transistionmap*gmera2[:emptyisometry]*gmera2[:emptyisometry]'*transistionmap'
    transistionmap = transistionmap*gmera2[:courierisometry]

    @info ("3rd gmera step...")
    gmera3 = @time gmerastep(rgblockedcorrelations,gmera2[:correlations],thirdcenterlist, modeselection1stbycountthenbythreshold(1,0.0001))
    @info ("3rd gmera approximation to correlations...")
    gmera3approx = transistionmap*gmera3[:emptyisometry]*gmera3[:emptyisometry]'*transistionmap'
    transistionmap = transistionmap*gmera3[:courierisometry]

    @info ("final gmera step...")
    gmera4 = @time gmerastep(rgblockedcorrelations,gmera3[:correlations],finalcenterlist, modeselection1stbycountthenbythreshold(1,0.0001))
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

rg1loca3size24 = gmera(blockedcorrelationslocal3size24,idmap(blockedcorrelationslocal3size24|>getinspace))
rg2loca3size24 = gmera(rg1loca3size24[:correlations],rg1loca3size24[:transistionmap])
rg2loca3size24[:correlations]

corelocal3size24 = @time distillation(rg2loca3size24[:correlations]|>crystalspectrum, :filled=> v -> v < 1e-5, :empty => v -> v > 1e-5)
coreemptyprojectorlocal3size24 = corelocal3size24[:empty]|>crystalprojector

fullfouriermaplocal3size24 = fullftmap(blockedcorrelationslocal3size24)

rg1local3size24approx = rg1loca3size24[:gmera1stapprox]+rg1loca3size24[:gmera2ndapprox]+rg1loca3size24[:gmera3rdapprox]+rg1loca3size24[:gmera4thapprox]
rg2local3size24approx = rg2loca3size24[:gmera1stapprox]+rg2loca3size24[:gmera2ndapprox]+rg2loca3size24[:gmera3rdapprox]+rg2loca3size24[:gmera4thapprox]

errormatinklocal3size24 = blockedcorrelationslocal3size24-rg1local3size24approx-rg2local3size24approx-rg2loca3size24[:transistionmap]*coreemptyprojectorlocal3size24*rg2loca3size24[:transistionmap]'
errormatlocal3size24 = fullfouriermaplocal3size24'*errormatinklocal3size24*fullfouriermaplocal3size24

# sum(abs(FockMap(errormatinklocal3size24))|>rep)/1152
errorloacl3size24 = sum(abs(FockMap(errormatlocal3size24))|>rep)/1152

errorlocal3size24 = sum(abs(FockMap(errormatlocal3size24))|>rep)/1152