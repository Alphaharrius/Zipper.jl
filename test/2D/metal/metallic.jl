using LinearAlgebra
using Zipper, Plots
plotlyjs()

usecrystaldensemap()
setmaxthreads(Threads.nthreads())

square = euclidean(RealSpace, 2)
point = square*(0, 0)
spatialsnappingcalibration([point])

kspace = convert(MomentumSpace, square)

c4 = pointgrouptransform([0 -1; 1 0], localspace=square)
c2 = pointgrouptransform([-1 0; 0 -1], localspace=square)
m135 = pointgrouptransform([0 -1; -1 0], localspace=square)
m45 = pointgrouptransform([0 1; 1 0], localspace=square)

unitcell = Subset(point)
bc = BoundaryCondition(square, 64, 64)
crystal = Crystal(unitcell, bc)
reciprocalhashcalibration(bc.bounds)

m = quantize(unitcell, 1)|>first

t_n = ComplexF64(-1.)
t_nn = ComplexF64(0.0)
bonds::FockMap = bondmap([
    (m, m|>setattr(:r=>[1, 0]∈square))=>t_n,
    (m, m|>setattr(:r=>[0, 1]∈square))=>t_n,
    (m, m|>setattr(:r=>[1, 1]∈square))=>t_nn,
    (m, m|>setattr(:r=>[1, -1]∈square))=>t_nn,])

@info "Computing energy spectrum..."
energyspectrum = @time computeenergyspectrum(bonds, crystal=crystal)
energyspectrum|>visualize
@info "Resolving ground states..."
groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=0.5)
@info "Computing ground state correlations..."
groundstateprojector = @time groundstates|>crystalprojector
correlations = @time 1 - groundstateprojector
correlations|>crystalspectrum|>visualize
@info "Performing initial blocking..."
initscale = Scale([2 0; 0 2], square)
initblock = @time initscale * (correlations|>getoutspace)
inputcorrelations = @time initblock * correlations * initblock'
inputcorrelations|>crystalspectrum|>visualize

function mainbranch(correlations)
    scale = Scale([2 0; 0 2], correlations|>getoutspace|>getcrystal|>getspace)
    block = scale*(correlations|>getoutspace)
    blockedcorrelations = block*correlations*block'

    blockedcorrelations|>crystalspectrum|>visualize
    blockedcrystalfock = blockedcorrelations|>getoutspace
    blockedcrystal = blockedcrystalfock|>getcrystal
    blockedspace = blockedcrystal|>getspace

    stripregion = getcrosssection(crystal=blockedcrystal, normalvector=blockedspace*(1, 1), radius=0.5) #.+ blockedspace*(0, 1/4)

    normalvector = blockedspace*(1, 1)
    contractvector = blockedspace*(0, getboundsize(blockedcrystal)|>first)
    scaledspace = affinespace(normalvector, contractvector)
    extendedscale = Scale(scaledspace|>rep, blockedcrystal|>getspace)
    extendedrestrict = ExtendedRestrict(extendedscale, normalvector, stripregion)
    restrict = extendedrestrict*blockedcrystalfock
    stripcorrelations = restrict*blockedcorrelations*restrict'
    stripspectrum = stripcorrelations|>crystalspectrum
    stripspectrum|>linespectrum|>visualize

    @info "Computing strip frozen correlations..."
    frozenthreshold = 0.01
    stripfrozenstates = groupbands(stripspectrum, :frozen=>(v -> v < frozenthreshold || v > 1-frozenthreshold))[:frozen]
    stripfrozenprojector = stripfrozenstates|>crystalprojector
    stripfrozencorrelations = 1 - stripfrozenprojector
    stripfrozencorrelations|>crystalspectrum|>linespectrum|>visualize

    @info "Computing strip truncation restricted Fourier transform..."
    stripunitcell = stripspectrum|>getcrystal|>getunitcell
    truncateregion = sum(stripunitcell.+normalvector*n for n in 0:1)
    truncatedcorrelations = truncatetoregion(stripfrozencorrelations, truncateregion)
    truncatedspectrum = truncatedcorrelations|>crystalspectrum
    truncatedspectrum|>linespectrum|>visualize

    @info "Extracting strip filled states..."
    stripfilledstates = groupbands(truncatedspectrum, bandgrouping=[8])|>first
    stripfilledprojector = stripfilledstates|>crystalprojector
    stripfilledcorrelations = 1 - stripfilledprojector

    visualize(blockedcrystal|>getunitcell, stripregion)

    genregion = getcrosssection(crystal=blockedcrystal, normalvector=blockedspace*(1, 1)*0.85, radius=0.5, minbottomheight=0.25)
    genregion = genregion .- getcenter(genregion)
    visualize(blockedcrystal|>getunitcell, genregion)

    visualize(blockedcrystal|>getunitcell, genregion, genregion.+blockedspace*(1/2, 1/2))

    localstates = getlocalstates(stripfilledcorrelations, region=genregion, count=4)
    localstates += getlocalstates(stripfilledcorrelations, region=genregion.+blockedspace*(1/2, 1/2), count=4)

    visualize(localstates, markersize=5, logscale=0.9)

    stripwannier = wannierprojection(stripfilledstates, localstates)
    wannierregion = sum(stripunitcell.+normalvector*n for n in 0:3) .- normalvector*2
    visualize(wannierregion)
    stripstates = getlocalstates(stripwannier, wannierregion)

    striplocal = stripstates|>FockMap
    transform = fourier(blockedcrystalfock, getregionfock(blockedcrystalfock, wannierregion))
    remapper = spatialremapper(striplocal|>getoutspace, blockedcrystalfock)
    striplocal = remapper * striplocal
    stripspanned = transform * striplocal
    stripdistiller = stripspanned .* stripspanned'
    stripdistiller|>crystalspectrum|>visualize

    stripbands = groupbands(stripdistiller|>crystalspectrum, :strip=>(v -> v > 0.01))[:strip]

    stripsdistiller = stripdistiller+c4*stripdistiller*c4
    stripsdistiller|>crystalspectrum|>visualize

    stripsdistilled = groupbands(stripsdistiller|>crystalspectrum, :frozen=>(v -> v > 0.01))

    courierbands = stripsdistilled[:others]
    courierprojector = courierbands|>crystalprojector
    courierprojector*blockedcorrelations*courierprojector|>crystalspectrum|>visualize
    seedingcorrelations = 1 - courierprojector

    @info "Searching seeds..."
    couriercrystalfock = seedingcorrelations|>getoutspace
    couriercrystal = couriercrystalfock|>getcrystal
    courierspace = couriercrystal|>getspace

    region = getsphericalregion(crystal=couriercrystal, radius=0.5, metricspace=courierspace)
    localstates = getlocalstates(seedingcorrelations, region=region, count=1)
    localstates += getlocalstates(seedingcorrelations, region=region.+courierspace*(0.5, 0.5), count=1)
    localstates += getlocalstates(seedingcorrelations, region=region.+courierspace*(0.5, 0), count=1)
    localstates += c4.*localstates[3]
    localstates = localstates|>normalize

    seedstates = localstates
    visualize(seedstates, markersize=5, logscale=0.9)

    courierwannier = wannierprojection(courierbands, seedstates)
    visualregion = getsphericalregion(crystal=couriercrystal, radius=2, metricspace=courierspace) .+ courierspace*(0.5, 0.5)
    courierlocalstates = getlocalstates(courierwannier, visualregion)
    visualize(courierlocalstates, markersize=6, logscale=0.5)

    correlations = courierwannier'*blockedcorrelations*courierwannier
    ee = entanglemententropy(correlations|>crystalspectrum) / (correlations|>getoutspace|>getcrystal|>vol)
    @info "Courier state entanglement entropy: $(ee)"

    correlations|>crystalspectrum|>visualize

    sp = stripbands|>crystalprojector
    c4sp = c4*sp*c4
    stripprojector = (1-c4sp)*sp*(1-c4sp)
    c4stripprojector = (1-sp)*c4sp*(1-sp)
    frozenprojector = (1-courierprojector)*(1-stripprojector)*(1-c4stripprojector)

    seedingcorrelations = 1 - stripprojector
    seedingcrystal = seedingcorrelations|>getoutspace|>getcrystal
    seedingspace = seedingcrystal|>getspace
    region = getsphericalregion(crystal=seedingcrystal, radius=0.5, metricspace=seedingspace)
    localstates = getlocalstates(seedingcorrelations, region=region, count=1)
    localstates += getlocalstates(seedingcorrelations, region=region.+seedingspace*(0.5, 0.5), count=1)
    localstates += getlocalstates(seedingcorrelations, region=region.+seedingspace*(0.5, 0), count=1)
    localstates += m45.*localstates[3]
    localstates = localstates|>normalize
    visualize(localstates, markersize=5, logscale=0.9)

    stripbands = groupbands(stripprojector|>crystalspectrum, :strip=>(v -> v > 0.01))[:strip]
    stripwannier = wannierprojection(stripbands, localstates)
    visualregion = getsphericalregion(crystal=couriercrystal, radius=2, metricspace=courierspace) .+ courierspace*(0.5, 0.5)
    striplocalstates = getlocalstates(stripwannier, visualregion)
    visualize(striplocalstates, markersize=6, logscale=0.9)

    stripcorrelations = stripwannier'*blockedcorrelations*stripwannier
    stripcorrelations|>crystalspectrum|>visualize

    seedingcorrelations = 1 - c4stripprojector
    seedingcrystal = seedingcorrelations|>getoutspace|>getcrystal
    seedingspace = seedingcrystal|>getspace
    region = getsphericalregion(crystal=seedingcrystal, radius=0.5, metricspace=seedingspace)
    localstates = getlocalstates(seedingcorrelations, region=region, count=1)
    localstates += getlocalstates(seedingcorrelations, region=region.+seedingspace*(-0.5, 0.5), count=1)
    localstates += getlocalstates(seedingcorrelations, region=region.+seedingspace*(0, 0.5), count=1)
    localstates += m135.*localstates[3]
    visualize(localstates, markersize=5, logscale=0.9)

    c4stripbands = groupbands(c4stripprojector|>crystalspectrum, :strip=>(v -> v > 0.01))[:strip]
    c4stripwannier = wannierprojection(c4stripbands, localstates)
    visualregion = getsphericalregion(crystal=seedingcrystal, radius=2, metricspace=seedingspace) .+ seedingspace*(-0.5, 0.5)
    c4striplocalstates = getlocalstates(c4stripwannier, visualregion)
    visualize(c4striplocalstates, markersize=6, logscale=0.9)

    c4stripcorrelations = c4stripwannier'*blockedcorrelations*c4stripwannier
    c4stripcorrelations|>crystalspectrum|>visualize

    seedingcorrelations = 1 - frozenprojector
    seedingcrystal = seedingcorrelations|>getoutspace|>getcrystal
    seedingspace = seedingcrystal|>getspace
    region = getsphericalregion(crystal=seedingcrystal, radius=0.5, metricspace=seedingspace)
    localstates = getlocalstates(seedingcorrelations, region=region, count=2)
    localstates += getlocalstates(seedingcorrelations, region=region.+seedingspace*(0.5, 0.5), count=2)
    localstates = localstates|>normalize
    visualize(localstates, markersize=5, logscale=0.9)

    frozenbands = groupbands(frozenprojector|>crystalspectrum, :frozen=>(v -> v > 0.01))[:frozen]
    frozenwannier = wannierprojection(frozenbands, localstates)
    visualregion = getsphericalregion(crystal=seedingcrystal, radius=2, metricspace=seedingspace) .+ seedingspace*(0.5, 0.5)
    frozenlocalstates = getlocalstates(frozenwannier, visualregion)
    visualize(frozenlocalstates, markersize=6, logscale=0.9)

    frozencorrelations = frozenwannier'*blockedcorrelations*frozenwannier
    frozencorrelations|>crystalspectrum|>visualize

    return Dict(
        :correlations=>correlations, 
        :stripcorrelations=>stripcorrelations, 
        :c4stripcorrelations=>c4stripcorrelations, 
        :frozencorrelations=>frozencorrelations,
        
        :courierlocalstates=>courierlocalstates,
        :striplocalstates=>striplocalstates,
        :c4striplocalstates=>c4striplocalstates,
        :frozenlocalstates=>frozenlocalstates,
        
        :courierbands=>courierbands,
        :stripbands=>stripbands,
        :c4stripbands=>c4stripbands,
        :frozenbands=>frozenbands)
end

rg1 = mainbranch(inputcorrelations)
inputcorrelations = roundingpurification(rg1[:correlations]|>crystalspectrum)|>crystalfockmap
rg2 = mainbranch(inputcorrelations)
inputcorrelations = roundingpurification(rg2[:correlations]|>crystalspectrum)|>crystalfockmap
rg3 = mainbranch(inputcorrelations)
inputcorrelations = roundingpurification(rg3[:correlations]|>crystalspectrum)|>crystalfockmap

entanglemententropy(rg3[:correlations]|>crystalspectrum) / (rg3[:correlations]|>getoutspace|>dimension)

visualize(rg3[:courierlocalstates], markersize=6, logscale=0.9)
rg3[:correlations]|>crystalspectrum|>visualize
