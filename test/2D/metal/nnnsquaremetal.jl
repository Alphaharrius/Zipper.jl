using LinearAlgebra
using Zipper, Plots
plotlyjs()

setmaxthreads(Threads.nthreads())

square = euclidean(RealSpace, 2)
point = [1/2, 1/2] ∈ square
spatialsnappingcalibration([point])

kspace = convert(MomentumSpace, square)

c4 = pointgrouptransform([0 -1; 1 0], localspace=square)
c2 = pointgrouptransform([-1 0; 0 -1], localspace=square)
m135 = pointgrouptransform([0 -1; -1 0], localspace=square)
m45 = pointgrouptransform([0 1; 1 0], localspace=square)

unitcell = Subset(point)
crystal = Crystal(unitcell, [128, 128])
reciprocalhashcalibration(crystal.sizes)

m = quantize(unitcell, 1)|>first

t_n = ComplexF64(-1.)
t_nn = ComplexF64(0.1)
bonds::FockMap = bondmap([
    (m, m|>setattr(:r=>[1, 0]∈square))=>t_n,
    (m, m|>setattr(:r=>[0, 1]∈square))=>t_n,
    (m, m|>setattr(:r=>[1, 1]∈square))=>t_nn,
    (m, m|>setattr(:r=>[1, -1]∈square))=>t_nn,])

@info "Computing energy spectrum..."
energyspectrum = @time computeenergyspectrum(bonds, crystal=crystal)
@info "Resolving ground states..."
groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=0.5)
@info "Computing ground state correlations..."
groundstateprojector = groundstates |> crystalprojector
correlations = idmap(groundstateprojector|>getoutspace) - groundstateprojector
@info "Performing initial blocking..."
initscale = Scale([2 0; 0 2], square)
initblock = initscale * (correlations|>getoutspace)
inputcorrelations = initblock * correlations * initblock'

function zer0d(correlations)
    @info "Starting 0D ZER..."
    @info "Performing blocking..."
    scale = Scale([2 0; 0 2], square)
    block = scale * (correlations|>getoutspace)
    blockedcorrelations = block * correlations * block'
    blockedcrystalfock = blockedcorrelations|>getoutspace
    blockedcrystal = blockedcrystalfock|>getcrystal
    blockedspace = blockedcrystal|>getspace

    @info "Searching for frozen states..."
    region = blockedcrystal|>getunitcell
    regionfock = getregionfock(blockedcrystalfock, region)
    transform = fourier(blockedcrystalfock, regionfock)
    localcorrelations = transform'*blockedcorrelations*transform/(blockedcrystal|>vol)
    localstates = getregionstates(localcorrelations=localcorrelations, grouping=[1, 14, 1])
    frozenlocalstates = localstates[1] + localstates[3]
    region = region .- blockedspace*(0.5, 0.5)
    regionfock = getregionfock(blockedcrystalfock, region)
    transform = fourier(blockedcrystalfock, regionfock)
    localcorrelations = transform'*blockedcorrelations*transform/(blockedcrystal|>vol)
    localstates = getregionstates(localcorrelations=localcorrelations, grouping=[1, 14, 1])
    frozenlocalstates += localstates[1] + localstates[3]

    @info "Computing frozen projector..."
    localfrozen = frozenlocalstates|>FockMap
    transform = fourier(blockedcrystalfock, localfrozen|>getoutspace)
    crystalfrozen = transform * localfrozen
    frozenprojector = crystalfrozen .* crystalfrozen'

    @info "Distilling..."
    bands = groupbands(frozenprojector|>crystalspectrum, :frozen=>(v -> v > 0.1))

    @info "Computing metallic correlations..."
    metalicbands = bands[:others]
    metalprojector = metalicbands|>crystalprojector
    metalcorrelations = 1 - metalprojector

    @info "Wannierizing metallic bands..."
    metalcrystalfock = metalcorrelations|>getoutspace
    metalcrystal = metalcrystalfock|>getcrystal
    metalspace = metalcrystal|>getspace

    @info "Searching metallic local seeds..."
    region = getsphericalregion(crystal=metalcrystal, radius=0.25, metricspace=metalspace) .+ metalspace*(0.25, 0.25)
    regionfock = getregionfock(metalcrystalfock, region)
    restrict = fourier(blockedcrystalfock, regionfock)
    localcorrelations = restrict'*metalcorrelations*restrict/(metalcrystal|>vol)
    localstates = getregionstates(localcorrelations=localcorrelations, grouping=[3])[1]|>normalize

    seedstates = sum(g.*localstates for g in c4|>recenter(metalspace*(0.5, 0.5))|>pointgroupelements)
    localseeds = seedstates|>FockMap
    transform = fourier(metalcorrelations|>getoutspace, localseeds|>getoutspace)
    crystalseeds = transform * localseeds
    crystalseeds = Dict(k=>crystalseeds[k, :] for k in metalcrystal|>brillouinzone)
    pseudoidens = (u' * u for (_, u) in crystalseeds)
    lineardepmetric = (v for id in pseudoidens for (_, v) in id|>eigvalsh)|>minimum
    @info "Local seeds linear-dependence metric: $lineardepmetric"

    @info "Wannierizing metal bands..."
    wanniermetalisometry = wannierprojection(
        crystalisometries=metalicbands|>geteigenvectors, crystal=metalcrystal, crystalseeds=crystalseeds)
    visualregion = getsphericalregion(crystal=metalcrystal, radius=2, metricspace=metalspace) .+ metalspace*(0.5, 0.5)
    visualfock = getregionfock(metalcrystalfock, visualregion)
    leftrestrict = fourier(metalcrystalfock, visualfock)
    rightfock = wanniermetalisometry|>getinspace|>unitcellfock|>RegionFock
    rightrestrict = fourier(wanniermetalisometry|>getinspace, rightfock)
    metallocalstates = leftrestrict' * wanniermetalisometry * rightrestrict

    couriercorrelations = wanniermetalisometry'*blockedcorrelations*wanniermetalisometry
    couriercorrelations = roundingpurification(couriercorrelations|>crystalspectrum)|>CrystalFockMap

    return Dict(
        :couriercorrelations=>couriercorrelations, 
        :wanniermetalisometry=>wanniermetalisometry,
        :courierlocalstates=>metallocalstates)
end

function stripdistillation(couriercorrelations, scaledspace, normalvector)
    @info "Starting strip distillation..."
    @info normalvector
    couriercrystalfock = couriercorrelations|>getoutspace
    couriercrystal = couriercrystalfock|>getcrystal
    courierspace = couriercrystal|>getspace

    @info "Performing extended restrictions..."
    extendedscale = Scale(scaledspace|>rep, courierspace)
    extendedrestrict = ExtendedRestrict(extendedscale, normalvector, 0.5)
    restrict = extendedrestrict * couriercrystalfock
    stripcorrelations = restrict * couriercorrelations * restrict'
    stripspectrum = stripcorrelations|>crystalspectrum

    @info "Computing strip frozen correlations..."
    threshold = 0.05
    stripfrozenstates = groupbands(stripspectrum, :frozen=>(v -> v < threshold || v > 1-threshold))[:frozen]
    stripfrozenprojector = stripfrozenstates|>crystalprojector
    stripfrozencorrelations = idmap(stripfrozenprojector|>getoutspace) - stripfrozenprojector

    @info "Computing strip truncation restricted Fourier transform..."
    stripunitcell = stripspectrum|>getcrystal|>getunitcell
    truncateregion = sum(stripunitcell.+normalvector*n for n in 0:2)
    truncatefock = getregionfock(stripfrozencorrelations|>getoutspace, truncateregion)
    restrict = fourier(stripfrozencorrelations|>getoutspace, truncatefock)
    localcorrelations = restrict'*stripfrozencorrelations*restrict/(stripfrozencorrelations|>getoutspace|>getcrystal|>vol)
    localindices = Iterators.product(localcorrelations|>getoutspace, localcorrelations|>getoutspace)
    pruningindices = Dict((getattr(to, :r) - getattr(from, :r), from|>removeattr(:r), to|>removeattr(:r))=>(from, to) for (from, to) in localindices)
    pruningindices = (index for (_, index) in pruningindices)
    prunedcorrelations = extractindices(localcorrelations, pruningindices)
    truncatedstripfrozencorrelations = restrict * prunedcorrelations .* restrict'
    truncatedspectrum = truncatedstripfrozencorrelations|>crystalspectrum

    @info "Extracting strip filled states..."
    stripfilledstates = groundstatespectrum(truncatedspectrum, perunitcellfillings=4)
    stripfilledprojector = stripfilledstates|>crystalprojector
    stripfilledcorrelations = idmap(stripfilledprojector|>getoutspace) - stripfilledprojector

    @info "Searching strip filled seeds..."
    scaledcrystal = stripfilledstates|>getcrystal

    region = getcrosssection(crystal=couriercrystal, normalvector=normalvector, radius=0.5, minbottomheight=0.5)
    regionfock = getregionfock(stripfilledcorrelations|>getoutspace, region)
    restrict = fourier(stripfilledcorrelations|>getoutspace, regionfock)
    localcorrelations = restrict'*stripfilledcorrelations*restrict/(scaledcrystal|>vol)
    localcorrelations|>eigspech|>visualize
    localstates = getregionstates(localcorrelations=localcorrelations, grouping=[2])[1]
    localstates = m45 * localstates
    localstates = m135 * localstates

    seedstates = localstates

    region = region .- normalvector*0.5
    regionfock = getregionfock(stripfilledcorrelations|>getoutspace, region)
    restrict = fourier(stripfilledcorrelations|>getoutspace, regionfock)
    localcorrelations = restrict'*stripfilledcorrelations*restrict/(scaledcrystal|>vol)
    localcorrelations|>eigspech|>visualize
    localstates = getregionstates(localcorrelations=localcorrelations, grouping=[2])[1]
    localstates = m45 * localstates
    localstates = m135 * localstates

    seedstates += localstates

    localseeds = seedstates|>FockMap
    transform = fourier(stripfilledcorrelations|>getoutspace, localseeds|>getoutspace) / (scaledcrystal|>vol|>sqrt)
    crystalseeds = transform * localseeds
    crystalseeds = Dict(k=>crystalseeds[k, :] for k in scaledcrystal|>brillouinzone)
    pseudoidens = (u' * u for (_, u) in crystalseeds)
    lineardepmetric = (v for id in pseudoidens for (_, v) in id|>eigvalsh if v < 0.2)|>minimum
    @info "Local seeds linear-dependence metric: $lineardepmetric"

    @info "Wannierizing strip filled states..."
    wannierfilledisometry = wannierprojection(
        crystalisometries=stripfilledstates|>geteigenvectors, crystal=scaledcrystal, crystalseeds=crystalseeds)
    wanniercrystal = wannierfilledisometry|>getinspace|>getcrystal
    rightrestrict = (
        fourier(wannierfilledisometry|>getinspace, wannierfilledisometry|>getinspace|>unitcellfock|>RegionFock) / (wanniercrystal|>vol|>sqrt))
    wannierlocalisometry = transform' * wannierfilledisometry * rightrestrict
    stripfilledlocalstates = wannierlocalisometry|>RegionState

    @info "Computing strip filled projector..."
    remapper = spatialremapper(
        wannierlocalisometry|>getoutspace, getsphericalregion(crystal=couriercrystal, radius=2, metricspace=courierspace))
    wannierlocalisometry = remapper * wannierlocalisometry
    transform = fourier(couriercrystalfock, wannierlocalisometry|>getoutspace)
    wanniercrystalisometry = transform * wannierlocalisometry
    wanniermetalprojector = wanniercrystalisometry .* wanniercrystalisometry'

    return Dict(
        :wanniermetalprojector=>wanniermetalprojector, 
        :stripfilledlocalstates=>stripfilledlocalstates)
end

rg0d = zer0d(inputcorrelations)
visualize(rg0d[:courierlocalstates]|>RegionState|>normalize, markersize=5, logscale=0.9)

metalcorrelations = rg0d[:couriercorrelations]
metalcrystal = metalcorrelations|>getoutspace|>getcrystal
metalspace = metalcrystal|>getspace

normalvector = metalspace*(1, 1)
contractbasis = metalspace*(0, metalcrystal|>size|>first)
scaledspace = affinespace(normalvector, contractbasis)

stripdistilled = stripdistillation(metalcorrelations, scaledspace, normalvector)
stripdistiller = stripdistilled[:wanniermetalprojector]
stripdistiller|>crystalspectrum|>visualize
stripbands = groupbands(stripdistiller|>crystalspectrum, :stripmetal=>(v -> v > 0.1))
stripmetalbands = stripbands[:stripmetal]
stripmetalprojector = stripmetalbands|>crystalprojector
stripmetalcorrelations = stripmetalprojector*metalcorrelations*stripmetalprojector
stripmetalcorrelations|>crystalspectrum|>visualize

c4stripdistilled = stripdistillation(metalcorrelations, c4*scaledspace, c4*normalvector)
c4stripdistiller = c4stripdistilled[:wanniermetalprojector]
c4stripbands = groupbands(c4stripdistiller|>crystalspectrum, :stripmetal=>(v -> v > 0.1))
c4stripmetalbands = c4stripbands[:stripmetal]
c4stripmetalprojector = c4stripmetalbands|>crystalprojector
c4stripmetalcorrelations = c4stripmetalprojector*metalcorrelations*c4stripmetalprojector
c4stripmetalcorrelations|>crystalspectrum|>visualize

combineddistiller = stripdistiller + c4stripdistiller
combineddistiller|>crystalspectrum|>visualize
combinedbands = groupbands(combineddistiller|>crystalspectrum, :stripmetal=>(v -> v > 0.01))

courierbands = combinedbands[:others]
courierprojector = courierbands|>crystalprojector
couriercorrelations = 1 - courierprojector

courierprojector*metalcorrelations*courierprojector|>crystalspectrum|>visualize

couriercorrelations|>getoutspace|>unitcellfock|>showmodes

@info "Searching courier seeds..."
couriercrystalfock = couriercorrelations|>getoutspace
couriercrystal = couriercrystalfock|>getcrystal
courierspace = couriercrystal|>getspace
region = courierspace*(0.25, 0.25)|>Subset
regionfock = getregionfock(couriercorrelations|>getoutspace, region)
restrict = fourier(couriercrystalfock, regionfock)
localcorrelations = restrict'*couriercorrelations*restrict/(couriercrystal|>vol)
localcorrelations|>eigspech|>visualize
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[1])[1]
seedstates = sum(g.*localstates for g in c4|>recenter(courierspace*(0.5, 0.5))|>pointgroupelements)

localseeds = seedstates|>FockMap
visualize(localseeds|>RegionState, markersize=5, logscale=0.5)
transform = fourier(couriercrystalfock, localseeds|>getoutspace)
crystalseeds = transform * localseeds
crystalseeds = Dict(k=>crystalseeds[k, :] for k in couriercrystal|>brillouinzone)
pseudoidens = (u' * u for (_, u) in crystalseeds)
lineardepmetric = (v for id in pseudoidens for (_, v) in id|>eigvalsh)|>minimum
@info "Local seeds linear-dependence metric: $lineardepmetric"

@info "Wannierizing courier bands..."
wanniercourierisometry = wannierprojection(
    crystalisometries=courierbands|>geteigenvectors, crystal=couriercrystal, crystalseeds=crystalseeds)
visualregion = getsphericalregion(crystal=couriercrystal, radius=2, metricspace=courierspace) .+ courierspace*(0.5, 0.5)
visualfock = getregionfock(couriercrystalfock, visualregion)
leftrestrict = fourier(couriercrystalfock, visualfock)
rightfock = wanniercourierisometry|>getinspace|>unitcellfock|>RegionFock
rightrestrict = fourier(wanniercourierisometry|>getinspace, rightfock)
courierlocalstates = transform' * wanniercourierisometry * rightrestrict
visualize(courierlocalstates|>RegionState|>normalize, markersize=5, logscale=0.9)

outputcorrelations = wanniercourierisometry'*metalcorrelations*wanniercourierisometry
outputcorrelations = roundingpurification(outputcorrelations|>crystalspectrum, tolerance=0.35)|>CrystalFockMap
outputcorrelations|>crystalspectrum|>visualize

correlations = inputcorrelations
@info "Starting regional distillation..."
@info "Performing blocking..."
scale = Scale([2 0; 0 2], square)
block = scale * (correlations|>getoutspace)
blockedcorrelations = block * correlations * block'
blockedcrystalfock = blockedcorrelations|>getoutspace
blockedcrystal = blockedcrystalfock|>getcrystal
blockedspace = blockedcrystal|>getspace

region = blockedcrystal|>getunitcell
regionfock = getregionfock(blockedcrystalfock, region)
transform = fourier(blockedcrystalfock, regionfock)
localcorrelations = transform'*blockedcorrelations*transform/(blockedcrystal|>vol)
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[1, 14, 1])
frozenlocalstates = localstates[1] + localstates[3]
region = region .- blockedspace*(0.5, 0.5)
regionfock = getregionfock(blockedcrystalfock, region)
transform = fourier(blockedcrystalfock, regionfock)
localcorrelations = transform'*blockedcorrelations*transform/(blockedcrystal|>vol)
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[1, 14, 1])
frozenlocalstates += localstates[1] + localstates[3]
localfrozen = frozenlocalstates|>FockMap
transform = fourier(blockedcrystalfock, localfrozen|>getoutspace)
crystalfrozen = transform * localfrozen
frozenprojector = crystalfrozen .* crystalfrozen'
frozenprojector|>crystalspectrum|>visualize

bands = groupbands(frozenprojector|>crystalspectrum, :frozen=>(v -> v > 0.1))

metalicbands = bands[:others]
metalprojector = metalicbands|>crystalprojector
metalcorrelations = 1 - metalprojector
metalcrystalfock = metalcorrelations|>getoutspace
metalcrystal = metalcrystalfock|>getcrystal
metalspace = metalcrystal|>getspace

region = getsphericalregion(crystal=metalcrystal, radius=0.25, metricspace=metalspace) .+ metalspace*(0.25, 0.25)
visualize(metalcrystal|>getunitcell, region)
regionfock = getregionfock(metalcrystalfock, region)
restrict = fourier(blockedcrystalfock, regionfock)
localcorrelations = restrict'*metalcorrelations*restrict/(metalcrystal|>vol)
localcorrelations|>eigspech|>visualize
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[3])[1]|>normalize
visualize(localstates, markersize=5, logscale=0.5)

seedstates = sum(g.*localstates for g in c4|>recenter(metalspace*(0.5, 0.5))|>pointgroupelements)
localseeds = seedstates|>FockMap
visualize(localseeds|>RegionState, markersize=5, logscale=0.5)
transform = fourier(metalcorrelations|>getoutspace, localseeds|>getoutspace)
crystalseeds = transform * localseeds
crystalseeds = Dict(k=>crystalseeds[k, :] for k in metalcrystal|>brillouinzone)
pseudoidens = (u' * u for (_, u) in crystalseeds)
lineardepmetric = (v for id in pseudoidens for (_, v) in id|>eigvalsh)|>minimum
@info "Local seeds linear-dependence metric: $lineardepmetric"

@info "Wannierizing metal bands..."
wanniermetalisometry = wannierprojection(
    crystalisometries=metalicbands|>geteigenvectors, crystal=metalcrystal, crystalseeds=crystalseeds)
visualregion = getsphericalregion(crystal=metalcrystal, radius=2, metricspace=metalspace) .+ metalspace*(0.5, 0.5)
visualfock = getregionfock(metalcrystalfock, visualregion)
leftrestrict = fourier(metalcrystalfock, visualfock)
rightfock = wanniermetalisometry|>getinspace|>unitcellfock|>RegionFock
rightrestrict = fourier(wanniermetalisometry|>getinspace, rightfock)
metallocalstates = leftrestrict' * wanniermetalisometry * rightrestrict
visualize(metallocalstates|>RegionState|>normalize, markersize=5, logscale=0.9)

couriercorrelations = wanniermetalisometry'*blockedcorrelations*wanniermetalisometry
couriercorrelations = couriercorrelations|>crystalspectrum|>roundingpurification|>CrystalFockMap
couriercorrelations|>crystalspectrum|>visualize

couriercorrelations = blockedcorrelations

couriercrystalfock = couriercorrelations|>getoutspace
couriercrystal = couriercrystalfock|>getcrystal
courierspace = couriercrystal|>getspace

@info "Starting strip distillation..."
@info "Performing extended restrictions..."
normalvector = courierspace*(1, 1)
contractbasis = (courierspace*(0, couriercrystal|>size|>first))
extendedscale = Scale(affinespace(normalvector, contractbasis)|>rep, courierspace)
extendedrestrict = ExtendedRestrict(extendedscale, normalvector, 0.5)
restrict = extendedrestrict * couriercrystalfock
stripcorrelations = restrict * couriercorrelations * restrict'
stripspectrum = stripcorrelations|>crystalspectrum
stripspectrum|>linespectrum|>visualize

stripspectrum|>getcrystal|>getunitcell|>visualize

@info "Computing strip frozen correlations..."
threshold = 0.05
stripfrozenstates = groupbands(stripspectrum, :frozen=>(v -> v < threshold || v > 1-threshold))[:frozen]
stripfrozenprojector = stripfrozenstates|>crystalprojector
stripfrozencorrelations = idmap(stripfrozenprojector|>getoutspace) - stripfrozenprojector

@info "Computing strip truncation restricted Fourier transform..."
stripunitcell = stripspectrum|>getcrystal|>getunitcell
truncateregion = sum(stripunitcell.+normalvector*n for n in 0:2)
truncatefock = getregionfock(stripfrozencorrelations|>getoutspace, truncateregion)
restrict = fourier(stripfrozencorrelations|>getoutspace, truncatefock)
localcorrelations = restrict'*stripfrozencorrelations*restrict/(stripfrozencorrelations|>getoutspace|>getcrystal|>vol)
localindices = Iterators.product(localcorrelations|>getoutspace, localcorrelations|>getoutspace)
pruningindices = Dict((getattr(to, :r) - getattr(from, :r), from|>removeattr(:r), to|>removeattr(:r))=>(from, to) for (from, to) in localindices)
pruningindices = (index for (_, index) in pruningindices)
prunedcorrelations = extractindices(localcorrelations, pruningindices)
stripfouriers = (restrict[subspace, :] for (_, subspace) in restrict|>getoutspace|>crystalsubspaces)
stripcrystal = restrict|>getoutspace|>getcrystal
truncatedstripfrozencorrelations = restrict * prunedcorrelations .* restrict'
truncatedspectrum = truncatedstripfrozencorrelations|>crystalspectrum
truncatedspectrum|>linespectrum|>visualize

@info "Extracting strip filled states..."
stripfilledstates = groundstatespectrum(truncatedspectrum, perunitcellfillings=4)
stripfilledprojector = stripfilledstates|>crystalprojector
stripfilledcorrelations = idmap(stripfilledprojector|>getoutspace) - stripfilledprojector

@info "Searching strip filled seeds..."
scaledcrystal = stripfilledstates|>getcrystal
scaledspace = scaledcrystal|>getspace

region = getcrosssection(crystal=couriercrystal, normalvector=normalvector, radius=0.5, minbottomheight=0.5)
visualize(stripcrystal|>getunitcell, region)
regionfock = getregionfock(stripfilledcorrelations|>getoutspace, region)
restrict = fourier(stripfilledcorrelations|>getoutspace, regionfock)
localcorrelations = restrict'*stripfilledcorrelations*restrict/(scaledcrystal|>vol)
localcorrelations|>eigspech|>visualize
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[2])[1]
localstates = m45 * localstates
localstates = m135 * localstates
visualize(localstates|>normalize, markersize=5, logscale=0.5)

seedstates = localstates

region = region .- normalvector*0.5
visualize(stripcrystal|>getunitcell, region)
regionfock = getregionfock(stripfilledcorrelations|>getoutspace, region)
restrict = fourier(stripfilledcorrelations|>getoutspace, regionfock)
localcorrelations = restrict'*stripfilledcorrelations*restrict/(scaledcrystal|>vol)
localcorrelations|>eigspech|>visualize
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[2])[1]
localstates = m45 * localstates
localstates = m135 * localstates
visualize(localstates|>normalize, markersize=5, logscale=0.5)

seedstates += localstates

localseeds = seedstates|>FockMap
visualize(localseeds|>RegionState, markersize=5, logscale=0.5)
transform = fourier(stripfilledcorrelations|>getoutspace, localseeds|>getoutspace) / (scaledcrystal|>vol|>sqrt)
crystalseeds = transform * localseeds
crystalseeds = Dict(k=>crystalseeds[k, :] for k in scaledcrystal|>brillouinzone)
pseudoidens = (u' * u for (_, u) in crystalseeds)
lineardepmetric = (v for id in pseudoidens for (_, v) in id|>eigvalsh if v < 0.2)|>minimum
@info "Local seeds linear-dependence metric: $lineardepmetric"

@info "Wannierizing strip filled states..."
wannierfilledisometry = wannierprojection(
    crystalisometries=stripfilledstates|>geteigenvectors, crystal=scaledcrystal, crystalseeds=crystalseeds)
wanniercrystal = wannierfilledisometry|>getinspace|>getcrystal
rightrestrict = (
    fourier(wannierfilledisometry|>getinspace, wannierfilledisometry|>getinspace|>unitcellfock|>RegionFock) / (wanniercrystal|>vol|>sqrt))
wannierlocalisometry = transform' * wannierfilledisometry * rightrestrict
stripfilledlocalstates = wannierlocalisometry|>RegionState
visualize(stripfilledlocalstates, markersize=5, logscale=0.8)

@info "Computing strip filled projector..."
remapper = spatialremapper(
    wannierlocalisometry|>getoutspace, getsphericalregion(crystal=couriercrystal, radius=2, metricspace=courierspace))
wannierlocalisometry = remapper * wannierlocalisometry
transform = fourier(couriercrystalfock, wannierlocalisometry|>getoutspace)
wanniercrystalisometry = transform * wannierlocalisometry
wanniermetalprojector = wanniercrystalisometry .* wanniercrystalisometry'
wanniermetalprojector|>crystalspectrum|>visualize

c4metalprojector = wanniermetalprojector
bands = groupbands(wanniermetalprojector|>crystalspectrum, :stripmetal=>(v -> v > 0.1))
c4metalbands = bands[:stripmetal]
c4stripmetalprojector = c4metalbands|>crystalprojector

c4stripmetalprojector*couriercorrelations*c4stripmetalprojector|>crystalspectrum|>visualize

globaldistiller = wanniermetalprojector + c4metalprojector

globaldistiller|>crystalspectrum|>visualize
distillerbands = groupbands(globaldistiller|>crystalspectrum, :frozen=>(v -> v > 0.05))
projector = distillerbands[:frozen]|>crystalprojector
projector*couriercorrelations*projector|>crystalspectrum|>visualize

bands = groupbands(wanniermetalprojector|>crystalspectrum, :stripmetal=>(v -> v > 0.1))

stripprojector = bands[:stripmetal]|>crystalprojector
stripprojector*couriercorrelations*stripprojector|>crystalspectrum|>visualize

c4projector = transform * projector * transform'
c4projector*couriercorrelations*c4projector|>crystalspectrum|>visualize

couriercorrelations = 1 - courierprojector

@info "Searching courier seeds..."
wannierregion = getsphericalregion(crystal=blockedcrystal, radius=0.25, metricspace=blockedcrystal|>getspace)
wannierregion = wannierregion .+ ([0.25, 0.25]∈getspace(wannierregion))
regionfock = quantize(wannierregion, 1)
restrict = fourier(blockedcrystalfock, regionfock) / (blockedcrystal|>vol|>sqrt)
localcorrelations = restrict' * couriercorrelations * restrict
regionseeds = getregionstates(localcorrelations=localcorrelations, grouping=[1])[1]
regionseeds = m45 * regionseeds
regionseeds = m135 * regionseeds

regionseedmap = regionseeds|>FockMap
for element in pointgroupelements(c4|>recenter([0.5, 0.5]∈getspace(wannierregion)))
    lefttransform = element * RegionFock(regionseedmap|>getoutspace)
    righttransform = element * RegionFock(regionseedmap|>getinspace)
    regionseedmap += lefttransform * regionseedmap * righttransform'
end
courierseeds = RegionFock(regionseedmap|>getoutspace) * regionseedmap * RegionFock(regionseedmap|>getinspace)

transform = fourier(blockedcrystalfock, courierseeds|>getoutspace) / (blockedcrystal|>vol|>sqrt)
crystalcourierseeds = transform * courierseeds
crystalcourierseeds = Dict(k=>crystalcourierseeds[k, :] for k in blockedcrystal|>brillouinzone)
pseudoidens = (u' * u for (_, u) in crystalcourierseeds)
lineardepmetric = (v for id in pseudoidens for (_, v) in id|>eigvalsh)|>minimum
@info "Courier seeds linear-dependence metric: $lineardepmetric"

@info "Wannierizing courier bands..."
wanniercourierisometry = wannierprojection(
    crystalisometries=courierbands|>geteigenvectors, crystal=blockedcrystal, crystalseeds=crystalcourierseeds)
visualregion = getsphericalregion(crystal=blockedcrystal, radius=2, metricspace=blockedcrystal|>getspace)
visualfock = quantize(visualregion, 1)
leftrestrict = fourier(blockedcrystalfock, visualfock) / (blockedcrystal|>vol|>sqrt)
rightrestrict = fourier(wanniercourierisometry|>getinspace, wanniercourierisometry|>getinspace|>unitcellfock|>RegionFock)
courierlocalstates = leftrestrict' * wanniercourierisometry * rightrestrict
visualize(courierlocalstates|>RegionState|>normalize, markersize=5, logscale=1)

wanniercourierisometry'*blockedcorrelations*wanniercourierisometry|>crystalspectrum|>visualize
resultcorrelations = wanniercourierisometry'*blockedcorrelations*wanniercourierisometry
resultcorrelations = resultcorrelations|>crystalspectrum|>roundingpurification|>CrystalFockMap

stripmetalprojector = bands[:metal]|>crystalprojector
stripmetalcorrelations = 1 - stripmetalprojector

region = getcrosssection(crystal=blockedcrystal, normalvector=normalvector*0.45, radius=0.5, minbottomheight=0.5)
region = region + (region.+normalvector*0.25) + (region.+normalvector*0.5) + (region.+normalvector*0.75)
region|>visualize
regionfock = getregionfock(stripmetalcorrelations|>getoutspace, region)
restrict = fourier(stripmetalcorrelations|>getoutspace, regionfock)
localcorrelations = restrict'*stripmetalcorrelations*restrict/(blockedcrystal|>vol)
localcorrelations|>eigspech|>visualize
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[6])[1]
localstates = m45 * localstates
localstates = m135 * localstates
visualize(localstates|>normalize, markersize=5, logscale=0.5)

seedstates = localstates

region = getcrosssection(crystal=blockedcrystal, normalvector=c4*normalvector*0.45, radius=0.5, minbottomheight=0.5)
region = region + (region.+c4*normalvector*0.25) + (region.+c4*normalvector*0.5) + (region.+c4*normalvector*0.75)
region|>visualize
regionfock = getregionfock(stripmetalcorrelations|>getoutspace, region)
restrict = fourier(stripmetalcorrelations|>getoutspace, regionfock)
localcorrelations = restrict'*stripmetalcorrelations*restrict/(blockedcrystal|>vol)
localcorrelations|>eigspech|>visualize
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[6])[1]
localstates = m45 * localstates
localstates = m135 * localstates
visualize(localstates|>normalize, markersize=5, logscale=0.5)

seedstates += localstates

localseeds = seedstates|>FockMap
transform = fourier(stripmetalcorrelations|>getoutspace, localseeds|>getoutspace)
crystalseeds = transform * localseeds
crystalseeds = Dict(k=>crystalseeds[k, :] for k in blockedcrystal|>brillouinzone)
pseudoidens = (u' * u for (_, u) in crystalseeds)
lineardepmetric = (v for id in pseudoidens for (_, v) in id|>eigvalsh)|>minimum
@info "Local seeds linear-dependence metric: $lineardepmetric"

@info "Wannierizing strip metal bands..."
wanniermetalisometry = wannierprojection(
    crystalisometries=bands[:metal]|>geteigenvectors, crystal=blockedcrystal, crystalseeds=crystalseeds)

rightrestrict = (
    fourier(wanniermetalisometry|>getinspace, wanniermetalisometry|>getinspace|>unitcellfock|>RegionFock))
metallocalstates = transform' * wanniermetalisometry * rightrestrict
visualize(metallocalstates|>RegionState|>normalize, markersize=5, logscale=0.9)

metallocalstates'*metallocalstates|>visualize

wanniermetalprojector = wanniermetalisometry * wanniermetalisometry'
wanniermetalprojector*globaldistiller*wanniermetalprojector|>crystalspectrum|>visualize
c4wanniermetalprojector = c4 * wanniermetalprojector * c4
wanniermetalprojector*blockedcorrelations*wanniermetalprojector|>crystalspectrum|>visualize
wanniermetalprojector|>crystalspectrum|>visualize

wanniermetalisometry'*blockedcorrelations*wanniermetalisometry|>crystalspectrum

wanniermetalprojector*(blockedcorrelations)*wanniermetalprojector|>crystalspectrum|>visualize

ret = wanniermetalisometry'*(frozenprojector)*wanniermetalisometry
ret|>crystalspectrum|>visualize
blk = ret.blocks|>first|>last
nfock = FockSpace(m for fock in sparsegrouping(blk|>getoutspace, :b) for m in fock)
nblk = nfock*blk*nfock
nblk|>getoutspace|>showmodes
nblk|>visualize
wanniermetalprojector*blockedcorrelations*wanniermetalprojector|>crystalspectrum|>visualize

fiodir("/Users/alphaharrius/ZERData/squaremetal/rg1")
rg1distiller = fioload("globaldistiller")
rg1distiller|>crystalspectrum|>visualize
groupbands(rg1distiller|>crystalspectrum, :frozen=>(v -> v > 0.05))

frozendeleter = 1 - frozenprojector
ret = frozendeleter*rg1distiller*frozendeleter|>crystalspectrum
ret|>visualize
groupbands(ret, :frozen=>(v -> v > 0.05))

metalbands = bands[:others]
metalprojector = metalbands|>crystalprojector
metalcorrelations = 1 - metalprojector

templateregion = getsphericalregion(crystal=blockedcrystal, radius=0.25, metricspace=blockedspace)

region = templateregion .+ blockedspace*(0.25, 0.25)
regionfock = getregionfock(metalcorrelations|>getoutspace, region)
restrict = fourier(blockedcrystalfock, regionfock)
localcorrelations = restrict'*metalcorrelations*restrict/(blockedcrystal|>vol)
localcorrelations|>eigspech|>visualize
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[3])[1]|>normalize
visualize(localstates, markersize=5, logscale=0.5)

seedstates = sum(g.*localstates for g in c4|>recenter(blockedspace*(0.5, 0.5))|>pointgroupelements)
seedlocal = seedstates|>FockMap
transform = fourier(metalcorrelations|>getoutspace, seedlocal|>getoutspace)
crystalseeds = transform * seedlocal
crystalseeds = Dict(k=>crystalseeds[k, :] for k in blockedcrystal|>brillouinzone)
pseudoidens = (u' * u for (_, u) in crystalseeds)
lineardepmetric = (v for id in pseudoidens for (_, v) in id|>eigvalsh)|>minimum
@info "Seeds linear-dependence metric: $lineardepmetric"

@info "Wannierizing metal bands..."
wanniermetalisometry = wannierprojection(
    crystalisometries=metalbands|>geteigenvectors, crystal=blockedcrystal, crystalseeds=crystalseeds)
visualregion = getsphericalregion(crystal=blockedcrystal, radius=2, metricspace=blockedspace) .+ blockedspace*(0.5, 0.5)
visualfock = getregionfock(blockedcrystalfock, visualregion)
leftrestrict = fourier(blockedcrystalfock, visualfock)
rightfock = wanniermetalisometry|>getinspace|>unitcellfock|>RegionFock
rightrestrict = fourier(wanniermetalisometry|>getinspace, rightfock)
metallocalstates = leftrestrict' * wanniermetalisometry * rightrestrict
visualize(metallocalstates|>RegionState|>normalize, markersize=5, logscale=0.9)

metalcorrelations = wanniermetalisometry' * blockedcorrelations * wanniermetalisometry
metalcorrelations|>crystalspectrum|>visualize

metalcrystalfock = metalcorrelations|>getoutspace
metalcrystal = metalcrystalfock|>getcrystal
metalspace = metalcrystal|>getspace

@info "Performing extended restrictions..."
normalvector = metalspace*(4, 4)
contractbasis = metalspace*(0, blockedcrystal|>size|>first)
extendedscale = Scale(affinespace(normalvector, contractbasis)|>rep, metalspace)
extendedrestrict = ExtendedRestrict(extendedscale, normalvector, 1)
restrict = extendedrestrict * metalcrystalfock
stripcorrelations = restrict * metalcorrelations * restrict'
stripspectrum = stripcorrelations|>crystalspectrum
stripspectrum|>linespectrum|>visualize

stripspectrum|>getcrystal|>getunitcell|>visualize

@info "Computing strip frozen correlations..."
stripfrozenstates = groupbands(stripspectrum, :frozen=>(v -> v < 1e-3 || v > 1-1e-3))[:frozen]
stripfrozenprojector = stripfrozenstates|>crystalprojector
stripfrozencorrelations = idmap(stripfrozenprojector|>getoutspace) - stripfrozenprojector

@info "Computing strip truncation restricted Fourier transform..."
stripunitcell = getcrosssection(crystal=metalcrystal, normalvector=normalvector*3, radius=1)
truncatedstripfrozencorrelations = truncatetoregion(stripfrozencorrelations, stripunitcell)
truncatedstripfrozencorrelationspectrum = truncatedstripfrozencorrelations|>crystalspectrum

truncatedstripfrozencorrelationspectrum|>linespectrum|>visualize

@info "Extracting strip filled states..."
stripfilledstates = groupbands(truncatedstripfrozencorrelationspectrum, bandgrouping=[6])[1]
stripfilledprojector = stripfilledstates|>crystalprojector
stripfilledcorrelations = idmap(stripfilledprojector|>getoutspace) - stripfilledprojector

stripfilledstates|>linespectrum|>visualize

@info "Searching strip filled seeds..."
scaledcrystal = stripfilledstates|>getcrystal
scaledspace = scaledcrystal|>getspace

wannierregion = getcrosssection(crystal=metalcrystal, normalvector=normalvector, radius=1, minbottomheight=1.25)
visualize(stripfilledcorrelations|>getoutspace|>getcrystal|>getunitcell, wannierregion)
regionfock = getregionfock(stripfilledcorrelations|>getoutspace, wannierregion)
restrict = fourier(stripfilledcorrelations|>getoutspace, regionfock)
localcorrelations = restrict'*stripfilledcorrelations*restrict/(scaledcrystal|>vol)
localcorrelations|>eigspech|>visualize
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[3])[1]
localstates = m45 * localstates
localstates = m135 * localstates
visualize(localstates|>normalize, markersize=5, logscale=0.5)

seedstates = localstates

wannierregion = wannierregion .- normalvector*0.5
regionfock = getregionfock(stripfilledcorrelations|>getoutspace, wannierregion)
restrict = fourier(stripfilledcorrelations|>getoutspace, regionfock)
localcorrelations = restrict'*stripfilledcorrelations*restrict/(scaledcrystal|>vol)
localcorrelations|>eigspech|>visualize
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[3])[1]
localstates = m45 * localstates
localstates = m135 * localstates
visualize(localstates|>normalize, markersize=5, logscale=0.5)

seedstates += localstates

localseeds = seedstates|>FockMap
visualize(localseeds|>RegionState, markersize=5, logscale=0.5)
transform = fourier(stripfilledcorrelations|>getoutspace, localseeds|>getoutspace)
crystalseeds = transform * localseeds
crystalseeds = Dict(k=>crystalseeds[k, :] for k in scaledcrystal|>brillouinzone)
pseudoidens = (u' * u for (_, u) in crystalseeds)
lineardepmetric = (v for id in pseudoidens for (_, v) in id|>eigvalsh)|>minimum
@info "Local seeds linear-dependence metric: $lineardepmetric"

@info "Wannierizing strip filled states..."
wannierfilledisometry = wannierprojection(
    crystalisometries=stripfilledstates|>geteigenvectors, crystal=scaledcrystal, crystalseeds=crystalseeds)
wanniercrystal = wannierfilledisometry|>getinspace|>getcrystal
rightrestrict = (
    fourier(wannierfilledisometry|>getinspace, wannierfilledisometry|>getinspace|>unitcellfock|>RegionFock) / (wanniercrystal|>vol|>sqrt))
wannierlocalisometry = transform' * wannierfilledisometry * rightrestrict
stripfilledlocalstates = wannierlocalisometry|>RegionState
visualize(stripfilledlocalstates, markersize=2, logscale=1)

@info "Computing strip filled projector..."
remapper = spatialremapper(wannierlocalisometry|>getoutspace, metalcrystalfock)
wannierlocalisometry = remapper * wannierlocalisometry
transform = fourier(metalcrystalfock, wannierlocalisometry|>getoutspace)
wanniercrystalisometry = transform * wannierlocalisometry
wanniermetalprojector = wanniercrystalisometry .* wanniercrystalisometry'

wanniermetalprojector|>crystalspectrum|>visualize

@info "Computing global distiller..."
transform = c4 * metalcrystalfock
globaldistiller = wanniermetalprojector + transform * wanniermetalprojector * transform'
globaldistillerspectrum = globaldistiller|>crystalspectrum
globaldistillerspectrum|>visualize

@info "Performing strip distillation..."
bands = groupbands(wanniermetalprojector|>crystalspectrum, :courier=>(v -> v < 1e-3))

courierbands = bands[:others]
courierprojector = courierbands|>crystalprojector
couriercorrelations = 1 - courierprojector

@info "Searching courier seeds..."
wannierregion = getsphericalregion(crystal=blockedcrystal, radius=0.25, metricspace=blockedcrystal|>getspace)
wannierregion = wannierregion .+ ([0.25, 0.25]∈getspace(wannierregion))
regionfock = quantize(wannierregion, 1)
restrict = fourier(blockedcrystalfock, regionfock) / (blockedcrystal|>vol|>sqrt)
localcorrelations = restrict' * couriercorrelations * restrict
regionseeds = getregionstates(localcorrelations=localcorrelations, grouping=[1])[1]

regionseedmap = regionseeds|>FockMap
for element in pointgroupelements(c4|>recenter([0.5, 0.5]∈getspace(wannierregion)))
    lefttransform = element * RegionFock(regionseedmap|>getoutspace)
    righttransform = element * RegionFock(regionseedmap|>getinspace)
    regionseedmap += lefttransform * regionseedmap * righttransform'
end
courierseeds = RegionFock(regionseedmap|>getoutspace) * regionseedmap * RegionFock(regionseedmap|>getinspace)

transform = fourier(blockedcrystalfock, courierseeds|>getoutspace) / (blockedcrystal|>vol|>sqrt)
crystalcourierseeds = transform * courierseeds
crystalcourierseeds = Dict(k=>crystalcourierseeds[k, :] for k in blockedcrystal|>brillouinzone)
pseudoidens = (u' * u for (_, u) in crystalcourierseeds)
lineardepmetric = (v for id in pseudoidens for (_, v) in id|>eigvalsh)|>minimum
@info "Courier seeds linear-dependence metric: $lineardepmetric"

@info "Wannierizing courier bands..."
wanniercourierisometry = wannierprojection(
    crystalisometries=courierbands|>geteigenvectors, crystal=blockedcrystal, crystalseeds=crystalcourierseeds)
visualregion = getsphericalregion(crystal=blockedcrystal, radius=2, metricspace=blockedcrystal|>getspace)
visualfock = quantize(visualregion, 1)
leftrestrict = fourier(blockedcrystalfock, visualfock) / (blockedcrystal|>vol|>sqrt)
rightrestrict = fourier(wanniercourierisometry|>getinspace, wanniercourierisometry|>getinspace|>unitcellfock|>RegionFock)
courierlocalstates = leftrestrict' * wanniercourierisometry * rightrestrict
visualize(courierlocalstates|>RegionState|>normalize, markersize=5, logscale=0.5)

wanniercourierisometry'*blockedcorrelations*wanniercourierisometry|>crystalspectrum|>visualize
resultcorrelations = wanniercourierisometry'*blockedcorrelations*wanniercourierisometry
resultcorrelations = resultcorrelations|>crystalspectrum|>roundingpurification|>CrystalFockMap

@info "Searching courier seeds..."
region = metalspace*(0.25, 0.25)|>Subset
regionfock = getregionfock(metalcrystalfock, region)
restrict = fourier(metalcrystalfock, regionfock)
localcorrelations = restrict'*couriercorrelations*restrict/(metalcrystal|>vol)
localcorrelations|>eigspech|>visualize
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[1])[1]|>normalize

seedstates = sum(g.*localstates for g in c4|>recenter(metalspace*(0.5, 0.5))|>pointgroupelements)
seedlocal = seedstates|>FockMap
transform = fourier(metalcrystalfock, seedlocal|>getoutspace)
crystalseeds = transform * seedlocal
crystalseeds = Dict(k=>crystalseeds[k, :] for k in metalcrystal|>brillouinzone)
pseudoidens = (u' * u for (_, u) in crystalseeds)
lineardepmetric = (v for id in pseudoidens for (_, v) in id|>eigvalsh)|>minimum
@info "Courier seeds linear-dependence metric: $lineardepmetric"

@info "Wannierizing courier bands..."
wanniercourierisometry = wannierprojection(
    crystalisometries=courierbands|>geteigenvectors, crystal=metalcrystal, crystalseeds=crystalseeds)
visualregion = getsphericalregion(crystal=metalcrystal, radius=2, metricspace=metalspace) .+ metalspace*(0.5, 0.5)
visualfock = getregionfock(metalcrystalfock, visualregion)
leftrestrict = fourier(metalcrystalfock, visualfock)
rightrestrict = fourier(wanniercourierisometry|>getinspace, wanniercourierisometry|>getinspace|>unitcellfock|>RegionFock)
courierlocalstates = leftrestrict' * wanniercourierisometry * rightrestrict
visualize(courierlocalstates|>RegionState|>normalize, markersize=5, logscale=0.9)

metalcorrelations|>crystalspectrum|>visualize
purifiedcorrelations = metalcorrelations|>crystalspectrum|>roundingpurification
couriercorrelations = wanniercourierisometry' * metalcorrelations * wanniercourierisometry
couriercorrelations|>crystalspectrum|>visualize
