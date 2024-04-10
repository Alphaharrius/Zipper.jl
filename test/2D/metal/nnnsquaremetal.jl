using LinearAlgebra
using Zipper, Plots
plotlyjs()

setmaxthreads(Threads.nthreads())

square = euclidean(RealSpace, 2)
point = square*(1/2, 1/2)
spatialsnappingcalibration([point])

kspace = convert(MomentumSpace, square)

c4 = pointgrouptransform([0 -1; 1 0], localspace=square)
c2 = pointgrouptransform([-1 0; 0 -1], localspace=square)
m135 = pointgrouptransform([0 -1; -1 0], localspace=square)
m45 = pointgrouptransform([0 1; 1 0], localspace=square)

unitcell = Subset(point)
bc = BoundaryCondition(square, 128, 128)
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
@info "Resolving ground states..."
groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=0.5)
@info "Computing ground state correlations..."
groundstateprojector = groundstates |> crystalprojector
correlations = idmap(groundstateprojector|>getoutspace) - groundstateprojector
correlations|>crystalspectrum|>visualize
@info "Performing initial blocking..."
initscale = Scale([2 0; 0 2], square)
initblock = initscale * (correlations|>getoutspace)
inputcorrelations = initblock * correlations * initblock'
inputcorrelations|>crystalspectrum|>visualize

ss = Scale([1 0; 0 2], inputcorrelations|>getoutspace|>getcrystal|>getspace)
bl = ss * getoutspace(inputcorrelations)
bl * inputcorrelations * bl'|>crystalspectrum|>visualize

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
    metallocalstates = leftrestrict'*wanniermetalisometry*rightrestrict|>RegionState|>normalize

    couriercorrelations = wanniermetalisometry'*blockedcorrelations*wanniermetalisometry
    couriercorrelations = roundingpurification(couriercorrelations|>crystalspectrum)|>CrystalFockMap

    return Dict(
        :couriercorrelations=>couriercorrelations, 
        :wanniermetalisometry=>wanniermetalisometry,
        :courierlocalstates=>metallocalstates,
        :block=>block)
end

function stripdistillation(couriercorrelations, scaledspace, normalvector, frozenthreshold)
    @info "Starting strip distillation..."
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
    stripfrozenstates = groupbands(stripspectrum, :frozen=>(v -> v < frozenthreshold || v > 1-frozenthreshold))[:frozen]
    stripfrozenprojector = stripfrozenstates|>crystalprojector
    stripfrozencorrelations = idmap(stripfrozenprojector|>getoutspace) - stripfrozenprojector

    @info "Computing strip truncation restricted Fourier transform..."
    stripunitcell = stripspectrum|>getcrystal|>getunitcell
    truncateregion = sum(stripunitcell.+normalvector*n for n in 0:2)
    truncatedcorrelations = truncatetoregion(stripfrozencorrelations, truncateregion)
    truncatedspectrum = truncatedcorrelations|>crystalspectrum

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
    lineardepmetric = (v for id in pseudoidens for (_, v) in id|>eigvalsh)|>minimum
    @info "Local seeds linear-dependence metric: $lineardepmetric"

    @info "Wannierizing strip filled states..."
    wannierfilledisometry = wannierprojection(
        crystalisometries=stripfilledstates|>geteigenvectors, crystal=scaledcrystal, crystalseeds=crystalseeds)
    wanniercrystal = wannierfilledisometry|>getinspace|>getcrystal
    rightrestrict = (
        fourier(wannierfilledisometry|>getinspace, wannierfilledisometry|>getinspace|>unitcellfock|>RegionFock) / (wanniercrystal|>vol|>sqrt))
    wannierlocalisometry = transform' * wannierfilledisometry * rightrestrict
    stripfilledlocalstates = wannierlocalisometry|>RegionState|>normalize

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

function statewannierization(courierbands, couriercorrelations)
    @info "Searching seeds..."
    couriercrystalfock = couriercorrelations|>getoutspace
    couriercrystal = couriercrystalfock|>getcrystal
    courierspace = couriercrystal|>getspace

    region = courierspace*(0.25, 0.25)|>Subset
    regionfock = getregionfock(couriercorrelations|>getoutspace, region)
    restrict = fourier(couriercrystalfock, regionfock)
    localcorrelations = restrict'*couriercorrelations*restrict/(couriercrystal|>vol)
    localstates = getregionstates(localcorrelations=localcorrelations, grouping=[1])[1]
    localstates = localstates|>FockMap|>abs|>RegionState|>normalize

    seedstates = localstates

    region = courierspace*(0.25, 0.75)|>Subset
    regionfock = getregionfock(couriercorrelations|>getoutspace, region)
    restrict = fourier(couriercrystalfock, regionfock)
    localcorrelations = restrict'*couriercorrelations*restrict/(couriercrystal|>vol)
    localstates = getregionstates(localcorrelations=localcorrelations, grouping=[1])[1]
    localstates = localstates|>FockMap|>abs|>RegionState|>normalize

    seedstates += localstates

    region = courierspace*(0.75, 0.25)|>Subset
    regionfock = getregionfock(couriercorrelations|>getoutspace, region)
    restrict = fourier(couriercrystalfock, regionfock)
    localcorrelations = restrict'*couriercorrelations*restrict/(couriercrystal|>vol)
    localstates = getregionstates(localcorrelations=localcorrelations, grouping=[1])[1]
    localstates = localstates|>FockMap|>abs|>RegionState|>normalize

    seedstates += localstates

    region = courierspace*(0.75, 0.75)|>Subset
    regionfock = getregionfock(couriercorrelations|>getoutspace, region)
    restrict = fourier(couriercrystalfock, regionfock)
    localcorrelations = restrict'*couriercorrelations*restrict/(couriercrystal|>vol)
    localstates = getregionstates(localcorrelations=localcorrelations, grouping=[1])[1]
    localstates = localstates|>FockMap|>abs|>RegionState|>normalize

    seedstates += localstates
    
    localseeds = seedstates|>normalize|>FockMap
    transform = fourier(couriercrystalfock, localseeds|>getoutspace)
    crystalseeds = transform * localseeds
    crystalseeds = Dict(k=>crystalseeds[k, :] for k in couriercrystal|>brillouinzone)
    pseudoidens = (u' * u for (_, u) in crystalseeds)
    lineardepmetric = (v for id in pseudoidens for (_, v) in id|>eigvalsh)|>minimum
    @info "Local seeds linear-dependence metric: $lineardepmetric"

    @info "Wannierizing bands..."
    wanniercourierisometry = wannierprojection(
        crystalisometries=courierbands|>geteigenvectors, crystal=couriercrystal, crystalseeds=crystalseeds)
    visualregion = getsphericalregion(crystal=couriercrystal, radius=4, metricspace=courierspace) .+ courierspace*(0.5, 0.5)
    visualfock = getregionfock(couriercrystalfock, visualregion)
    leftrestrict = fourier(couriercrystalfock, visualfock)
    rightfock = wanniercourierisometry|>getinspace|>unitcellfock|>RegionFock
    rightrestrict = fourier(wanniercourierisometry|>getinspace, rightfock)
    courierlocalstates = leftrestrict'*wanniercourierisometry*rightrestrict|>RegionState|>normalize

    return Dict(:localstates=>courierlocalstates, :courierisometry=>wanniercourierisometry)
end

function renormalize(correlations; frozenthreshold=0.01)
    rg0d = zer0d(correlations)

    metalcorrelations = rg0d[:couriercorrelations]
    metalcrystal = metalcorrelations|>getoutspace|>getcrystal
    metalspace = metalcrystal|>getspace

    normalvector = metalspace*(1, 1)
    contractbasis = metalspace*(0, getboundsize(metalcrystal)|>first)
    @info "Contract basis: $contractbasis"
    scaledspace = affinespace(normalvector, contractbasis)

    stripdistilled = stripdistillation(metalcorrelations, scaledspace, normalvector, frozenthreshold)
    stripdistiller = stripdistilled[:wanniermetalprojector]
    stripbands = groupbands(stripdistiller|>crystalspectrum, :stripmetal=>(v -> v > 0.01))
    stripmetalbands = stripbands[:stripmetal]
    stripmetalprojector = stripmetalbands|>crystalprojector
    stripmetalcorrelations = stripmetalprojector*metalcorrelations*stripmetalprojector

    c4stripdistilled = stripdistillation(metalcorrelations, c4*scaledspace, c4*normalvector, frozenthreshold)
    c4stripdistiller = c4stripdistilled[:wanniermetalprojector]
    c4stripbands = groupbands(c4stripdistiller|>crystalspectrum, :stripmetal=>(v -> v > 0.1))
    c4stripmetalbands = c4stripbands[:stripmetal]
    c4stripmetalprojector = c4stripmetalbands|>crystalprojector
    c4stripmetalcorrelations = c4stripmetalprojector*metalcorrelations*c4stripmetalprojector

    combineddistiller = stripdistiller + c4stripdistiller
    combinedbands = groupbands(combineddistiller|>crystalspectrum, :stripmetal=>(v -> v > 0.01))

    courierbands = combinedbands[:others]
    courierprojector = courierbands|>crystalprojector
    couriercorrelations = 1 - courierprojector

    couriers = statewannierization(courierbands, couriercorrelations)
    courierisometry = couriers[:courierisometry]
    outputcorrelations = courierisometry'*metalcorrelations*courierisometry
    outputcorrelations = roundingpurification(outputcorrelations|>crystalspectrum, tolerance=0.42)|>CrystalFockMap

    stripmetals = statewannierization(stripmetalbands, stripmetalcorrelations)
    stripmetalisometry = stripmetals[:courierisometry]
    outputstripcorrelations = stripmetalisometry'*metalcorrelations*stripmetalisometry
    outputstripcorrelations = roundingpurification(outputstripcorrelations|>crystalspectrum, tolerance=0.4)|>CrystalFockMap

    c4stripmetals = statewannierization(c4stripmetalbands, c4stripmetalcorrelations)
    c4stripmetalisometry = c4stripmetals[:courierisometry]
    outputc4stripcorrelations = c4stripmetalisometry'*metalcorrelations*c4stripmetalisometry
    outputc4stripcorrelations = roundingpurification(outputc4stripcorrelations|>crystalspectrum, tolerance=0.4)|>CrystalFockMap

    return Dict(
        :block=>rg0d[:block],
        :metalcorrelations=>metalcorrelations,
        :metallicstates=>rg0d[:courierlocalstates],
        :metallicisometry=>rg0d[:wanniermetalisometry],
        :stripdistiller=>stripdistiller,
        :stripmetalprojector=>stripmetalprojector,
        :stripmetalcorrelations=>outputstripcorrelations,
        :stripmetalstates=>stripdistilled[:stripfilledlocalstates],
        :c4stripdistiller=>c4stripdistiller,
        :c4stripmetalprojector=>c4stripmetalprojector,
        :c4stripmetalcorrelations=>outputc4stripcorrelations,
        :globaldistiller=>combineddistiller,
        :couriercorrelations=>outputcorrelations, 
        :courierlocalstates=>couriers[:localstates],
        :courierisometry=>courierisometry,)
end

fiodir("/Users/alphaharrius/ZERData/squaremetal")
fiosave(energyspectrum|>CrystalFockMap, name="hamiltonian")
fiosave(correlations, name="groundstatecorrelations")
fiosave(initblock, name="initblock")
fiosave(inputcorrelations, name="inputcorrelations")
fioload("inputcorrelations")|>crystalspectrum|>visualize

rg1 = renormalize(inputcorrelations, frozenthreshold=0.01)
rg1[:couriercorrelations]|>crystalspectrum|>visualize
rg1[:c4stripmetalcorrelations]|>crystalspectrum|>visualize
rg1[:globaldistiller]|>crystalspectrum|>visualize
roundingpurification(rg1[:stripmetalcorrelations]|>crystalspectrum, tolerance=0.4)|>visualize
visualize(rg1[:courierlocalstates], markersize=5, logscale=1)
visualize(rg1[:metallicstates], markersize=5, logscale=1)
visualize(rg1[:stripmetalstates], markersize=5, logscale=1)

fiodir("/Users/alphaharrius/ZERData/squaremetal/RG1")
for (k, v) in rg1
    fiosave(v, name=k|>String)
end

rg2 = renormalize(rg1[:couriercorrelations], frozenthreshold=0.01)
rg2[:couriercorrelations]|>crystalspectrum|>visualize
rg2[:globaldistiller]|>crystalspectrum|>visualize
visualize(rg2[:courierlocalstates], markersize=5, logscale=1)
visualize(rg2[:metallicstates], markersize=5, logscale=1)

rg3 = renormalize(rg2[:couriercorrelations], frozenthreshold=0.01)
rg3[:couriercorrelations]|>crystalspectrum|>visualize
rg3[:globaldistiller]|>crystalspectrum|>visualize
visualize(rg3[:metallicstates], markersize=5, logscale=1)
visualize(rg3[:courierlocalstates], markersize=5, logscale=1)

rg4 = renormalize(rg3[:couriercorrelations], frozenthreshold=0.01)
rg4[:couriercorrelations]|>crystalspectrum|>visualize
rg4[:globaldistiller]|>crystalspectrum|>visualize
visualize(rg4[:courierlocalstates], markersize=5, logscale=1)
visualize(rg4[:metallicstates], markersize=5, logscale=1)
