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
@info "Performing initial blocking..."
initscale = Scale([2 0; 0 2], square)
initblock = initscale * (correlations|>getoutspace)
inputcorrelations = initblock * correlations * initblock'

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

bands = groupbands(frozenprojector|>crystalspectrum, :frozen=>(v -> v > 0.05))
frozenbands = bands[:frozen]
frozenprojector = frozenbands|>crystalprojector

@info "Starting strip distillation..."
@info "Performing extended restrictions..."
normalvector = (blockedcrystal|>getspace)*(1, 1)
contractbasis = (blockedcrystal|>getspace)*(0, blockedcrystal|>size|>first)
extendedscale = Scale(affinespace(normalvector, contractbasis)|>rep, blockedcrystal|>getspace)
extendedrestrict = ExtendedRestrict(extendedscale, normalvector, 0.5)
restrict = extendedrestrict * blockedcrystalfock
stripcorrelations = restrict * blockedcorrelations * restrict'
stripspectrum = stripcorrelations|>crystalspectrum

@info "Computing strip frozen correlations..."
stripfrozenstates = groupbands(stripspectrum, :frozen=>(v -> v < 0.003 || v > 0.99))[:frozen]
stripfrozenprojector = stripfrozenstates|>crystalprojector
stripfrozencorrelations = idmap(stripfrozenprojector|>getoutspace) - stripfrozenprojector

@info "Computing strip truncation restricted Fourier transform..."
stripunitcell = getcrosssection(crystal=blockedcrystal, normalvector=normalvector*3, radius=0.5)
truncatedstripfrozencorrelations = truncatetoregion(stripfrozencorrelations, stripunitcell)
truncatedstripfrozencorrelationspectrum = truncatedstripfrozencorrelations|>crystalspectrum

truncatedstripfrozencorrelationspectrum|>linespectrum|>visualize

@info "Extracting strip filled states..."
stripfilledstates = groundstatespectrum(truncatedstripfrozencorrelationspectrum, perunitcellfillings=8)
stripfilledprojector = stripfilledstates|>crystalprojector
stripfilledcorrelations = idmap(stripfilledprojector|>getoutspace) - stripfilledprojector

@info "Searching strip filled seeds..."
scaledcrystal = stripfilledstates|>getcrystal
scaledspace = scaledcrystal|>getspace

region = getcrosssection(crystal=blockedcrystal, normalvector=normalvector*0.875, radius=0.5, minbottomheight=0.15)
regionfock = getregionfock(stripfilledcorrelations|>getoutspace, region)
restrict = fourier(stripfilledcorrelations|>getoutspace, regionfock)
localcorrelations = restrict'*stripfilledcorrelations*restrict/(scaledcrystal|>vol)
localcorrelations|>eigspech|>visualize
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[3])[1]
localstates = m45 * localstates
localstates = m135 * localstates
visualize(localstates|>normalize, markersize=5, logscale=0.5)

seedstates = localstates

region = region .- normalvector*0.5
regionfock = getregionfock(stripfilledcorrelations|>getoutspace, region)
restrict = fourier(stripfilledcorrelations|>getoutspace, regionfock)
localcorrelations = restrict'*stripfilledcorrelations*restrict/(scaledcrystal|>vol)
localcorrelations|>eigspech|>visualize
localstates = getregionstates(localcorrelations=localcorrelations, grouping=[5])[1]
localstates = m45 * localstates
localstates = m135 * localstates
visualize(localstates|>normalize, markersize=5, logscale=0.5)

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
    wannierlocalisometry|>getoutspace, getsphericalregion(crystal=blockedcrystal, radius=2, metricspace=blockedcrystal|>getspace))
wannierlocalisometry = remapper * wannierlocalisometry
transform = fourier(blockedcrystalfock, wannierlocalisometry|>getoutspace)
wanniercrystalisometry = transform * wannierlocalisometry
wanniermetalprojector = wanniercrystalisometry .* wanniercrystalisometry'

@info "Computing global distiller..."
transform = c4 * blockedcrystalfock
globaldistiller = wanniermetalprojector + transform * wanniermetalprojector * transform'
globaldistillerspectrum = globaldistiller|>crystalspectrum
globaldistillerspectrum|>visualize

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

courierbands = bands[:courier]
courierprojector = courierbands|>crystalprojector
couriercorrelations = 1 - courierprojector

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
