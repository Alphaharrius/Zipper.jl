using Plotly, SmithNormalForm, LinearAlgebra, OrderedCollections, SparseArrays, Combinatorics
using Zipper

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')
kspace = convert(MomentumSpace, triangular)

pa = triangular & [1/3, 2/3]
pb = triangular & [2/3, 1/3]
pc = (pa + pb) / 2
spatialsnappingcalibration((pa, pb, pc))

c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
c3 = c6^2

unitcell = Subset(pa, pb)
crystal = Crystal(unitcell, [24, 24])
reciprocalhashcalibration(crystal.sizes)

modes::Subset{Mode} = quantize(:pos, unitcell, 1)
m0, m1 = members(modes)

t_a = ComplexF64(-0.3)
t_b = ComplexF64(-0.7)
t_n = ComplexF64(-0.5)
bonds::FockMap = bondmap([
    (m0, m0) => t_a,
    (m1, m1) => t_b,
    (m0, m1) => t_n,
    (m0, setattr(m1, :offset => Point([-1, 0], triangular))) => t_n,
    (m0, setattr(m1, :offset => Point([0, 1], triangular))) => t_n])

energyspectrum = computeenergyspectrum(bonds, crystal=crystal)
energyspectrum |> visualize

groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=1)
groundstates |> visualize
groundstateprojector = groundstates |> crystalprojector

C = idmap(groundstateprojector.outspace) - groundstateprojector
correlationspectrum = C |> crystalspectrum
visualize(correlationspectrum, title="Physical Correlation")

correlations = C

crystalfock = correlations.outspace

scale = Scale([2 0; 0 2])
blockresult = blocking(:action => scale, :correlations => correlations, :crystal => crystalfock |> getcrystal)

blockedcrystal::Crystal = blockresult[:crystal]
blockedcorrelations::FockMap = blockresult[:correlations]

function circularregionmodes(origin::Offset, physicalmodes::Subset{Mode}, radius::Number)::Subset{Mode}
    currentspace::RealSpace = correlations.outspace |> getcrystal |> getspace |> orthogonalspace
    physicalnorm = m -> lineartransform(currentspace, m |> getpos) |> norm
    return filter(p -> physicalnorm(p - origin) < radius, physicalmodes)
end

crystalpoints::Subset{Offset} = latticepoints(blockedcrystal)
blockedmodes::Subset{Mode} = quantize(:pos, blockedcrystal.unitcell, 1)
physicalmodes::Subset{Mode} = spanoffset(blockedmodes, crystalpoints)
scaledtriangular = scale*triangular 

# local corr
localsitemodes::Subset{Mode} = circularregionmodes(scaledtriangular&[-1,-1], physicalmodes, 2.0) 
localregion::Subset{Offset} = Subset(m |> getpos for m in localsitemodes)
visualize(localregion, title="Frozen Seeding Region", visualspace=euclidean(RealSpace, 2))
localfock::FockSpace = FockSpace{Region}(localsitemodes)
localeigspec = regioncorrelations(blockedcorrelations,localfock) |> eigspech 
visualize(localeigspec)
localfrozenmodes = Subset(m for (m,v) in localeigspec |> geteigenvalues if v <0.01 || v>1-0.01)
localcouriermodes = Subset(m for (m,v) in localeigspec |> geteigenvalues if  0.01<v<1-0.01 )
localfrozenbasis = (localeigspec |> geteigenvectors)[:,localfrozenmodes|> FockSpace]
localcourierbasis = (localeigspec |> geteigenvectors)[:,localcouriermodes|> FockSpace]

# courier seed
c6recenter = c6 |> recenter(scaledtriangular&[-1,-1])
c3recenter = c3 |> recenter(scaledtriangular&[-1,-1])


courierseedingcenterA::Offset = (blockedmodes |> getspace) & [-1+2/3, -1+1/3]
courierseedingmodesA::Subset{Mode} = circularregionmodes(courierseedingcenterA, physicalmodes, 0.8)
courierseedingregionA::Subset{Offset} = Subset(m |> getpos for m in courierseedingmodesA)
visualize(courierseedingregionA, title="Frozen Seeding Region", visualspace=euclidean(RealSpace, 2))
courierseedingfockA::FockSpace = FockSpace{Region}(courierseedingmodesA)
courierseedingfockB = FockSpace{Region}(m for m in c6recenter*courierseedingfockA |> getoutspace)

regioncorrelations(blockedcorrelations,courierseedingfockB) |> eigspech |>visualize
courierseedA = findlocalspstates(statecorrelations = blockedcorrelations, regionfock = courierseedingfockA, 
spectrumextractpredicate = v -> 0.2 < v < 0.4, symmetry = identitytransform(2))[2]

courierseedB = findlocalspstates(statecorrelations = blockedcorrelations, regionfock = courierseedingfockB, 
spectrumextractpredicate = v -> 0.6 < v < 0.8, symmetry = identitytransform(2))[2]

courierseedAs = courierseedA + (c3recenter*courierseedA.outspace)*courierseedA*(c3recenter*courierseedA.inspace)' + (c3recenter^2*courierseedA.outspace)*courierseedA*(c3recenter^2*courierseedA.inspace)'
courierseedBs = courierseedB + (c3recenter*courierseedB.outspace)*courierseedB*(c3recenter*courierseedB.inspace)' + (c3recenter^2*courierseedB.outspace)*courierseedB*(c3recenter^2*courierseedB.inspace)'

courierseeds = courierseedAs+courierseedBs

Subset(m |> getpos for m in courierseeds|> getoutspace)|> visualize
localcourierbasis|> getinspace

courierseeds|> getoutspace |> modeattrs
localcourierbasis|> getoutspace |> modeattrs

u,s,vd = svd(localcourierbasis'*courierseeds)
(localcourierbasis*u*vd)[:,1] |> Zipper.columnspec |> visualize
wannierizedlocalcourier = localcourierbasis*u*vd

wannierizedlocalcourier = wannierizedlocalcourier*spatialmap(wannierizedlocalcourier)'

# frozen seed
c6recenter = c6 |> recenter(scaledtriangular&[-1,-1])
c3recenter = c3 |> recenter(scaledtriangular&[-1,-1])


frozenseedingcenterA::Offset = (blockedmodes |> getspace) & [-1+2/3, -1+1/3]
frozenseedingmodesA::Subset{Mode} = circularregionmodes(frozenseedingcenterA, physicalmodes, 0.8)
frozenseedingregionA::Subset{Offset} = Subset(m |> getpos for m in frozenseedingmodesA)
visualize(frozenseedingregionA, title="Frozen Seeding Region", visualspace=euclidean(RealSpace, 2))
frozenseedingfockA::FockSpace = FockSpace{Region}(frozenseedingmodesA)
frozenseedingfockB = FockSpace{Region}(m for m in c6recenter*frozenseedingfockA |> getoutspace)

regioncorrelations(blockedcorrelations,frozenseedingfockA) |> eigspech |>visualize
frozenseedA = findlocalspstates(statecorrelations = blockedcorrelations, regionfock = frozenseedingfockA, 
spectrumextractpredicate = v -> 0.8 < v, symmetry = identitytransform(2))[1]

frozenseedB = findlocalspstates(statecorrelations = blockedcorrelations, regionfock = frozenseedingfockB, 
spectrumextractpredicate = v -> v < 0.2, symmetry = identitytransform(2))[1]

frozenseedAs = frozenseedA + (c3recenter*frozenseedA.outspace)*frozenseedA*(c3recenter*frozenseedA.inspace)' + (c3recenter^2*frozenseedA.outspace)*frozenseedA*(c3recenter^2*frozenseedA.inspace)'
frozenseedBs = frozenseedB + (c3recenter*frozenseedB.outspace)*frozenseedB*(c3recenter*frozenseedB.inspace)' + (c3recenter^2*frozenseedB.outspace)*frozenseedB*(c3recenter^2*frozenseedB.inspace)'

frozenseedingcenter::Offset = (blockedmodes |> getspace) & [-1, -1]
frozenseedingmodes::Subset{Mode} = circularregionmodes(frozenseedingcenter, physicalmodes, 1)
frozenseedingregion::Subset{Offset} = Subset(m |> getpos for m in frozenseedingmodes)
visualize(frozenseedingregion, title="Frozen Seeding Region", visualspace=euclidean(RealSpace, 2))
frozenseedingfock::FockSpace = FockSpace{Region}(frozenseedingmodes)
regioncorrelations(blockedcorrelations,frozenseedingfock) |> eigspech |>visualize

function _findlocalspstates(;
    statecorrelations::FockMap, regionfock::FockSpace,
    symmetry::AffineTransform = identitytransform(statecorrelations |> getcrystal |> dimension),
    spectrumextractpredicate::Function = v -> v < 1e-2,
    linearindependencethreshold::Real = 1E-2,
    degeneracythreshold::Real = 1e-7)

    function lineardependencefilter(spstate::FockMap)::Bool
        crystalspstates::Dict{Momentum, FockMap} = crystalisometries(localisometry=spstate, crystalfock=statecorrelations.outspace)
        crystalspstate::FockMap = directsum(v for (_, v) in crystalspstates)
        pseudoidentity::FockMap = (crystalspstate' * crystalspstate)
        mineigenvalue = minimum(v for (_, v) in pseudoidentity |> eigvalsh)
        return mineigenvalue > linearindependencethreshold
    end

    localcorrelations::FockMap = regioncorrelations(statecorrelations, regionfock)
    localspectrum::EigenSpectrum = eigspech(localcorrelations, groupingthreshold=degeneracythreshold)
    groupeigenvalues::Base.Generator = (
        subset => (localspectrum |> geteigenvalues)[subset |> first]
        for subset in localspectrum |> geteigenvectors |> getinspace |> sparsegrouping(:eigenindex) |> rep)
    selectedgroups = Iterators.filter(p -> p.second |> spectrumextractpredicate, groupeigenvalues)

    selectedisometries = ((localspectrum |> geteigenvectors)[:, group.first |> FockSpace] for group in selectedgroups)
    orthogonalspstates = Iterators.filter(lineardependencefilter, selectedisometries)
    symmetricspstates = (state * symmetricmap(symmetry, state) for state in orthogonalspstates)
    spstates = (state * spatialmap(state)' for state in symmetricspstates)

    return (state |> getinspace |> dimension => state for state in spstates)
end

frozenseedsorig = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = blockedcorrelations, regionfock = frozenseedingfock, 
spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))

frozenseeds = frozenseedsorig+frozenseedAs+frozenseedBs

Subset(m |> getpos for m in frozenseeds|> getoutspace)|> visualize
localfrozenbasis|> getinspace

frozenseeds|> getoutspace |> modeattrs
localfrozenbasis|> getoutspace |> modeattrs

u,s,vd = svd(localfrozenbasis'*frozenseeds)
(localfrozenbasis*u*vd)[:,1] |> Zipper.columnspec |> visualize
wannierizedlocalfrozen= localfrozenbasis*u*vd

wannierizedlocalfrozen = wannierizedlocalfrozen*spatialmap(wannierizedlocalfrozen)'

# local unitary

ft = fourier(blockedcorrelations |> getoutspace, spanoffset(blockedcorrelations |> getinspace |> unitcellfock |> orderedmodes, blockedcrystal |> latticepoints)|> FockSpace)/sqrt(blockedcrystal |> vol)
blockedcorrelationsRS = ft'*blockedcorrelations*ft
localunit1 = (wannierizedlocalcourier+localfrozenbasis)

localunitary = wannierizedlocalfrozen+wannierizedlocalcourier
wannierizedlocalcourier|>getinspace|>modeattrs

localunitaryphysical = ((((blockedcorrelationsRS|> getoutspace) - (localunitary |> getoutspace)) |> idmap) + localunitary)


transformedblockedcorrelationsRS = localunitaryphysical'*blockedcorrelationsRS*localunitaryphysical

function circularfilter(mode::Mode, radius::Real = 1.8)::Bool
    return norm((mode |> getpos |> euclidean))<=radius
end

localreuniregion::Subset{Offset} = Subset(m |> getpos for m in localunitary |> getinspace)

localrefmodes::Subset{Mode} = Subset(m for m in transformedblockedcorrelationsRS |> getoutspace |> orderedmodes if circularfilter(m))

#localrefmodes::Subset{Mode} = Iterators.filter(circularfilter,m for m in transformedblockedcorrelationsRS |> getoutspace) |> Subset
localrefregion::Subset{Offset} = Subset(m |> getpos for m in localrefmodes)
visualize(localrefregion, title="l", visualspace=euclidean(RealSpace, 2))

visualize(eigspech(transformedblockedcorrelationsRS[localrefmodes |> FockSpace, localrefmodes |> FockSpace]))


localregion::Subset{Offset} = Subset(m |> getpos for m in localrefmodes)
visualize(localregion, title="Frozen Seeding Region", visualspace=euclidean(RealSpace, 2))


globaldistiller = globaldistillerhamiltonian(
    correlations=blockresult[:correlations],
    restrictspace=frozenseedingfock,
    localisometryselectionstrategy=frozenselectionbycount(6),
    symmetry=c3)

globaldistillerspectrum = globaldistiller |> crystalspectrum
globaldistillerspectrum |> visualize

distillresult = distillation(globaldistillerspectrum, :courier => v -> abs(v) < 1e-5, :filled => v -> v > 1e-5, :empty => v -> v < -1e-5)

filledseedingcenter::Offset = (blockedmodes |> getspace) & [2/3, 1/3]
courierseedingmodes::Subset{Mode} = circularregionmodes(courierseedingcenter, physicalmodes, 1.8)
courierseedingregion::Subset{Offset} = Subset(m |> pos for m in courierseedingmodes)
visualize(courierseedingregion, courierseedingcenter |> Subset, title="Courier Seeding Region", visualspace=euclidean(RealSpace, 2))
courierseedingfock::FockSpace = FockSpace(courierseedingregion)

c3 = c6^2 |> recenter(courierseedingcenter)

blockedcourierprojector = distillresult[:courier] |> crystalprojector
blockedcouriercorrelation = idmap(blockedcourierprojector.outspace, blockedcourierprojector.outspace) - blockedcourierprojector

localcourierseed = [regionalwannierseeding(blockedcouriercorrelation, courierseedingfock, symmetry=c3)...][1]
fullcourierseed = localcourierseed + (c6 * localcourierseed.outspace) * localcourierseed * (c6 * localcourierseed.inspace)'

crystalcourierseeds = crystalisometries(localisometry=fullcourierseed, crystalfock=blockedcorrelations.outspace, addinspacemomentuminfo=true)

wanniercourierisometry = wannierprojection(
    crystalisometries=distillresult[:courier].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalcourierseeds)

couriercorrelations = wanniercourierisometry' * blockedcorrelations * wanniercourierisometry
couriercorrelations |> crystalspectrum |> visualize
couriercorrelationspectrum = couriercorrelations |> crystalspectrum
purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
couriercorrelations = purifiedcorrelationspectrum |> FockMap

commutation(c6 * couriercorrelations.outspace, couriercorrelations) |> maximum

