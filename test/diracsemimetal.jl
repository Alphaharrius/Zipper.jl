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

modes::Subset{Mode} = quantize(:b, unitcell, 1)
m0, m1 = members(modes)

tₙ = ComplexF64(-1.)
bonds::FockMap = bondmap([
    (m0, m1) => tₙ,
    (m0, m1 |> setattr(:r => triangular & [-1, 0])) => tₙ,
    (m0, m1 |> setattr(:r => triangular & [0, 1])) => tₙ])

energyspectrum = computeenergyspectrum(bonds, crystal=crystal)
energyspectrum |> visualize

groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=1)
groundstates |> visualize
groundstateprojector = groundstates |> crystalprojector

C = idmap(groundstateprojector.outspace) - groundstateprojector

correlations = C

crystalfock = correlations.outspace

scale = Scale([2 0; 0 2])
blockresult = blocking(:action => scale, :correlations => correlations, :crystal => crystalfock |> getcrystal)

blockedcrystal::Crystal = blockresult[:crystal]
blockedcorrelations::FockMap = blockresult[:correlations]


frozenseedingmodes::Subset{Mode} = circularregionmodes(triangular |> origin, physicalmodes, 2.0) 
frozenseedingregion::Subset{Offset} = Subset(m |> pos for m in frozenseedingmodes)
visualize(frozenseedingregion, title="Frozen Seeding Region", visualspace=euclidean(RealSpace, 2))
frozenseedingfock::FockSpace = FockSpace(frozenseedingmodes)
localevec = regioncorrelations(blockedcorrelations,frozenseedingfock) |> eigvecsh 
filledmodes = (localevec |> getinspace |> orderedmodes |> collect)[1:3]
emptymodes = (localevec |> getinspace |> orderedmodes |> collect)[22:24]
filledmodes |> Subset |> FockSpace 

fillediso = columns(localevec, filledmodes |> Subset |> FockSpace ) 
filledproj = fillediso * fillediso'

emptyiso = columns(localevec, emptymodes |> Subset |> FockSpace ) 
emptyproj = emptyiso * emptyiso'
locproj = filledproj+emptyproj
couriermodes = [((locproj |> eigspech) |> groupbyeigenvalues)...][1].second
courieriso = columns(locproj |> eigvecsh, couriermodes|> FockSpace )

courierseedingcenter::Offset = (blockedmodes |> getspace) & [2/3, 1/3]
courierseedingmodes::Subset{Mode} = circularregionmodes(courierseedingcenter, physicalmodes, 0.8)
courierseedingregion::Subset{Offset} = Subset(m |> pos for m in courierseedingmodes)
visualize(courierseedingregion, title="Frozen Seeding Region", visualspace=euclidean(RealSpace, 2))
courierseedingfock::FockSpace = FockSpace(courierseedingmodes)
locseedregionevec = regioncorrelations(blockedcorrelations,courierseedingfock) |> eigvecsh 
locseedvec = (locseedregionevec |> getinspace |> orderedmodes |> collect)[2:4]

localseed = columns(locseedregionevec, localseedvec |> Subset |> FockSpace )

# columns(locseedregionevec,(localseed |> getinspace |> orderedmodes |> collect)[3] |> Subset |> FockSpace ) |> Quantum.columnspec |> visualize

c3rep = c3 * (localseed |> getoutspace)
c3localseed = c3rep * localseed
c3localseed = FockMap(c3localseed,inspace=c3localseed |> getinspace |> orderedmodes |> setattr(:name=> "c3a")|> FockSpace,permute=false)
c3c3rep = c3^2 * (localseed |> getoutspace)
c3c3localseed = c3c3rep * localseed
c3c3localseed = FockMap(c3c3localseed,inspace=c3c3localseed |> getinspace |> orderedmodes |> setattr(:name=> "c3c3a")|> FockSpace,permute=false)
aseed = (localseed+c3localseed+c3c3localseed)

courierseedingcenter::Offset = (blockedmodes |> getspace) & [1/3, 2/3]
courierseedingmodes::Subset{Mode} = circularregionmodes(courierseedingcenter, physicalmodes, 0.8)
courierseedingregion::Subset{Offset} = Subset(m |> pos for m in courierseedingmodes)
visualize(courierseedingregion, title="Frozen Seeding Region", visualspace=euclidean(RealSpace, 2))
courierseedingfock::FockSpace = FockSpace(courierseedingmodes)
locseedregionevec = regioncorrelations(blockedcorrelations,courierseedingfock) |> eigvecsh 
[(regioncorrelations(blockedcorrelations,courierseedingfock) |> eigvalsh)...]
localseedvec = (locseedregionevec |> getinspace |> orderedmodes |> collect)[1:3]
localseed = columns(locseedregionevec, localseedvec |> Subset |> FockSpace )
localseed = FockMap(localseed,inspace=localseed |> getinspace |> orderedmodes |> setattr(:name=> "b")|> FockSpace,permute=false)

c3rep = c3 * (localseed |> getoutspace)
c3localseed = c3rep * localseed
c3localseed = FockMap(c3localseed,inspace=c3localseed |> getinspace |> orderedmodes |> setattr(:name=> "c3b")|> FockSpace,permute=false)
c3c3rep = c3^2 * (localseed |> getoutspace)
c3c3localseed = c3c3rep * localseed
c3c3localseed = FockMap(c3c3localseed,inspace=c3c3localseed |> getinspace |> orderedmodes |> setattr(:name=> "c3c3b")|> FockSpace,permute=false)
bseed = (localseed+c3localseed+c3c3localseed)
aseed

u,sval, vd = (aseed+bseed)'*courieriso |> svd
u
vd
[sval...]
u*vd
(u*vd) |> getinspace |> modeattrs
courieriso' |> getinspace |> modeattrs
wanniercourieriso = courieriso*(u*vd)'

columns(wanniercourieriso, (wanniercourieriso|> getinspace |> orderedmodes |> collect)[6] |> Subset |> FockSpace ) |> Quantum.columnspec |> visualize

[regionalwannierseeding(blockedcorrelations,courierseedingfock,symmetry=c3)...]


function circularregionmodes(origin::Offset, physicalmodes::Subset{Mode}, radius::Number)::Subset{Mode}
    currentspace::RealSpace = correlations.outspace |> getcrystal |> getspace |> orthogonalspace
    physicalnorm = m -> lineartransform(currentspace, m |> getpos) |> norm
    return filter(p -> physicalnorm(p - origin) < radius, physicalmodes)
end

crystalpoints::Subset{Offset} = latticepoints(blockedcrystal)
blockedmodes::Subset{Mode} = quantize(:b, blockedcrystal.unitcell, 1)
physicalmodes::Subset{Mode} = spanoffset(blockedmodes, crystalpoints)

frozenseedingmodes::Subset{Mode} = circularregionmodes(triangular |> getorigin, physicalmodes, 2.0)
frozenseedingregion::Subset{Offset} = Subset(m |> getpos for m in frozenseedingmodes)
visualize(frozenseedingregion, title="Frozen Seeding Region", visualspace=euclidean(RealSpace, 2))
frozenseedingfock::FockSpace = RegionFock(frozenseedingmodes)

globaldistiller = globaldistillerhamiltonian(
    correlations=blockresult[:correlations],
    regionfock=frozenseedingfock,
    localisometryselectionstrategy=frozenselectionbycount(3),
    symmetry=c6)

globaldistillerspectrum = globaldistiller |> crystalspectrum
visualize(globaldistillerspectrum, title="Global Distiller")

distillresult = groupbands(globaldistillerspectrum, :courier => v -> abs(v) < 1e-5, :filled => v -> v > 1e-5, :empty => v -> v < -1e-5)

courierseedingcenter::Offset = (blockedmodes |> getspace) & [2/3, 1/3]
courierseedingmodes::Subset{Mode} = circularregionmodes(courierseedingcenter, physicalmodes, 1.8)
courierseedingregion::Subset{Offset} = Subset(m |> getpos for m in courierseedingmodes)
visualize(courierseedingregion, courierseedingcenter |> Subset, title="Courier Seeding Region", visualspace=euclidean(RealSpace, 2))
courierseedingfock::RegionFock = RegionFock(courierseedingmodes)

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
purifiedcorrelationspectrum |> visualize
couriercorrelations = purifiedcorrelationspectrum |> FockMap

using ColorTypes
function visualizeregionstate2d(state::RegionState{2}; title::String = "")
    mapdata::SparseMatrixCSC = state |> rep |> rep

    function generatestateplot(spstate::FockMap)
        spmode::Mode = spstate |> getinspace |> first
        columnspectrum::Base.Generator = (m => spstate[m, spmode] for m in spstate |> getoutspace |> orderedmodes)
        positions::Vector{Offset} = [v.first |> pos for v in columnspectrum]
        mesh::Matrix{Float64} = hcat(map(p -> p |> euclidean |> pos, positions)...)
        markersizes::Vector{Float64} = [v.second |> abs for v in columnspectrum]
        normalizedmarkersizes::Vector{Float64} = markersizes / norm(markersizes) * 120
        markercolors::Vector = [convert(RGB{Float64}, HSV(angle(v.second) / 2π * 360, 1, 1)) for v in columnspectrum]
        return scatter(
            x=mesh[1, :], y=mesh[2, :], mode="markers",
            marker=attr(
                symbol="circle",
                size=normalizedmarkersizes,
                color=markercolors))
    end

    scatters::Vector = [spstate |> generatestateplot for spstate in state]
    fig = make_subplots(rows=1, cols=scatters |> length)
    for (n, scatter) in enumerate(scatters)
        add_trace!(fig, scatter, row=1, col=n)
    end
    relayout!(fig, title_text=title)
    fig
end

commutation(c6 * couriercorrelations.outspace, couriercorrelations) |> maximum

blockedfilledprojector = distillresult[:filled] |> crystalprojector
blockedfilledcorrelation = idmap(blockedfilledprojector.outspace, blockedfilledprojector.outspace) - blockedfilledprojector
filledseed = [regionalwannierseeding(blockedfilledcorrelation, frozenseedingfock, symmetry=c6, seedsgroupingprecision=1e-3)...][1]

crystalfilledseeds = crystalisometries(localisometry=filledseed, crystalfock=blockedcorrelations.outspace, addinspacemomentuminfo=true)
wannierfilledisometry = wannierprojection(
    crystalisometries=distillresult[:filled].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalfilledseeds)

state = regionalrestriction(wannierfilledisometry, frozenseedingfock, frozenseedingfock, frozenseedingfock)
state |> visualizeregionstate2d

blockedemptyprojector = distillresult[:empty] |> crystalprojector
blockedemptycorrelation = idmap(blockedemptyprojector.outspace, blockedemptyprojector.outspace) - blockedemptyprojector
emptyseed = [regionalwannierseeding(blockedemptycorrelation, frozenseedingfock, symmetry=c6, seedsgroupingprecision=1e-3)...][1]

crystalemptyseeds = crystalisometries(localisometry=emptyseed, crystalfock=blockedcorrelations.outspace, addinspacemomentuminfo=true)
wannieremptyisometry = wannierprojection(
    crystalisometries=distillresult[:empty].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalemptyseeds)

filledC = wannierfilledisometry' * blockedcorrelations * wannierfilledisometry
filledC |> crystalspectrum |> visualize
emptyC = wannieremptyisometry' * blockedcorrelations * wannieremptyisometry
emptyC |> crystalspectrum |> visualize

filledP = wannierfilledisometry * wannierfilledisometry'
courierP = wanniercourierisometry * wanniercourierisometry'

filledC = idmap(filledP |> getinspace) - filledP
courierC = idmap(courierP |> getinspace) - courierP

evals = regioncorrelations(filledC, frozenseedingfock) |> eigvalsh
plot(scatter(y=[v for (_, v) in evals] |> sort, mode="markers"))

somemodes = spanoffset(blockedhamiltonian |> getoutspace |> unitcellfock |> orderedmodes, blockedcrystal |> latticepoints)

Ft = fourier(filledC.outspace, somemodes |> FockSpace)

realfilledC = Ft' * filledC * Ft

U = regioncorrelations(filledC, frozenseedingfock) |> eigvecsh

M = idmap(somemodes |> FockSpace) - idmap(frozenseedingfock) + FockMap(U, inspace=U |> getoutspace, permute=false)

R = M * realfilledC * M'

visualize(R - realfilledC, rowrange=1:64, colrange=1:64)
R - realfilledC |> maximum
courierR = restrict(R, courierseedingfock, courierseedingfock)
plot(scatter(y=[v for (_, v) in courierR |> eigvalsh] |> sort, mode="markers"))

translatedfrozenmodes = circularregionmodes(triangular & [1, 0], physicalmodes, 2.0)
visualize(frozenseedingregion, Subset(m |> pos for m in translatedfrozenmodes), visualspace=euclidean(RealSpace, 2))

TfrozenR = restrict(R, FockSpace(translatedfrozenmodes), FockSpace(translatedfrozenmodes))
plot(scatter(y=[v for (_, v) in TfrozenR |> eigvalsh] |> sort, mode="markers"))

frozenR = restrict(R, frozenseedingfock, frozenseedingfock)
plot(scatter(y=[v for (_, v) in frozenR |> eigvalsh] |> sort, mode="markers"))

localcourierinfilled = regioncorrelations(courierC, courierseedingfock)
plot(scatter(y=[v for (_, v) in localcourierinfilled |> eigvalsh] |> sort, mode="markers"))
