using Plotly, SmithNormalForm, LinearAlgebra, OrderedCollections, SparseArrays, Combinatorics
using Revise
include("../src/Zipper.jl")
using .Zipper

# function findeigenfunction(transformation::AffineTransform; eigenvalue::Number = 1)::BasisFunction
#     eigenvaluekey::Tuple = hashablecomplex(eigenvalue |> Complex, transformation.eigenvaluehashdenominator)
#     if !haskey(transformation.eigenfunctions, eigenvaluekey)
#         error("No eigenfunction is found for eigenvalue $(eigenvalue)!")
#     end
#     return transformation.eigenfunctions[eigenvaluekey]
# end

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')
kspace = convert(MomentumSpace, triangular)

pa = triangular & [1/3, 2/3]
pb = triangular & [2/3, 1/3]
pc = (pa + pb) / 2
spatialsnappingcalibration((pa, pb, pc))

c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])

unitcell = Subset(pa, pb)
crystal = Crystal(unitcell, [96, 96])
reciprocalhashcalibration(crystal.sizes)

modes::Subset{Mode} = quantize(:pos, unitcell, 1)
m0, m1 = members(modes)

tₙ = -1 + 0im
tₕ = 0.1im

nearestneighbor = [
    (m0, m1) => tₙ,
    (m0, m1 |> setattr(:offset => triangular & [-1, 0])) => tₙ,
    (m0, m1 |> setattr(:offset => triangular & [0, 1])) => tₙ]

haldane = [
    (m0, m0 |> setattr(:offset => triangular & [1, 1])) => tₕ,
    (m0, m0 |> setattr(:offset => triangular & [-1, 0])) => tₕ,
    (m0, m0 |> setattr(:offset => triangular & [0, -1])) => tₕ,
    (m1, m1 |> setattr(:offset => triangular & [1, 1])) => -tₕ,
    (m1, m1 |> setattr(:offset => triangular & [-1, 0])) => -tₕ,
    (m1, m1 |> setattr(:offset => triangular & [0, -1])) => -tₕ]

bonds::FockMap = bondmap([nearestneighbor..., haldane...])

energyspectrum = computeenergyspectrum(bonds, crystal=crystal)
energyspectrum |> visualize

groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=1)
groundstates |> visualize
groundstateprojector = groundstates |> crystalprojector

# correlation = 1 - projector
C = idmap(groundstateprojector.outspace) - groundstateprojector

C |> crystalspectrum |> visualize

function entanglemententropy(correlationspectrum::CrystalSpectrum)
    sumresult = 0
    for (_, v) in correlationspectrum |> geteigenvalues
        if isapprox(v, 0, atol=1e-7) || isapprox(v, 1, atol=1e-7)
            continue
        end
        sumresult += v * log(v) + (1 - v) * log(1 - v)
    end
    return -sumresult
end

entanglemententropy(C |> crystalspectrum) 

# function zipper0(correlations::FockMap)
#     crystalfock = correlations.outspace     # map to the physical space

#     scale = Scale([2 0; 0 2], crystalfock |> getcrystal |> getspace)    
#     blockresult = blocking(:action => scale, :correlations => correlations, :crystal => crystalfock |> getcrystal)      
    
#     blockedcrystal::Crystal = blockresult[:crystal]
#     blockedcorrelations::FockMap = blockresult[:correlations]
#     blocker = blockresult[:transformer]

#     # get modes in a circular region of given radius
#     function circularregionmodes(origin::Offset, physicalmodes::Subset{Mode}, radius::Number)::Subset{Mode}
#         currentspace::RealSpace = correlations.outspace |> getcrystal |> getspace |> orthogonalspace
#         physicalnorm = m -> lineartransform(currentspace, m |> getpos) |> norm
#         return filter(p -> physicalnorm(p - origin) < radius, physicalmodes)
#     end

#     crystalpoints::Subset{Offset} = latticepoints(blockedcrystal)
#     blockedmodes::Subset{Mode} = quantize(:pos, blockedcrystal.unitcell, 1)
#     physicalmodes::Subset{Mode} = spanoffset(blockedmodes, crystalpoints)

#     frozenseedingmodes::Subset{Mode} = circularregionmodes(triangular |> getorigin, physicalmodes, 2.0)         # specify modes in the distill region
#     # frozenseedingregion::Subset{Offset} = Subset(m |> getpos for m in frozenseedingmodes)
#     # visualize(frozenseedingregion, title="Frozen Seeding Region", visualspace=euclidean(RealSpace, 2))
#     frozenseedingfock::FockSpace = FockSpace{Region}(frozenseedingmodes)


#     globaldistiller = globaldistillerhamiltonian(
#         correlations=blockresult[:correlations],
#         restrictspace=frozenseedingfock,
#         localisometryselectionstrategy=frozenselectionbycount(3),   # manually choose the number of filled/empty modes
#         symmetry=c6)

#     globaldistillerspectrum = globaldistiller |> crystalspectrum
#     #visualize(globaldistillerspectrum, title="Global Distiller")

#     # create a Dict of courier/filled/empty modes from h_distill
#     distillresult = distillation(globaldistillerspectrum, :courier => v -> abs(v) < 1e-5, :empty => v -> v > 1e-5, :filled => v -> v < -1e-5)

#     # choose the courier seeding of B and get the corresponding projector  
#     courierseedingcenter::Offset = (blockedmodes |> getspace) & [2/3, 1/3]
#     courierseedingmodes::Subset{Mode} = circularregionmodes(courierseedingcenter, physicalmodes, 1.8)
#     #courierseedingregion::Subset{Offset} = Subset(m |> getpos for m in courierseedingmodes)
#     #visualize(courierseedingregion, courierseedingcenter |> Subset, title="Courier Seeding Region", visualspace=euclidean(RealSpace, 2))
#     courierseedingfock::FockSpace{Region} = FockSpace{Region}(courierseedingmodes)      # create a Fock space out of the courier modes around center of B sublattice

#     c3 = c6^2 |> recenter(courierseedingcenter)
#     # projector constructed from the courier modes in h_distill
#     blockedcourierprojector = distillresult[:courier] |> crystalprojector
#     blockedcouriercorrelation = idmap(blockedcourierprojector.outspace, blockedcourierprojector.outspace) - blockedcourierprojector

#     localcourierseed = findlocalspstates(statecorrelations=blockedcouriercorrelation, regionfock=courierseedingfock, symmetry=c3, spectrumextractpredicate=v -> v < 5e-2)[1]
#     fullcourierseed = localcourierseed + (c6 * localcourierseed.outspace) * localcourierseed * (c6 * localcourierseed.inspace)'

#     crystalcourierseeds = crystalisometries(localisometry=fullcourierseed, crystalfock=blockedcorrelations.outspace, addinspacemomentuminfo=true)

#     wanniercourierisometry = wannierprojection(
#         crystalisometries=distillresult[:courier].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalcourierseeds)

#     # regionalrestriction(wanniercourierisometry, courierseedingfock) |> visualize

#     couriercorrelations = wanniercourierisometry' * blockedcorrelations * wanniercourierisometry
#     # couriercorrelations |> crystalspectrum |> visualize
#     couriercorrelationspectrum = couriercorrelations |> crystalspectrum
#     purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification        # legally a projector
#     # purifiedcorrelationspectrum |> visualize
#     couriercorrelations = purifiedcorrelationspectrum |> FockMap

#     blockedfilledprojector = distillresult[:filled] |> crystalprojector
#     blockedfilledcorrelation = idmap(blockedfilledprojector.outspace, blockedfilledprojector.outspace) - blockedfilledprojector
#     filledseed = findlocalspstates(statecorrelations=blockedfilledcorrelation, regionfock=frozenseedingfock, symmetry=c6, spectrumextractpredicate=v -> v < 1e-2, degeneracythreshold=1e-3)[3]

#     crystalfilledseeds = crystalisometries(localisometry=filledseed, crystalfock=blockedcorrelations.outspace, addinspacemomentuminfo=true)
#     wannierfilledisometry = wannierprojection(
#         crystalisometries=distillresult[:filled].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalfilledseeds)

#     blockedemptyprojector = distillresult[:empty] |> crystalprojector
#     blockedemptycorrelation = idmap(blockedemptyprojector.outspace, blockedemptyprojector.outspace) - blockedemptyprojector
#     emptyseed = findlocalspstates(statecorrelations=blockedemptycorrelation, regionfock=frozenseedingfock, symmetry=c6, spectrumextractpredicate=v -> v < 1e-2, degeneracythreshold=1e-3)[3]

#     crystalemptyseeds = crystalisometries(localisometry=emptyseed, crystalfock=blockedcorrelations.outspace, addinspacemomentuminfo=true)
#     wannieremptyisometry = wannierprojection(
#         crystalisometries=distillresult[:empty].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalemptyseeds)

#     filledcorrelations = wannierfilledisometry' * blockedcorrelations * wannierfilledisometry
#     emptycorrelations = wannieremptyisometry' * blockedcorrelations * wannieremptyisometry
    
#     isometry = wannierfilledisometry + wanniercourierisometry + wannieremptyisometry    # the total isometry (unitary)
#     u_zipper = blocker' * isometry

#     C_total = u_zipper' * correlations * u_zipper # =filledcorrelations + couriercorrelations + emptycorrelations      # total correlation matrix of the resulting lattice

#     return Dict(
#         :blocker => blocker,
#         :blockedcorrelations => blockedcorrelations,
#         :correlations => couriercorrelations, 
#         :courierzipper => wanniercourierisometry' * blocker, 
#         :filledzipper => wannierfilledisometry' * blocker,
#         :emptyzipper => wannieremptyisometry' * blocker,
#         :globaldistiller => globaldistiller, 
#         :filledcorrelations => filledcorrelations, 
#         :emptycorrelations => emptycorrelations,
#         :entanglemententropy => entanglemententropy(couriercorrelationspectrum),
#         :wanniercourierisometry => wanniercourierisometry,
#         :wannierfilledisometry => wannierfilledisometry,
#         :wannieremptyisometry => wannieremptyisometry, 
#         :C_total => C_total,
#         :isometry => isometry,
#         :u_zipper => u_zipper,
#         :frozenseedingfock => frozenseedingfock,
#         :courierseedingfock => courierseedingfock,
#         :distillresult => distillresult)
# end

function zipper(correlations::FockMap)
    crystalfock = correlations.outspace     # map to the physical space

    scale = Scale([2 0; 0 2], crystalfock |> getcrystal |> getspace)    
    blockresult = blocking(:action => scale, :correlations => correlations, :crystal => crystalfock |> getcrystal)      
    
    blockedcrystal::Crystal = blockresult[:crystal]
    blockedcorrelations::FockMap = blockresult[:correlations]
    blocker = blockresult[:transformer]

    # get modes in a circular region of given radius
    function circularregionmodes(origin::Offset, physicalmodes::Subset{Mode}, radius::Number)::Subset{Mode}
        currentspace::RealSpace = correlations.outspace |> getcrystal |> getspace |> orthogonalspace
        physicalnorm = m -> lineartransform(currentspace, m |> getpos) |> norm
        return filter(p -> physicalnorm(p - origin) < radius, physicalmodes)
    end

    crystalpoints::Subset{Offset} = latticepoints(blockedcrystal)
    blockedmodes::Subset{Mode} = quantize(:pos, blockedcrystal.unitcell, 1)
    physicalmodes::Subset{Mode} = spanoffset(blockedmodes, crystalpoints)

    frozenseedingmodes::Subset{Mode} = circularregionmodes(triangular |> getorigin, physicalmodes, 2.0)         # specify modes in the distill region
    # frozenseedingregion::Subset{Offset} = Subset(m |> getpos for m in frozenseedingmodes)
    # visualize(frozenseedingregion, title="Frozen Seeding Region", visualspace=euclidean(RealSpace, 2))
    frozenseedingfock::FockSpace = FockSpace{Region}(frozenseedingmodes)


    globaldistiller = globaldistillerhamiltonian(
        correlations=blockresult[:correlations],
        restrictspace=frozenseedingfock,
        localisometryselectionstrategy=frozenselectionbycount(3),   # manually choose the number of filled/empty modes
        symmetry=c6)

    globaldistillerspectrum = globaldistiller |> crystalspectrum
    #visualize(globaldistillerspectrum, title="Global Distiller")

    # create a Dict of courier/filled/empty modes from h_distill
    distillresult = distillation(globaldistillerspectrum, :courier => v -> abs(v) < 1e-5, :empty => v -> v > 1e-5, :filled => v -> v < -1e-5)

    # choose the courier seeding of B and get the corresponding projector  
    courierseedingcenter::Offset = (blockedmodes |> getspace) & [2/3, 1/3]
    courierseedingmodes::Subset{Mode} = circularregionmodes(courierseedingcenter, physicalmodes, 1.8)
    #courierseedingregion::Subset{Offset} = Subset(m |> getpos for m in courierseedingmodes)
    #visualize(courierseedingregion, courierseedingcenter |> Subset, title="Courier Seeding Region", visualspace=euclidean(RealSpace, 2))
    courierseedingfock::FockSpace{Region} = FockSpace{Region}(courierseedingmodes)      # create a Fock space out of the courier modes around center of B sublattice

    c3 = c6^2 |> recenter(courierseedingcenter)
    # projector constructed from the courier modes in h_distill
    blockedcourierprojector = distillresult[:courier] |> crystalprojector
    blockedcouriercorrelation = idmap(blockedcourierprojector.outspace, blockedcourierprojector.outspace) - blockedcourierprojector

    localcourierseed = findlocalspstates(statecorrelations=blockedcouriercorrelation, regionfock=courierseedingfock, symmetry=c3, spectrumextractpredicate=v -> v < 5e-2)[1]
    fullcourierseed = localcourierseed + (c6 * localcourierseed.outspace) * localcourierseed * (c6 * localcourierseed.inspace)'

    crystalcourierseeds = crystalisometries(localisometry=fullcourierseed, crystalfock=blockedcorrelations.outspace, addinspacemomentuminfo=true)

    wanniercourierisometry = wannierprojection(
        crystalisometries=distillresult[:courier].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalcourierseeds)

    # regionalrestriction(wanniercourierisometry, courierseedingfock) |> visualize

    couriercorrelations = wanniercourierisometry' * blockedcorrelations * wanniercourierisometry
    # couriercorrelations |> crystalspectrum |> visualize
    couriercorrelationspectrum = couriercorrelations |> crystalspectrum
    purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification        # legally a projector
    # purifiedcorrelationspectrum |> visualize
    couriercorrelations = purifiedcorrelationspectrum |> FockMap

    blockedfilledprojector = distillresult[:filled] |> crystalprojector
    blockedfilledcorrelation = idmap(blockedfilledprojector.outspace, blockedfilledprojector.outspace) - blockedfilledprojector
    filledseed = findlocalspstates(statecorrelations=blockedfilledcorrelation, regionfock=frozenseedingfock, symmetry=c6, spectrumextractpredicate=v -> v < 1e-2, degeneracythreshold=1e-3)[3]

    crystalfilledseeds = crystalisometries(localisometry=filledseed, crystalfock=blockedcorrelations.outspace, addinspacemomentuminfo=true)
    wannierfilledisometry = wannierprojection(
        crystalisometries=distillresult[:filled].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalfilledseeds)

    blockedemptyprojector = distillresult[:empty] |> crystalprojector
    blockedemptycorrelation = idmap(blockedemptyprojector.outspace, blockedemptyprojector.outspace) - blockedemptyprojector
    emptyseed = findlocalspstates(statecorrelations=blockedemptycorrelation, regionfock=frozenseedingfock, symmetry=c6, spectrumextractpredicate=v -> v < 1e-2, degeneracythreshold=1e-3)[3]

    crystalemptyseeds = crystalisometries(localisometry=emptyseed, crystalfock=blockedcorrelations.outspace, addinspacemomentuminfo=true)
    wannieremptyisometry = wannierprojection(
        crystalisometries=distillresult[:empty].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalemptyseeds)

    filledcorrelations = wannierfilledisometry' * blockedcorrelations * wannierfilledisometry
    emptycorrelations = wannieremptyisometry' * blockedcorrelations * wannieremptyisometry
    
    isometry = wannierfilledisometry + wanniercourierisometry + wannieremptyisometry    # the total isometry (unitary)
    u_zipper = blocker' * isometry

    C_total = u_zipper' * correlations * u_zipper # =filledcorrelations + couriercorrelations + emptycorrelations      # total correlation matrix of the resulting lattice
    
    #modifiedTotalCorrelationMat = modified_filledcorrelations + correlations + modified_emptycorrelations
    
    swavefockspace = FockSpace(mode for mode in filledcorrelations |> getinspace if swave == (mode |> getorbital(swave)))
    swavefockspace = FockSpace(swavefockspace, reflected=wannierfilledisometry |> getinspace |> getcrystal)
    swavefilledisometry = wannierfilledisometry[:, swavefockspace]
    
    swavefilledprojector = swavefilledisometry * swavefilledisometry'
    swavefilledcorrelations = idmap(swavefilledprojector |> getinspace) - swavefilledprojector
    filledcorrelations_tilde = wannierfilledisometry' * swavefilledcorrelations * wannierfilledisometry
    
    C_tilde = filledcorrelations_tilde + couriercorrelations + emptycorrelations 

    C_previous_tilde = u_zipper * C_tilde * u_zipper'
    C_previous= u_zipper * C_total * u_zipper'

    return Dict(
        :blocker => blocker,
        :blockedcorrelations => blockedcorrelations,
        :correlations => couriercorrelations, 
        :courierzipper => wanniercourierisometry' * blocker, 
        :filledzipper => wannierfilledisometry' * blocker,
        :emptyzipper => wannieremptyisometry' * blocker,
        :globaldistiller => globaldistiller, 
        :filledcorrelations => filledcorrelations, 
        :emptycorrelations => emptycorrelations,
        :entanglemententropy => entanglemententropy(couriercorrelationspectrum),
        :wanniercourierisometry => wanniercourierisometry,
        :wannierfilledisometry => wannierfilledisometry,
        :wannieremptyisometry => wannieremptyisometry, 
        :C_total => C_total,
        :filledcorrelations_tilde => filledcorrelations_tilde,
        :C_tilde => C_tilde,
        :isometry => isometry,
        :u_zipper => u_zipper,
        :C_previous_tilde => C_previous_tilde,
        :C_previous => C_previous,
        :frozenseedingfock => frozenseedingfock,
        :courierseedingfock => courierseedingfock,
        :distillresult => distillresult)
end 

function zipper_even(correlations::FockMap)
    crystalfock = correlations.outspace     # map to the physical space

    scale = Scale([2 0; 0 2], crystalfock |> getcrystal |> getspace)    
    blockresult = blocking(:action => scale, :correlations => correlations, :crystal => crystalfock |> getcrystal)      
    
    blockedcrystal::Crystal = blockresult[:crystal]
    blockedcorrelations::FockMap = blockresult[:correlations]
    blocker = blockresult[:transformer]

    # get modes in a circular region of given radius
    function circularregionmodes(origin::Offset, physicalmodes::Subset{Mode}, radius::Number)::Subset{Mode}
        currentspace::RealSpace = correlations.outspace |> getcrystal |> getspace |> orthogonalspace
        physicalnorm = m -> lineartransform(currentspace, m |> getpos) |> norm
        return filter(p -> physicalnorm(p - origin) < radius, physicalmodes)
    end

    crystalpoints::Subset{Offset} = latticepoints(blockedcrystal)
    blockedmodes::Subset{Mode} = quantize(:pos, blockedcrystal.unitcell, 1)
    physicalmodes::Subset{Mode} = spanoffset(blockedmodes, crystalpoints)

    frozenseedingmodes::Subset{Mode} = circularregionmodes(triangular |> getorigin, physicalmodes, 2.0)         # specify modes in the distill region
    # frozenseedingregion::Subset{Offset} = Subset(m |> getpos for m in frozenseedingmodes)
    # visualize(frozenseedingregion, title="Frozen Seeding Region", visualspace=euclidean(RealSpace, 2))
    frozenseedingfock::FockSpace = FockSpace{Region}(frozenseedingmodes)


    globaldistiller = globaldistillerhamiltonian(
        correlations=blockresult[:correlations],
        restrictspace=frozenseedingfock,
        localisometryselectionstrategy=frozenselectionbycount(3),   # manually choose the number of filled/empty modes
        symmetry=c6)

    globaldistillerspectrum = globaldistiller |> crystalspectrum
    #visualize(globaldistillerspectrum, title="Global Distiller")

    # create a Dict of courier/filled/empty modes from h_distill
    distillresult = distillation(globaldistillerspectrum, :courier => v -> abs(v) < 1e-5, :empty => v -> v > 1e-5, :filled => v -> v < -1e-5)

    # choose the courier seeding of B and get the corresponding projector  
    courierseedingcenter::Offset = (blockedmodes |> getspace) & [2/3, 1/3]
    courierseedingmodes::Subset{Mode} = circularregionmodes(courierseedingcenter, physicalmodes, 1.8)
    #courierseedingregion::Subset{Offset} = Subset(m |> getpos for m in courierseedingmodes)
    #visualize(courierseedingregion, courierseedingcenter |> Subset, title="Courier Seeding Region", visualspace=euclidean(RealSpace, 2))
    courierseedingfock::FockSpace{Region} = FockSpace{Region}(courierseedingmodes)      # create a Fock space out of the courier modes around center of B sublattice

    c3 = c6^2 |> recenter(courierseedingcenter)
    # projector constructed from the courier modes in h_distill
    blockedcourierprojector = distillresult[:courier] |> crystalprojector
    blockedcouriercorrelation = idmap(blockedcourierprojector.outspace, blockedcourierprojector.outspace) - blockedcourierprojector

    localcourierseed = findlocalspstates(statecorrelations=blockedcouriercorrelation, regionfock=courierseedingfock, symmetry=c3, spectrumextractpredicate=v -> v < 5e-2)[1]
    fullcourierseed = localcourierseed + (c6 * localcourierseed.outspace) * localcourierseed * (c6 * localcourierseed.inspace)'

    crystalcourierseeds = crystalisometries(localisometry=fullcourierseed, crystalfock=blockedcorrelations.outspace, addinspacemomentuminfo=true)

    wanniercourierisometry = wannierprojection(
        crystalisometries=distillresult[:courier].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalcourierseeds)

    # regionalrestriction(wanniercourierisometry, courierseedingfock) |> visualize

    couriercorrelations = wanniercourierisometry' * blockedcorrelations * wanniercourierisometry
    # couriercorrelations |> crystalspectrum |> visualize
    couriercorrelationspectrum = couriercorrelations |> crystalspectrum
    purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification        # legally a projector
    # purifiedcorrelationspectrum |> visualize
    couriercorrelations = purifiedcorrelationspectrum |> FockMap

    blockedfilledprojector = distillresult[:filled] |> crystalprojector
    blockedfilledcorrelation = idmap(blockedfilledprojector.outspace, blockedfilledprojector.outspace) - blockedfilledprojector
    filledseed = findlocalspstates(statecorrelations=blockedfilledcorrelation, regionfock=frozenseedingfock, symmetry=c6, spectrumextractpredicate=v -> v < 1e-2, degeneracythreshold=1e-3)[3]

    crystalfilledseeds = crystalisometries(localisometry=filledseed, crystalfock=blockedcorrelations.outspace, addinspacemomentuminfo=true)
    wannierfilledisometry = wannierprojection(
        crystalisometries=distillresult[:filled].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalfilledseeds)

    blockedemptyprojector = distillresult[:empty] |> crystalprojector
    blockedemptycorrelation = idmap(blockedemptyprojector.outspace, blockedemptyprojector.outspace) - blockedemptyprojector
    emptyseed = findlocalspstates(statecorrelations=blockedemptycorrelation, regionfock=frozenseedingfock, symmetry=c6, spectrumextractpredicate=v -> v < 1e-2, degeneracythreshold=1e-3)[3]

    crystalemptyseeds = crystalisometries(localisometry=emptyseed, crystalfock=blockedcorrelations.outspace, addinspacemomentuminfo=true)
    wannieremptyisometry = wannierprojection(
        crystalisometries=distillresult[:empty].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalemptyseeds)

    filledcorrelations = wannierfilledisometry' * blockedcorrelations * wannierfilledisometry
    emptycorrelations = wannieremptyisometry' * blockedcorrelations * wannieremptyisometry
    
    isometry = wannierfilledisometry + wanniercourierisometry + wannieremptyisometry    # the total isometry (unitary)
    u_zipper = blocker' * isometry

    C_total = u_zipper' * correlations * u_zipper # =filledcorrelations + couriercorrelations + emptycorrelations      # total correlation matrix of the resulting lattice
    
    #modifiedTotalCorrelationMat = modified_filledcorrelations + correlations + modified_emptycorrelations
    
    swavefockspace = FockSpace(mode for mode in filledcorrelations |> getinspace if (mode|> getorbital(swave) |> rep) == [-1+0im])
    swavefockspace = FockSpace(swavefockspace, reflected=wannierfilledisometry |> getinspace |> getcrystal)
    swavefilledisometry = wannierfilledisometry[:, swavefockspace]
    
    swavefilledprojector = swavefilledisometry * swavefilledisometry'
    swavefilledcorrelations = idmap(swavefilledprojector |> getinspace) - swavefilledprojector
    filledcorrelations_tilde = wannierfilledisometry' * swavefilledcorrelations * wannierfilledisometry
    
    C_tilde = filledcorrelations_tilde + couriercorrelations + emptycorrelations 

    C_previous_tilde = u_zipper * C_tilde * u_zipper'
    C_previous= u_zipper * C_total * u_zipper'

    return Dict(
        :blocker => blocker,
        :blockedcorrelations => blockedcorrelations,
        :correlations => couriercorrelations, 
        :courierzipper => wanniercourierisometry' * blocker, 
        :filledzipper => wannierfilledisometry' * blocker,
        :emptyzipper => wannieremptyisometry' * blocker,
        :globaldistiller => globaldistiller, 
        :filledcorrelations => filledcorrelations, 
        :emptycorrelations => emptycorrelations,
        :entanglemententropy => entanglemententropy(couriercorrelationspectrum),
        :wanniercourierisometry => wanniercourierisometry,
        :wannierfilledisometry => wannierfilledisometry,
        :wannieremptyisometry => wannieremptyisometry, 
        :C_total => C_total,
        :filledcorrelations_tilde => filledcorrelations_tilde,
        :C_tilde => C_tilde,
        :isometry => isometry,
        :u_zipper => u_zipper,
        :C_previous_tilde => C_previous_tilde,
        :C_previous => C_previous,
        :frozenseedingfock => frozenseedingfock,
        :courierseedingfock => courierseedingfock,
        :distillresult => distillresult)
end 

# rg_tree = []
# function zer(correlations::FockMap)
#     i=0
#     push!(rg_tree, zipper(correlations))
#     try push!(rg_tree, zipper(correlations))
#         correlations = zipper(correlations)[:correlations]
#         i += 1
#     catch
# end

@time begin
rg1 = zipper(C)  
end 
@time begin
rg2 = zipper_even(rg1[:correlations]) 
end 
@time begin
rg3 = zipper(rg2[:correlations]) 
end 
@time begin
rg4 = zipper_even(rg3[:correlations])  
end

#####25/11/23###########################################################################################################################
test_corr = rg4[:filledzipper]' * rg4[:filledcorrelations] * rg4[:filledzipper] + rg4[:courierzipper]' * rg4[:correlations] * rg4[:courierzipper] + rg4[:emptyzipper]' * rg4[:emptycorrelations] * rg4[:emptyzipper]
rg4[:C_previous] - test_corr |> visualize # valid or not???
# anyway, try it
C3_tilde = rg4[:filledzipper]' * rg4[:filledcorrelations_tilde] * rg4[:filledzipper] + rg4[:emptyzipper]' * rg4[:emptycorrelations] * rg4[:emptyzipper] + rg4[:courierzipper]' *  rg4[:correlations] * rg4[:courierzipper]
C2_tilde = rg3[:filledzipper]' * rg3[:filledcorrelations_tilde] * rg3[:filledzipper] + rg3[:emptyzipper]' * rg3[:emptycorrelations] * rg3[:emptyzipper] + rg3[:courierzipper]' *  C3_tilde * rg3[:courierzipper]
C1_tilde = rg2[:filledzipper]' * rg2[:filledcorrelations_tilde] * rg2[:filledzipper] + rg2[:emptyzipper]' * rg2[:emptycorrelations] * rg2[:emptyzipper] + rg2[:courierzipper]' *  C2_tilde * rg2[:courierzipper]
C_tilde = rg1[:filledzipper]' * rg1[:filledcorrelations_tilde] * rg1[:filledzipper] + rg1[:emptyzipper]' * rg1[:emptycorrelations] * rg1[:emptyzipper] + rg1[:courierzipper]' *  C1_tilde * rg1[:courierzipper]  




#####23/11/23###########################################################################################################################
# function UV_Correlations(IR_correlation::FockMap)::FockMap
#     fff
# end


# wrong, need to add filled and empty terms
C3_tilde = rg4[:C_previous_tilde]
C2_tilde = rg3[:courierzipper]' *  C3_tilde * rg3[:courierzipper]
C1_tilde = rg2[:courierzipper]' *  C2_tilde * rg2[:courierzipper]
C_tilde = rg1[:courierzipper]' *  C1_tilde * rg1[:courierzipper]  
C
# mat = C_tilde |> rep |> Matrix |> Hermitian  
# plot(scatter(y=mat |> eigvals)) 

C_tilde |> visualize 
visualize(C_tilde, colrange=1:1000, rowrange=1:1000)

#####22/11/23###########################################################################################################################
"""
rg1filled = rg1[:filledcorrelations] |> rep
rg2filled = rg2[:filledcorrelations] |> rep

filled_Fock1 = rg1[:filledcorrelations] |> getinspace 
list1 = [] 
for mode in filled_Fock1
    push!(list1, mode |> getorbital(swave))
    # getorbital(swave)
    # if mode|> getorbital(swave) == swave
    #     push!(list, mode)
    # end
end 

filled_Fock2 = rg2[:filledcorrelations] |> getinspace 
list2 = [] 
for mode in filled_Fock2
    push!(list2, mode |> getorbital(swave))
    #push!(list2, mode |> getorbital(swave))
    # getorbital(swave)
    # if mode|> getorbital(swave) == swave
    #     push!(list, mode)
    # end
end 

list2p = [] 
for mode in filled_Fock2
    if mode|> getorbital(swave) |> rep == [-1+0im]
        push!(list2p, mode)
    end
    #push!(list2, mode |> getorbital(swave))
    # getorbital(swave)
    # if mode|> getorbital(swave) == swave
    #     push!(list, mode)
    # end
end 

swavefockspace = FockSpace(mode for mode in rg2[:filledcorrelations] |> getinspace if (mode|> getorbital(swave) |> rep) == [-1+0im])
swavefockspace = FockSpace(swavefockspace, reflected=rg2[:wannierfilledisometry] |> getinspace |> getcrystal)
swavefilledisometry = rg2[:wannierfilledisometry][:, swavefockspace]

swavefilledprojector = swavefilledisometry * swavefilledisometry'
swavefilledcorrelations = idmap(swavefilledprojector |> getinspace) - swavefilledprojector
filledcorrelations_tilde = rg2[:wannierfilledisometry]' * swavefilledcorrelations * rg2[:wannierfilledisometry]
filledcorrelations_tilde |> visualize

C_tilde = filledcorrelations_tilde + rg2[:correlations] + rg2[:emptycorrelations] 
C_tilde |> visualize
C_previous_tilde = rg2[:u_zipper] * C_tilde * rg2[:u_zipper]'  
C_previous_tilde |> visualize
mat = C_previous_tilde |> rep |> Matrix |> Hermitian
plot(scatter(y=mat |> eigvals))

filled_Fock3 = rg3[:filledcorrelations] |> getinspace 
list3 = [] 
for mode in filled_Fock3
    if mode|> getorbital(swave) == swaveminus
        push!(list, mode)
    end
    #push!(list3, mode |> getorbital(swave))
    # getorbital(swave)
    # if mode|> getorbital(swave) == swave
    #     push!(list, mode)
    # end
end 

filled_Fock4 = rg4[:filledcorrelations] |> getinspace 
list4 = [] 
for mode in filled_Fock4
    push!(list4, mode |> getorbital(swave))
    # getorbital(swave)
    # if mode|> getorbital(swave) == swave
    #     push!(list, mode)
    # end
end 

mode for mode in filledcorrelations |> getinspace if swave == (mode |> getorbital(swave))

rg1[:C_previous] |> visualize

crystal |> sitepoints |> visualize

rg3[:correlations] |> getoutspace |> getcrystal |> sitepoints |> visualize

# rg1[:filledzipper]' * rg1[:filledcorrelations_tilde] * rg1[:filledzipper] |> visualize
# rg1[:emptyzipper]' * rg1[:emptycorrelations] * rg1[:emptyzipper] |> visualize
# rg1[:filledzipper]' * rg1[:filledcorrelations] * rg1[:filledzipper]+ rg1[:emptyzipper]' * rg1[:emptycorrelations] * rg1[:emptyzipper] + rg1[:courierzipper]' * rg1[:correlations] * rg1[:courierzipper] |> visualize
# rg1[:courierzipper]' * rg1[:correlations] * rg1[:courierzipper] |> visualize

# corr = rg1[:filledcorrelations_tilde] + rg1[:correlations] + rg1[:emptycorrelations]
# rg1[:u_zipper]
# physcorr |> visualize
# physcorr = rg1[:u_zipper] * corr * rg1[:u_zipper]'
# PhyscorrMat = physcorr |> rep |> Matrix
# import Pkg
# using Plots
# layout = Layout(yaxis=attr(autorange="reversed"))
# plot(heatmap(z=real(PhyscorrMat)), layout)

plot(scatter(y=physcorr |> rep |> Matrix |> Hermitian |> eigvals))

rg1[:C_previous_tilde] |> visualize

mat = rg1[:C_previous_tilde] |> rep |> Matrix |> Hermitian

###########################################################################################################################

#####22/11/23#####################################################################################
crystalfock = C.outspace     # map to the physical space

scale = Scale([2 0; 0 2], crystalfock |> getcrystal |> getspace)    
blockresult = blocking(:action => scale, :correlations => C, :crystal => crystalfock |> getcrystal)      
    
blockedcrystal::Crystal = blockresult[:crystal]

crystalpoints::Subset{Offset} = latticepoints(blockedcrystal)
blockedmodes::Subset{Mode} = quantize(:pos, blockedcrystal.unitcell, 1)
physicalmodes::Subset{Mode} = spanoffset(blockedmodes, crystalpoints)

function circularregionmodes(origin::Offset, physicalmodes::Subset{Mode}, radius::Number)::Subset{Mode}
    currentspace::RealSpace = C.outspace |> getcrystal |> getspace |> orthogonalspace
    physicalnorm = m -> lineartransform(currentspace, m |> getpos) |> norm
    return filter(p -> physicalnorm(p - origin) < radius, physicalmodes)
end

frozenseedingmodes::Subset{Mode} = circularregionmodes(triangular |> getorigin, physicalmodes, 4)         # specify modes in the distill region
    # frozenseedingregion::Subset{Offset} = Subset(m |> getpos for m in frozenseedingmodes)
    # visualize(frozenseedingregion, title="Frozen Seeding Region", visualspace=euclidean(RealSpace, 2))
frozenseedingfock::FockSpace = FockSpace{Region}(frozenseedingmodes)

# filled2physical = rg1[:courierzipper]' * rg2[:filledzipper]' * rg2[:filledcorrelations] * rg2[:filledzipper] * rg1[:courierzipper]

function visualize(state::RegionState{2}; title::String = "", markersizemultiplier::Real = 120)
    function generatestateplot(spstate::FockMap)
        spmode::Mode = spstate |> getinspace |> first
        columnspectrum::Base.Generator = (m => spstate[m, spmode] for m in spstate |> getoutspace |> orderedmodes)
        positions::Vector{Offset} = [v.first |> getpos for v in columnspectrum]
        mesh::Matrix{Float64} = hcat(map(p -> p |> euclidean |> vec, positions)...)
        markersizes::Vector{Float64} = [v.second |> abs for v in columnspectrum]
        normalizedmarkersizes::Vector{Float64} = markersizes / norm(markersizes) * markersizemultiplier
        markercolors::Vector = [convert(RGB{Float64}, HSV(angle(v.second) / 2π * 360, 1, 1)) for v in columnspectrum]
        return scatter(
            x=mesh[1, :], y=mesh[2, :], mode="markers",
            marker=attr(
                symbol="circle",
                size=normalizedmarkersizes,
                color=markercolors))
    end

    scatters::Vector = [spstate |> generatestateplot for (_, spstate) in state]
    subplotnames::Base.Generator = (mode |> string for mode in state |> getinspace |> orderedmodes |> indexmodes)
    fig = make_subplots(rows=1, cols=scatters |> length, subplot_titles=hcat(subplotnames...))
    for (n, scatter) in enumerate(scatters)
        add_trace!(fig, scatter, row=1, col=n)
    end
    relayout!(fig, title_text=title)
    fig
end

regionalrestriction(rg1[:wannierfilledisometry], frozenseedingfock) |> visualize

rg2filled_physical_isometry = rg1[:wanniercourierisometry]*rg2[:blocker]'*rg2[:wannierfilledisometry]
regionalrestriction(rg2filled_physical_isometry, frozenseedingfock) |> visualize

rg3filled_physical_isometry = rg1[:wanniercourierisometry]*rg2[:blocker]'*rg2[:wanniercourierisometry]*rg3[:blocker]'*rg3[:wannierfilledisometry]
regionalrestriction(rg3filled_physical_isometry, frozenseedingfock) |> visualize

##################################################################################################################################



c6.transformmatrix |> eigvals

pplus = findeigenfunction(c6, eigenvalue=(exp(im * π/3)))

swavefilledfockspace = FockSpace(mode for mode in wannierfilledunitary |> getinspace if (mode |> getorbital(swave)) == swave)

correlations[swavefilledfockspace, swavefilledfockspace]




H = energyspectrum |> FockMap
energyspectrum |> visualize

rg1[:blocker] * H * rg1[:blocker]' |> crystalspectrum |> visualize

rg1[:courierzipper] * H * rg1[:courierzipper]' |> crystalspectrum |> visualize
rg1[:filledzipper] * H * rg1[:filledzipper]' |> crystalspectrum |> visualize
rg1[:emptyzipper] * H * rg1[:emptyzipper]' |> crystalspectrum |> visualize

courierH = rg1[:courierzipper] * H * rg1[:courierzipper]'
rg2[:blocker] * courierH * rg2[:blocker]' |> crystalspectrum |> visualize

rg2[:courierzipper] * rg1[:courierzipper] * H * rg1[:courierzipper]' * rg2[:courierzipper]' |> crystalspectrum |> visualize
rg2[:filledzipper] * rg1[:courierzipper] * H * rg1[:courierzipper]' * rg2[:filledzipper]' |> crystalspectrum |> visualize
rg2[:emptyzipper] * rg1[:courierzipper] * H * rg1[:courierzipper]' * rg2[:emptyzipper]' |> crystalspectrum |> visualize

courierH = rg2[:courierzipper] * courierH * rg2[:courierzipper]'
rg3[:blocker] * courierH * rg3[:blocker]' |> crystalspectrum |> visualize

rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper] * H * rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' |> crystalspectrum |> visualize
rg3[:filledzipper] * rg2[:courierzipper] * rg1[:courierzipper] * H * rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:filledzipper]' |> crystalspectrum |> visualize
rg3[:emptyzipper] * rg2[:courierzipper] * rg1[:courierzipper] * H * rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:emptyzipper]' |> crystalspectrum |> visualize

courierH = rg3[:courierzipper] * courierH * rg3[:courierzipper]'
rg4[:blocker] * courierH * rg4[:blocker]' |> crystalspectrum |> visualize

rg4[:courierzipper] * rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper] * H * rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' * rg4[:courierzipper]' |> crystalspectrum |> visualize
rg4[:filledzipper] * rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper] * H * rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' * rg4[:filledzipper]' |> crystalspectrum |> visualize
rg4[:emptyzipper] * rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper] * H * rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' * rg4[:emptyzipper]' |> crystalspectrum |> visualize

function momentumoccupations(correlations::FockMap)::FockMap
    kcorrelations::Base.Generator = correlations |> crystalsubmaps
    crystal::Crystal = correlations |> getoutspace |> getcrystal
    center::Offset = crystal |> getspace |> getorigin
    tracecrystal::Crystal = Crystal(center |> Subset, crystal |> size)
    mode::Mode = Mode(:pos => center)

    function tracing(k::Momentum, corr::FockMap)::FockMap
        space::FockSpace = mode |> setattr(:offset => k) |> FockSpace
        return FockMap(space, space, [corr |> tr][:, :] |> SparseMatrixCSC) / dimension(corr |> getoutspace)
    end

    occupations::FockMap = directsum(tracing(k, corr) for (k, corr) in kcorrelations)
    fockspace::FockSpace = FockSpace(occupations |> getoutspace, reflected=tracecrystal)
    return FockMap(occupations, inspace=fockspace, outspace=fockspace)
end

using Compat
function Zipper.visualize(spectrum::CrystalSpectrum{2}; title="", toppadding::Bool = true)
    kspectrum::Dict{Momentum} = Dict(k => ([spectrum.eigenvalues[m] for m in modes] |> sort) for (k, modes) in spectrum.eigenmodes)
    mesh::Matrix{Momentum} = spectrum.crystal |> brillouinmesh
    plottingdata::Matrix{Vector} = map(k -> haskey(kspectrum, k) ? kspectrum[k] : [], mesh)
    bandcount::Integer = map(v -> v |> length, plottingdata) |> maximum
    function padding(v::Vector)::Vector
        return toppadding ? vcat(v, repeat([NaN], bandcount - length(v))) : vcat(repeat([NaN], bandcount - length(v)), v)
    end
    paddeddata::Matrix{Vector} = map(v -> v |> padding, plottingdata)
    plottingspectrum::Array = paddeddata |> stack
    layout::Layout = Layout(title=title)
    plot([contour(z=plottingspectrum[n, :, :]) for n in axes(plottingspectrum, 1)], layout)
end

filled = rg1[:filledzipper]' * rg1[:filledzipper]
empty = rg1[:emptyzipper]' * rg1[:emptyzipper]

filled1 = rg1[:courierzipper]' * rg2[:filledzipper]' * rg2[:filledzipper] * rg1[:courierzipper]
empty1 = rg1[:courierzipper]' * rg2[:emptyzipper]' * rg2[:emptyzipper] * rg1[:courierzipper]

filled2 = rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:filledzipper]' * rg3[:filledzipper] * rg2[:courierzipper] * rg1[:courierzipper]
empty2 = rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:emptyzipper]' * rg3[:emptyzipper] * rg2[:courierzipper] * rg1[:courierzipper]

filled3 = rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' * rg4[:filledzipper]' * rg4[:filledzipper] * rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper]
empty3 = rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' * rg4[:emptyzipper]' * rg4[:emptyzipper] * rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper]

core = rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' * rg4[:courierzipper]' * rg4[:courierzipper] * rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper]

occ1 = momentumoccupations(filled) + momentumoccupations(empty)
occ2 = momentumoccupations(filled1) + momentumoccupations(empty1)
occ3 = momentumoccupations(filled2) + momentumoccupations(empty2)
occ4 = momentumoccupations(filled3) + momentumoccupations(empty3)
occcore = momentumoccupations(core)

occ1 |> crystalspectrum |> visualize
occ2 |> crystalspectrum |> visualize
occ3 |> crystalspectrum |> visualize
occ4 |> crystalspectrum |> visualize
occcore |> crystalspectrum |> visualize

occ1 + occ2 + occ3 + occ4 + occcore |> crystalspectrum |> visualize

occ1 + occ2 + occ3 + occcore |> crystalspectrum |> visualize

frozenocc = momentumoccupations(rg1[:filledcorrelations] + rg1[:emptycorrelations]) |> crystalspectrum
frozenocc |> FockMap |> eigspech |> visualize

occ0 = momentumoccupations(C)
occ1 = momentumoccupations(rg1[:courierzipper]' * rg1[:correlations] * rg1[:courierzipper])
occ2 = momentumoccupations(rg1[:courierzipper]' * rg2[:courierzipper]' * rg2[:correlations] * rg2[:courierzipper] * rg1[:courierzipper]) - occ1
occ3 = momentumoccupations(rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' * rg3[:correlations] * rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper]) - occ2 - occ1

occ0 |> crystalspectrum |> visualize
occ1 |> crystalspectrum |> visualize
occ2 |> crystalspectrum |> visualize
occ3 |> crystalspectrum |> visualize

occ1 + occ2 + occ3 |> crystalspectrum |> visualize

"""
