using Plotly, SmithNormalForm, LinearAlgebra, OrderedCollections, SparseArrays, Combinatorics
using Zipper

setmaxthreads(Threads.nthreads())

function H_haldane(systemsize::Number,tₙ::Number,tₕ::Number)
    triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')
    kspace = convert(MomentumSpace, triangular)
    
    pa = [1/3, 2/3] ∈ triangular
    pb = [2/3, 1/3] ∈ triangular
    pc = (pa + pb) / 2
    spatialsnappingcalibration((pa, pb, pc))
    
    c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
    
    unitcell = Subset(pa, pb)
    crystal = Crystal(unitcell, [systemsize, systemsize])
    reciprocalhashcalibration(crystal.sizes)
    
    modes::Subset{Mode} = quantize(unitcell, 1)|>orderedmodes
    m0, m1 = members(modes)
    
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
    return H, correlations
end

H, correlations = H_haldane(48,ComplexF64(-1),0.1im)

function startingcorrelations(correlations::FockMap, unitcellscale::Number)
    crystalfock = correlations|>getoutspace
    scale = Scale([unitcellscale 0; 0 unitcellscale], crystalfock|>getcrystal|>getspace)
    @info("Performing unitcellblocking...")
    @info("Generating unitcellblocking transformation...")
    blocker = @time scale * crystalfock
    @info("Performing unitcellblocking on correlations...")
    blockedcorrelations = @time blocker * correlations * blocker'
    return blockedcorrelations
end

blockedcorrelations = startingcorrelations(correlations,3)
