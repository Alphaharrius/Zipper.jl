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
    gmera1 = @time gmerastep(rgblockedcorrelations,rgblockedcorrelations,firstcenterlist, modeselectionbythreshold(0.001))
    # gmera1 = @time gmerastep1(rgblockedcorrelations,firstcenterlist)
    @info ("1st gmera approximation to correlations...")
    gmera1approx = transistionmap*gmera1[:emptyisometry]*gmera1[:emptyisometry]'*transistionmap'

    if haskey(gmera1,:correlations)
        transistionmap = transistionmap*gmera1[:courierisometry]
        @info ("2nd gmera step...")
        gmera2 = @time gmerastep(rgblockedcorrelations,gmera1[:correlations],secondcenterlist,modeselectionbythreshold(0.001))
        @info ("2nd gmera approximation to correlations...")
        gmera2approx = transistionmap*gmera2[:emptyisometry]*gmera2[:emptyisometry]'*transistionmap'
        if haskey(gmera2,:correlations)
            transistionmap = transistionmap*gmera2[:courierisometry]
            @info ("3rd gmera step...")
            gmera3 = @time gmerastep(rgblockedcorrelations,gmera2[:correlations],thirdcenterlist, modeselectionbythreshold(0.001))
            @info ("3rd gmera approximation to correlations...")
            gmera3approx = transistionmap*gmera3[:emptyisometry]*gmera3[:emptyisometry]'*transistionmap'
            if haskey(gmera3,:correlations)
                transistionmap = transistionmap*gmera3[:courierisometry]
                @info ("final gmera step...")
                gmera4 = @time gmerastep(rgblockedcorrelations,gmera3[:correlations],finalcenterlist,modeselectionbythreshold(0.001))
                @info ("4th gmera approximation to correlations...")
                gmera4approx = transistionmap*gmera4[:emptyisometry]*gmera4[:emptyisometry]'*transistionmap'
                if haskey(gmera4,:correlations)
                    transistionmap = transistionmap*gmera4[:courierisometry]
                    @info("go through all gmera step")
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
                else
                    @info("terminate at 4th gmera step")
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
                        :gmera4thapprox => gmera4approx,
                        :transistionmap => transistionmap)
                end
            else
                @info("terminate at 3rd gmera step")
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
                    :gmera3rdapprox => gmera3approx,
                    :transistionmap => transistionmap)
            end
        else
            @info("terminate at 2nd gmera step")
            return Dict(
                :rgblockedmap => rgblocker,
                :gmera1stemptyisometry => gmera1[:emptyisometry],
                :gmera1stfilledisometry => gmera1[:filledisometry],
                :gmera1stcourierisometry => gmera1[:courierisometry],
                :gmera1stcorrelations => gmera1[:correlations],
                :gmera1stapprox => gmera1approx,
                :gmera2ndemptyisometry => gmera2[:emptyisometry],
                :gmera2ndfilledisometry => gmera2[:filledisometry],
                :gmera2ndapprox => gmera2approx,
                :transistionmap => transistionmap)
        end
    else
        @info("terminate at 1st gmera step")
        return Dict(
        :rgblockedmap => rgblocker,
        :gmera1stemptyisometry => gmera1[:emptyisometry],
        :gmera1stfilledisometry => gmera1[:filledisometry],
        :gmera1stapprox => gmera1approx,
        :transistionmap => transistionmap)
    end
end

rg1 = gmera(blockedcorrelations,idmap(blockedcorrelations|>getinspace))
rg2 = gmera(rg1[:correlations],rg1[:transistionmap])