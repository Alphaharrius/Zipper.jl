using Plotly, SmithNormalForm, LinearAlgebra, OrderedCollections, SparseArrays, Combinatorics
using Zipper

setmaxthreads(Threads.nthreads())

function focktraceL1norm(fockmap,systemsize)
    modeevalpairs = fockmap|>eigspech|>geteigenvalues
    return(sum([abs(modeevalpair[2]) for modeevalpair in modeevalpairs])/systemsize)
end

function focktraceL2norm(fockmap,systemsize)
    return(real(sqrt(tr(fockmap*fockmap'|>rep))/systemsize))
end

function rgapproximation(rgresults)
    finalrg = Symbol("rg"*string(length(rgresults)))
    if haskey(rgresults[finalrg],:correlations)
        core = @time groupbands(rgresults[finalrg][:correlations]|>crystalspectrum, :filled=> v -> v < 1e-5, :empty => v -> v > 1e-5)
        coreproj = core[:empty]|>crystalprojector
        coreapprox = rgresults[finalrg][:transistionmap]*coreproj*rgresults[finalrg][:transistionmap]'
        return sum([rgresult[2][:approximation] for rgresult in rgresults])+coreapprox
    else
        return sum([rgresult[2][:approximation] for rgresult in rgresults])
    end
end

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


function gmera(correlations,reftransistionmap,selectionstrategy)
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
    transistionmapdict = Dict(:orig=>transistionmap)

    firstcenterlist = [[1/2,0] ∈ rgblockedspace,[-1/2,0] ∈ rgblockedspace]
    secondcenterlist = [[0,1/2] ∈ rgblockedspace,[0,-1/2] ∈ rgblockedspace]
    thirdcenterlist = [[1/2,1/2] ∈ rgblockedspace,[-1/2,-1/2] ∈ rgblockedspace, [-1/2,1/2] ∈ rgblockedspace,[1/2,-1/2] ∈ rgblockedspace]
    finalcenterlist = [[0,0] ∈ rgblockedspace]

    @info ("1st gmera step...")
    gmera1 = @time gmerastep(rgblockedcorrelations,rgblockedcorrelations,firstcenterlist, selectionstrategy)
    # gmera1 = @time gmerastep1(rgblockedcorrelations,firstcenterlist)
    @info ("1st gmera approximation to correlations...")
    gmera1approx = transistionmap*gmera1[:emptyisometry]*gmera1[:emptyisometry]'*transistionmap'
    approximation = gmera1approx
    approxdict = Dict(:first=>gmera1approx)
    filledisometrydict = Dict(:first=>gmera1[:filledisometry])
    emptyisometrydict = Dict(:first=>gmera1[:emptyisometry])

    if haskey(gmera1,:correlations)
        transistionmap = transistionmap*gmera1[:courierisometry]
        correlationsdict = Dict(:first=>gmera1[:correlations])
        transistionmapdict[:first] = transistionmap
        @info ("2nd gmera step...")
        gmera2 = @time gmerastep(rgblockedcorrelations,gmera1[:correlations],secondcenterlist, selectionstrategy)
        @info ("2nd gmera approximation to correlations...")
        gmera2approx = transistionmap*gmera2[:emptyisometry]*gmera2[:emptyisometry]'*transistionmap'
        approximation = approximation+gmera2approx
        approxdict[:second] = gmera2approx
        filledisometrydict[:second] = gmera2[:filledisometry]
        emptyisometrydict[:second] = gmera2[:emptyisometry]
        if haskey(gmera2,:correlations)
            transistionmap = transistionmap*gmera2[:courierisometry]
            transistionmapdict[:second] = transistionmap
            correlationsdict[:second] = gmera2[:correlations]
            @info ("3rd gmera step...")
            gmera3 = @time gmerastep(rgblockedcorrelations,gmera2[:correlations],thirdcenterlist, selectionstrategy)
            @info ("3rd gmera approximation to correlations...")
            gmera3approx = transistionmap*gmera3[:emptyisometry]*gmera3[:emptyisometry]'*transistionmap'
            approximation = approximation+gmera3approx
            approxdict[:third] = gmera3approx
            filledisometrydict[:third] = gmera3[:filledisometry]
            emptyisometrydict[:third] = gmera3[:emptyisometry]
            if haskey(gmera3,:correlations)
                transistionmap = transistionmap*gmera3[:courierisometry]
                transistionmapdict[:third] = transistionmap
                correlationsdict[:third] = gmera3[:correlations]
                @info ("final gmera step...")
                gmera4 = @time gmerastep(rgblockedcorrelations,gmera3[:correlations],finalcenterlist, selectionstrategy)
                @info ("4th gmera approximation to correlations...")
                gmera4approx = transistionmap*gmera4[:emptyisometry]*gmera4[:emptyisometry]'*transistionmap'
                approximation = approximation+gmera4approx
                approxdict[:forth] = gmera4approx
                filledisometrydict[:forth] = gmera4[:filledisometry]
                emptyisometrydict[:forth] = gmera4[:emptyisometry]
                if haskey(gmera4,:correlations)
                    transistionmap = transistionmap*gmera4[:courierisometry]
                    transistionmapdict[:forth] = transistionmap
                    correlationsdict[:forth] = gmera4[:correlations]
                    @info("go through all gmera step")
                    return Dict(
                        :rgblockedmap => rgblocker,
                        :emptyisometryresults => emptyisometrydict,
                        :filledisometryresults => filledisometrydict,
                        :correlationsresults => correlationsdict,
                        :approximation => approximation,
                        :approxresults => approxdict,
                        :transistionmaps => transistionmapdict,
                        :correlations => gmera4[:correlations],
                        :transistionmap => transistionmap)
                else
                    @info("terminate at 4th gmera step")
                    return Dict(
                        :rgblockedmap => rgblocker,
                        :emptyisometryresults => emptyisometrydict,
                        :filledisometryresults => filledisometrydict,
                        :correlationsresults => correlationsdict,
                        :approximation => approximation,
                        :approxresults => approxdict,
                        :transistionmaps => transistionmapdict,
                        :transistionmap => transistionmap)
                end
            else
                @info("terminate at 3rd gmera step")
                return Dict(
                    :rgblockedmap => rgblocker,
                    :emptyisometryresults => emptyisometrydict,
                    :filledisometryresults => filledisometrydict,
                    :correlationsresults => correlationsdict,
                    :approximation => approximation,
                    :approxresults => approxdict,
                    :transistionmaps => transistionmapdict,
                    :transistionmap => transistionmap)
            end
        else
            @info("terminate at 2nd gmera step")
            return Dict(
                :rgblockedmap => rgblocker,
                :emptyisometryresults => emptyisometrydict,
                :filledisometryresults => filledisometrydict,
                :correlationsresults => correlationsdict,
                :approximation => approximation,
                :approxresults => approxdict,
                :transistionmaps => transistionmapdict,
                :transistionmap => transistionmap)
        end
    else
        @info("terminate at 1st gmera step")
        return Dict(
            :rgblockedmap => rgblocker,
            :emptyisometryresults => emptyisometrydict,
            :filledisometryresults => filledisometrydict,
            :aproximation => approximation,
            :approxresults => approxdict,
            :transistionmaps => transistionmapdict,
            :transistionmap => transistionmap)
    end
end

function rg(correlations,selectionstrategy)
    step = 1
    @info("implementing "*string(step)*" rg step")
    gmeraresult = gmera(correlations,idmap(correlations|>getinspace),selectionstrategy)
    rgresult = Dict(Symbol("rg"*string(step))=>gmeraresult)
    step+=1
    while (haskey(gmeraresult,:correlations)) 
        if ((gmeraresult[:correlations]|>getinspace|>getcrystal|>size)[1]==3) || ((gmeraresult[:correlations]|>getinspace|>getcrystal|>size)[1]==2)
            @info("cannot further do RG")
            return rgresult
        else 
            @info("implementing "*string(step)*" rg step and the remaining systemsize "*string(gmeraresult[:correlations]|>getinspace|>getcrystal))
            gmeraresult = gmera(gmeraresult[:correlations],gmeraresult[:transistionmaps][:forth],selectionstrategy)
            rgresult[Symbol("rg"*string(step))] = gmeraresult
            step+=1
        end
    end
    return rgresult
end

modeselectionstrategy = modeselection1stbycountthenbythreshold(1,0.001)

H, correlations = H_haldane(12,ComplexF64(-1),0.4im)

blockedcorrelationsloc2size12chern = startingcorrelations(correlations,2)
rgresultsloc2size12chern = rg(blockedcorrelationsloc2size12chern,modeselectionstrategy)
rgresultsloc2size12chern[:rg1][:correlations]
rgapproxloc2size12chern = rgapproximation(rgresultsloc2size12chern)
dataloc2size12chern = focktraceL2norm(FockMap(blockedcorrelationsloc2size12chern-rgapproxloc2size12chern),288)

blockedcorrelationsloc3size12chern = startingcorrelations(correlations,3)
rgresultsloc3size12chern = rg(blockedcorrelationsloc3size12chern,modeselectionstrategy)
rgresultsloc3size12chern[:rg1][:correlations]
rgapproxloc3size12chern = rgapproximation(rgresultsloc3size12chern)
dataloc3size12chern = focktraceL2norm(FockMap(blockedcorrelationsloc3size12chern-rgapproxloc3size12chern),288)

H, correlations = H_trivial(24,ComplexF64(-1),0.4im)

blockedcorrelationsloc2size24chern = startingcorrelations(correlations,2)
rgresultsloc2size24chern = rg(blockedcorrelationsloc2size24chern,modeselectionstrategy)
rgresultsloc2size24chern[:rg2][:correlations]
rgapproxloc2size24chern = rgapproximation(rgresultsloc2size24chern)
dataloc2size24chern = focktracenorm(FockMap(blockedcorrelationsloc2size24chern-rgapproxloc2size24chern),1152)

blockedcorrelationsloc3size24chern = startingcorrelations(correlations,3)
rgresultsloc3size24chern = rg(blockedcorrelationsloc3size24chern,modeselectionstrategy)
rgresultsloc3size24chern[:rg2][:correlations]
rgapproxloc3size24chern = rgapproximation(rgresultsloc3size24chern)
dataloc3size24chern = focktracenorm(FockMap(blockedcorrelationsloc3size24chern-rgapproxloc3size24chern),1152)

blockedcorrelationsloc4size24chern = startingcorrelations(correlations,4)
rgresultsloc4size24chern = rg(blockedcorrelationsloc4size24chern,modeselectionstrategy)
rgresultsloc4size24chern[:rg1][:correlations]
rgapproxloc4size24chern = rgapproximation(rgresultsloc4size24chern)
dataloc4size24chern = focktracenorm(FockMap(blockedcorrelationsloc4size24chern-rgapproxloc4size24chern),1152)

blockedcorrelationsloc6size24 = startingcorrelations(correlations,6)
rgresultsloc6size24 = rg(blockedcorrelationsloc6size24,modeselectionstrategy)
rgresultsloc6size24[:rg1][:correlations]
rgapproxloc6size24 = rgapproximation(rgresultsloc6size24)
dataloc6size24 = focktracenorm(FockMap(blockedcorrelationsloc6size24-rgapproxloc6size24),1152)

H, correlations = H_trivial(48,ComplexF64(-1),0.4im)

blockedcorrelationsloc2size48 = startingcorrelations(correlations,2)
rgresultsloc2size48 = rg(blockedcorrelationsloc2size48,modeselectionstrategy)
rgresultsloc2size48[:rg3]
rgapproxloc2size48 = rgapproximation(rgresultsloc2size48)
dataloc2size48 = focktracenorm(FockMap(blockedcorrelationsloc2size48-rgapproxloc2size48),4608)

blockedcorrelationsloc3size48 = startingcorrelations(correlations,3)
rgresultsloc3size48 = rg(blockedcorrelationsloc3size48,modeselectionstrategy)
rgresultsloc3size48[:rg3][:correlations]
rgapproxloc3size48 = rgapproximation(rgresultsloc3size48)
dataloc3size48 = focktracenorm(FockMap(blockedcorrelationsloc3size48-rgapproxloc3size48),4608)

blockedcorrelationsloc4size48 = startingcorrelations(correlations,4)
rgresultsloc4size48 = rg(blockedcorrelationsloc4size48,modeselectionstrategy)
rgresultsloc4size48[:rg2][:correlations]
rgapproxloc4size48 = rgapproximation(rgresultsloc4size48)
dataloc4size48 = focktracenorm(FockMap(blockedcorrelationsloc4size48-rgapproxloc4size48),4608)

blockedcorrelationsloc6size48 = startingcorrelations(correlations,6)
rgresultsloc6size48 = rg(blockedcorrelationsloc6size48,modeselectionstrategy)
rgresultsloc6size48[:rg2][:correlations]
rgapproxloc6size48 = rgapproximation(rgresultsloc6size48)
dataloc6size48 = focktracenorm(FockMap(blockedcorrelationsloc6size48-rgapproxloc6size48),4608)

blockedcorrelationsloc8size48 = startingcorrelations(correlations,8)
rgresultsloc8size48 = rg(blockedcorrelationsloc8size48,modeselectionstrategy)
rgresultsloc8size48[:rg1][:correlations]
rgapproxloc8size48 = rgapproximation(rgresultsloc8size48)
dataloc8size48 = focktracenorm(FockMap(blockedcorrelationsloc8size48-rgapproxloc8size48),4608)

blockedcorrelationsloc12size48 = startingcorrelations(correlations,12)
rgresultsloc12size48 = rg(blockedcorrelationsloc12size48,modeselectionstrategy)
rgresultsloc12size48[:rg1][:correlations]
rgapproxloc12size48 = rgapproximation(rgresultsloc12size48)
dataloc12size48 = focktracenorm(FockMap(blockedcorrelationsloc12size48-rgapproxloc12size48),4608)

H, correlations = H_trivial(96,ComplexF64(-1),-0.2,-0.8)

blockedcorrelationsloc2size96 = startingcorrelations(correlations,2)
rgresultsloc2size96 = rg(blockedcorrelationsloc2size96,modeselectionstrategy)
rgresultsloc2size96[:rg3]
rgapproxloc2size96 = rgapproximation(rgresultsloc2size96)
dataloc2size96 = focktracenorm(FockMap(blockedcorrelationsloc2size96-rgapproxloc2size96),4608)

blockedcorrelationsloc3size48 = startingcorrelations(correlations,3)
rgresultsloc3size48 = rg(blockedcorrelationsloc3size48,modeselectionstrategy)
rgresultsloc3size48[:rg3][:correlations]
rgapproxloc3size48 = rgapproximation(rgresultsloc3size48)
dataloc3size48 = focktracenorm(FockMap(blockedcorrelationsloc3size48-rgapproxloc3size48),4608)

blockedcorrelationsloc4size48 = startingcorrelations(correlations,4)
rgresultsloc4size48 = rg(blockedcorrelationsloc4size48,modeselectionstrategy)
rgresultsloc4size48[:rg2][:correlations]
rgapproxloc4size48 = rgapproximation(rgresultsloc4size48)
dataloc4size48 = focktracenorm(FockMap(blockedcorrelationsloc4size48-rgapproxloc4size48),4608)

blockedcorrelationsloc6size48 = startingcorrelations(correlations,6)
rgresultsloc6size48 = rg(blockedcorrelationsloc6size48,modeselectionstrategy)
rgresultsloc6size48[:rg2][:correlations]
rgapproxloc6size48 = rgapproximation(rgresultsloc6size48)
dataloc6size48 = focktracenorm(FockMap(blockedcorrelationsloc6size48-rgapproxloc6size48),4608)

blockedcorrelationsloc8size48 = startingcorrelations(correlations,8)
rgresultsloc8size48 = rg(blockedcorrelationsloc8size48,modeselectionstrategy)
rgresultsloc8size48[:rg1][:correlations]
rgapproxloc8size48 = rgapproximation(rgresultsloc8size48)
dataloc8size48 = focktracenorm(FockMap(blockedcorrelationsloc8size48-rgapproxloc8size48),4608)

blockedcorrelationsloc12size48 = startingcorrelations(correlations,12)
rgresultsloc12size48 = rg(blockedcorrelationsloc12size48,modeselectionstrategy)
rgresultsloc12size48[:rg1][:correlations]
rgapproxloc12size48 = rgapproximation(rgresultsloc12size48)
dataloc12size48 = focktracenorm(FockMap(blockedcorrelationsloc12size48-rgapproxloc12size48),4680)