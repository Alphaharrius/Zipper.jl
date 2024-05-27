using Pkg
Pkg.resolve()
Pkg.instantiate()
Pkg.activate("../../../../../../../")
Pkg.resolve()
Pkg.instantiate()
using Zipper
using LinearAlgebra

setmaxthreads(8)

function focktraceL1norm(fockmap,volume)
    @info("Calculating L1norm...")
    modeevalpairs = fockmap|>eigspech|>geteigenvalues
    return(sum([abs(modeevalpair[2]) for modeevalpair in modeevalpairs])/volume)
end

function focktraceL2norm(fockmap,volume)
    @info("Calculating L2norm...")
    return (real(sqrt(tr(fockmap*fockmap'|>rep)/volume)))
end

function rgapproximation(rgresults)
    finalrg = Symbol("rg"*string(length(rgresults)))
    if haskey(rgresults[finalrg],:correlations)
        core = @time groupbands(rgresults[finalrg][:correlations]|>crystalspectrum, :filled=> v -> v < 1e-5, :empty => v -> v > 1e-5)
        coreproj = core[:empty]|>crystalprojector
        coreapprox = rgresults[finalrg][:transistionmap]*coreproj*rgresults[finalrg][:transistionmap]'
        return sum([rgresult[2][:approximation] for rgresult in rgresults])+(coreapprox|>CrystalDenseMap)
    else
        return sum([rgresult[2][:approximation] for rgresult in rgresults])
    end
end

function saveisometry(rgresults)
    for l in range(1,length(rgresults))
        rgstep = Symbol("rg"*string(l))
        for (key,value) in rgresults[rgstep][:emptyisometryresults]
            fiosave(value,name="rg"*string(l)*string(key)*"emptyisometryresults")
        end
        for (key,value) in rgresults[rgstep][:filledisometryresults]
            fiosave(value,name="rg"*string(l)*string(key)*"filledisometryresults")
        end
        for (key,value) in rgresults[rgstep][:correlationsresults]
            fiosave(value,name="rg"*string(l)*string(key)*"correlationsresults")
        end
        for (key,value) in rgresults[rgstep][:transistionmaps]
            fiosave(value,name="rg"*string(l)*string(key)*"transistionmaps")
        end
    end
end



function countnoofcoremodes(rgresult)
    if haskey(rgresult[Symbol("rg"*string(length(rgresult)))],:correlations)
        noofcoremodes = rgresult[Symbol("rg"*string(length(rgresult)))][:correlations]|>getoutspace|>getmodes|>length
    else
        noofcoremodes = 0
    end
    return noofcoremodes
end

function H_trivial(systemsize::Number,tₙ::Number,t_a::Number,t_b::Number)
    triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')

    pa = [1/3, 2/3] ∈ triangular
    pb = [2/3, 1/3] ∈ triangular
    pc = (pa + pb) / 2
    spatialsnappingcalibration((pa, pb, pc))

    # c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])

    unitcell = Subset(pa, pb)
    crystal = Crystal(unitcell, [systemsize, systemsize])
    reciprocalhashcalibration(crystal.sizes)

    modes::Subset{Mode} = quantize(unitcell, 1)|>orderedmodes
    m0, m1 = members(modes)


    bonds::FockMap = bondmap([
        (m0, m0) => t_a,
        (m1, m1) => t_b,
        (m0, m1) => tₙ,
        (m0, setattr(m1, :r => Point([-1, 0], triangular))) => tₙ,
        (m0, setattr(m1, :r => Point([0, 1], triangular))) => tₙ])

    energyspectrum = @time computeenergyspectrum(bonds, crystal=crystal)

    groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=1)

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

function gmera(correlations,reftransistionmap,selectionstrategy,terminate)
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
    transistionmapdict = Dict(:orig=>(transistionmap))

    firstcenterlist = [[0,0] ∈ rgblockedspace]
    secondcenterlist = [[1/2,0] ∈ rgblockedspace,[-1/2,0] ∈ rgblockedspace]
    thirdcenterlist = [[0,1/2] ∈ rgblockedspace,[0,-1/2] ∈ rgblockedspace]
    finalcenterlist = [[1/2,1/2] ∈ rgblockedspace,[-1/2,-1/2] ∈ rgblockedspace, [-1/2,1/2] ∈ rgblockedspace,[1/2,-1/2] ∈ rgblockedspace]

    if terminate
        selectionstrategy1 = modeselectionbycount(224)
        selectionstrategy2 = modeselectionbycount(16)
        selectionstrategy3 = modeselectionbycount(14)
        selectionstrategy4 = modeselectionbycount(2)
        # selectionstrategy1 = modeselectionbycount(224)
        # selectionstrategy2 = modeselectionbycount(16)
        # selectionstrategy3 = modeselectionbycount(8)
        # selectionstrategy4 = modeselectionbycount(8)
    else
        selectionstrategy1 = modeselectionbycount(64)
        selectionstrategy2 = modeselectionbycount(64)
        selectionstrategy3 = modeselectionbycount(32)
        selectionstrategy4 = modeselectionbycount(32)
    end

    @info ("1st gmera step...")
    gmera1 = @time gmerastep(rgblockedcorrelations,rgblockedcorrelations,firstcenterlist, selectionstrategy1)
    # gmera1 = @time gmerastep1(rgblockedcorrelations,firstcenterlist)
    @info ("1st gmera approximation to correlations...")
    gmera1approx = transistionmap*gmera1[:emptyisometry]*gmera1[:emptyisometry]'*transistionmap'
    approximation = gmera1approx|>CrystalDenseMap
    approxdict = Dict(:first=>gmera1approx)
    filledisometrydict = Dict(:first=>(gmera1[:filledisometry])|>CrystalDenseMap)
    emptyisometrydict = Dict(:first=>(gmera1[:emptyisometry])|>CrystalDenseMap)

    if haskey(gmera1,:correlations)
        transistionmap = transistionmap*gmera1[:courierisometry]
        correlationsdict = Dict(:first=>(gmera1[:correlations])|>CrystalDenseMap)
        transistionmapdict[:first] = transistionmap
        @info ("2nd gmera step...")
        gmera2 = @time gmerastep(rgblockedcorrelations,gmera1[:correlations],secondcenterlist, selectionstrategy2)
        @info ("2nd gmera approximation to correlations...")
        gmera2approx = transistionmap*gmera2[:emptyisometry]*gmera2[:emptyisometry]'*transistionmap'
        approximation = approximation+(gmera2approx|>CrystalDenseMap)
        approxdict[:second] = gmera2approx
        filledisometrydict[:second] = gmera2[:filledisometry]|>CrystalDenseMap
        emptyisometrydict[:second] = gmera2[:emptyisometry]|>CrystalDenseMap
        if haskey(gmera2,:correlations)
            transistionmap = transistionmap*gmera2[:courierisometry]
            transistionmapdict[:second] = transistionmap
            correlationsdict[:second] = gmera2[:correlations]|>CrystalDenseMap
            @info ("3rd gmera step...")
            gmera3 = @time gmerastep(rgblockedcorrelations,gmera2[:correlations],thirdcenterlist, selectionstrategy3)
            @info ("3rd gmera approximation to correlations...")
            gmera3approx = transistionmap*gmera3[:emptyisometry]*gmera3[:emptyisometry]'*transistionmap'
            approximation = approximation+(gmera3approx|>CrystalDenseMap)
            approxdict[:third] = gmera3approx
            filledisometrydict[:third] = gmera3[:filledisometry]|>CrystalDenseMap
            emptyisometrydict[:third] = gmera3[:emptyisometry]|>CrystalDenseMap
            if haskey(gmera3,:correlations)
                transistionmap = transistionmap*gmera3[:courierisometry]
                transistionmapdict[:third] = transistionmap
                correlationsdict[:third] = gmera3[:correlations]|>CrystalDenseMap
                @info ("final gmera step...")
                gmera4 = @time gmerastep(rgblockedcorrelations,gmera3[:correlations],finalcenterlist, selectionstrategy4)
                @info ("4th gmera approximation to correlations...")
                gmera4approx = transistionmap*gmera4[:emptyisometry]*gmera4[:emptyisometry]'*transistionmap'
                approximation = approximation+(gmera4approx|>CrystalDenseMap)
                approxdict[:forth] = gmera4approx
                filledisometrydict[:forth] = gmera4[:filledisometry]|>CrystalDenseMap
                emptyisometrydict[:forth] = gmera4[:emptyisometry]|>CrystalDenseMap
                if haskey(gmera4,:correlations)
                    transistionmap = transistionmap*gmera4[:courierisometry]
                    transistionmapdict[:forth] = transistionmap
                    correlationsdict[:forth] = gmera4[:correlations]|>CrystalDenseMap
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

function rg(correlations,selectionstrategy,RGstep)
    step = 1
    @info("implementing "*string(step)*" rg step")
    gmeraresult = gmera(correlations,idmap(correlations|>getinspace),selectionstrategy,RGstep==step)
    rgresult = Dict(Symbol("rg"*string(step))=>gmeraresult)
    step+=1
    while (haskey(gmeraresult,:correlations)) 
        if ((gmeraresult[:correlations]|>getinspace|>getcrystal).sizes[1]==3) || ((gmeraresult[:correlations]|>getinspace|>getcrystal).sizes[1]==2)
            @info("cannot further do RG")
            return rgresult
        else 
            @info("implementing "*string(step)*" rg step and the remaining systemsize "*string(gmeraresult[:correlations]|>getinspace|>getcrystal))
            gmeraresult = gmera(gmeraresult[:correlations],gmeraresult[:transistionmaps][:forth],selectionstrategy,RGstep==step)
            rgresult[Symbol("rg"*string(step))] = gmeraresult
            step+=1
        end
    end
    return rgresult
end

fiodir("../../../../../../../ZERData/GMERA/trivial/t_a_neg_04_t_b_neg_06/size96/localsize8/RGstep1")

# hyperparemeter
modeselectionstrategy = modeselection1stbycountthenbythreshold(1,0.002)
t_a=-0.4
t_b=-0.6
size=96
localsize=8
RGstep=1

# norm(regioncorrelations(blockedcorrelations,quantize(blockedcorrelations|>getinspace|>getcrystal|>getunitcell,1))[quantize([offset for offset in blockedcorrelations|>getinspace|>getcrystal|>getunitcell][1],1),:])
# implement the process
H, correlations = H_trivial(size,ComplexF64(-1),t_a,t_b)
fiosave(H, name="hamiltonian")
fiosave(correlations, name="correlations")
blockedcorrelations = startingcorrelations(correlations,localsize)
# visualize(regioncorrelations(blockedcorrelations,quantize(blockedcorrelations|>getinspace|>getcrystal|>getunitcell,1))|>eigspech)

rgresults = rg(blockedcorrelations,modeselectionstrategy,RGstep)
rgapprox = rgapproximation(rgresults)
fiosave(rgapprox, name="rgapprox")
saveisometry(rgresults)
rgapproxCrystalFock = rgapprox|>CrystalFockMap
traceL1norm = focktraceL1norm(blockedcorrelations-rgapproxCrystalFock,size*size*2)
fiosave(traceL1norm, name="traceL1norm")

visualize(rgresults[:rg1][:correlationsresults][:first]|>getinspace|>getcrystal|>getunitcell)
