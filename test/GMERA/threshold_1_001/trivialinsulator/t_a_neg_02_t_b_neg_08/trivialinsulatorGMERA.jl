using Pkg
Pkg.resolve()
Pkg.instantiate()
# Pkg.rm("PlotlyJS")
Pkg.Registry.add(RegistrySpec(url="https://github.com/wildart/BoffinStuff.git"))
Pkg.add("SmithNormalForm")
Pkg.activate("../../../../../")
Pkg.resolve()
Pkg.instantiate()
using Zipper
using SmithNormalForm, LinearAlgebra, OrderedCollections, SparseArrays, Combinatorics
using Plots

setmaxthreads(8)

function focktraceL1norm(fockmap,systemsize)
    @info("Calculating L1norm...")
    modeevalpairs = fockmap|>eigspech|>geteigenvalues
    return(sum([abs(modeevalpair[2]) for modeevalpair in modeevalpairs])/systemsize)
end

function focktraceL2norm(fockmap,systemsize)
    @info("Calculating L2norm...")
    return (real(sqrt(tr(fockmap*fockmap'|>rep)/systemsize)))
end


function rgapproximation(rgresults)
    finalrg = Symbol("rg"*string(length(rgresults)))
    @info("Performing rg approximation...")
    if haskey(rgresults[finalrg],:correlations)
        core = @time groupbands(rgresults[finalrg][:correlations]|>crystalspectrum, :filled=> v -> v < 1e-5, :empty => v -> v > 1e-5)
        coreproj = core[:empty]|>crystalprojector
        coreapprox = rgresults[finalrg][:transistionmap]*coreproj*rgresults[finalrg][:transistionmap]'
        @info("Return rg approximation with core ...")
        return sum([rgresult[2][:approximation] for rgresult in rgresults])+coreapprox
    else
        @info("Return rg approximation without core ...")
        return sum([rgresult[2][:approximation] for rgresult in rgresults])
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
    kspace = convert(MomentumSpace, triangular)

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

t_a=-0.2
t_b=-0.8
modeselectionstrategy = modeselection1stbycountthenbythreshold(1,0.001)


# H_size12, correlations_size12 = H_trivial(12,ComplexF64(-1),t_a,t_b)
fiodir("../../../../../ZERData/GMERA/threshold_1_001/trivial/t_a_neg_02_t_b_neg_08")
# fiosave(H_size12, name="hamiltonian_size12")
# fiosave(correlations_size12, name="correlations_size12")

# blockedcorrelationsloc2size12 = startingcorrelations(correlations_size12,2)
# fiosave(blockedcorrelationsloc2size12, name="blockedcorrelationsloc2size12")
# rgresultsloc2size12 = rg(blockedcorrelationsloc2size12,modeselectionstrategy)
# noofrgsteploc2size12 = length(rgresultsloc2size12)
# fiosave(noofrgsteploc2size12, name="noofrgsteploc2size12")
# noofcoremodesloc2size12 = countnoofcoremodes(rgresultsloc2size12)
# fiosave(noofcoremodesloc2size12, name="noofcoremodesloc2size12")
# rgapproxloc2size12 = rgapproximation(rgresultsloc2size12)
# fiosave(rgapproxloc2size12, name="rgapproxloc2size12")

# # errorloc2size12L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc2size12-rgapproxloc2size12),288)
# # errorloc2size12L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc2size12-rgapproxloc2size12),288)
# # fiosave(errorloc2size12L1norm, name="errorloc2size12L1norm")
# # fiosave(errorloc2size12L2norm, name="errorloc2size12L2norm")


# blockedcorrelationsloc3size12 = startingcorrelations(correlations_size12,3)
# fiosave(blockedcorrelationsloc3size12, name="blockedcorrelationsloc3size12")
# rgresultsloc3size12 = rg(blockedcorrelationsloc3size12,modeselectionstrategy)
# noofrgsteploc3size12 = length(rgresultsloc3size12)
# fiosave(noofrgsteploc3size12, name="noofrgsteploc3size12")
# noofcoremodesloc3size12 = countnoofcoremodes(rgresultsloc3size12)
# fiosave(noofcoremodesloc3size12, name="noofcoremodesloc3size12")
# rgapproxloc3size12 = rgapproximation(rgresultsloc3size12)
# fiosave(rgapproxloc3size12, name="rgapproxloc3size12")

# # errorloc3size12L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc3size12-rgapproxloc3size12),288)
# # errorloc3size12L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc3size12-rgapproxloc3size12),288)
# # fiosave(errorloc3size12L1norm, name="errorloc3size12L1norm")
# # fiosave(errorloc3size12L2norm, name="errorloc3size12L2norm")


# H_size24, correlations_size24 = H_trivial(24,ComplexF64(-1),t_a,t_b)
# fiosave(H_size24, name="hamiltonian_size24")
# fiosave(correlations_size24, name="correlations_size24")

# blockedcorrelationsloc2size24 = startingcorrelations(correlations_size24,2)
# fiosave(blockedcorrelationsloc2size24, name="blockedcorrelationsloc2size24")
# rgresultsloc2size24 = rg(blockedcorrelationsloc2size24,modeselectionstrategy)
# noofrgsteploc2size24 = length(rgresultsloc2size24)
# fiosave(noofrgsteploc2size24, name="noofrgsteploc2size24")
# noofcoremodesloc2size24 = countnoofcoremodes(rgresultsloc2size24)
# fiosave(noofcoremodesloc2size24, name="noofcoremodesloc2size24")
# rgapproxloc2size24 = rgapproximation(rgresultsloc2size24)
# fiosave(rgapproxloc2size24, name="rgapproxloc2size24")

# # errorloc2size24L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc2size24-rgapproxloc2size24),1152)
# # errorloc2size24L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc2size24-rgapproxloc2size24),1152)
# # fiosave(errorloc2size24L1norm, name="errorloc2size24L1norm")
# # fiosave(errorloc2size24L2norm, name="errorloc2size24L2norm")

# blockedcorrelationsloc3size24 = startingcorrelations(correlations_size24,3)
# fiosave(blockedcorrelationsloc3size24, name="blockedcorrelationsloc3size24")
# rgresultsloc3size24 = rg(blockedcorrelationsloc3size24,modeselectionstrategy)
# noofrgsteploc3size24 = length(rgresultsloc3size24)
# fiosave(noofrgsteploc3size24, name="noofrgsteploc3size24")
# noofcoremodesloc3size24 = countnoofcoremodes(rgresultsloc3size24)
# fiosave(noofcoremodesloc3size24, name="noofcoremodesloc3size24")
# rgapproxloc3size24 = rgapproximation(rgresultsloc3size24)
# fiosave(rgapproxloc3size24, name="rgapproxloc3size24")

# # errorloc3size24L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc3size24-rgapproxloc3size24),1152)
# # errorloc3size24L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc3size24-rgapproxloc3size24),1152)
# # fiosave(errorloc3size24L1norm, name="errorloc3size24L1norm")
# # fiosave(errorloc3size24L2norm, name="errorloc3size24L2norm")

# blockedcorrelationsloc4size24 = startingcorrelations(correlations_size24,4)
# fiosave(blockedcorrelationsloc4size24, name="blockedcorrelationsloc4size24")
# rgresultsloc4size24 = rg(blockedcorrelationsloc4size24,modeselectionstrategy)
# noofrgsteploc4size24 = length(rgresultsloc4size24)
# fiosave(noofrgsteploc4size24, name="noofrgsteploc4size24")
# noofcoremodesloc4size24 = countnoofcoremodes(rgresultsloc4size24)
# fiosave(noofcoremodesloc4size24, name="noofcoremodesloc4size24")
# rgapproxloc4size24 = rgapproximation(rgresultsloc4size24)
# fiosave(rgapproxloc4size24, name="rgapproxloc4size24")

# # errorloc4size24L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc4size24-rgapproxloc4size24),1152)
# # errorloc4size24L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc4size24-rgapproxloc4size24),1152)
# # fiosave(errorloc4size24L1norm, name="errorloc4size24L1norm")
# # fiosave(errorloc4size24L2norm, name="errorloc4size24L2norm")

# blockedcorrelationsloc6size24 = startingcorrelations(correlations_size24,6)
# fiosave(blockedcorrelationsloc6size24, name="blockedcorrelationsloc6size24")
# rgresultsloc6size24 = rg(blockedcorrelationsloc6size24,modeselectionstrategy)
# noofrgsteploc6size24 = length(rgresultsloc6size24)
# fiosave(noofrgsteploc6size24, name="noofrgsteploc6size24")
# noofcoremodesloc6size24 = countnoofcoremodes(rgresultsloc6size24)
# fiosave(noofcoremodesloc6size24, name="noofcoremodesloc6size24")
# rgapproxloc6size24 = rgapproximation(rgresultsloc6size24)
# fiosave(rgapproxloc6size24, name="rgapproxloc6size24")

# # errorloc6size24L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc6size24-rgapproxloc6size24),1152)
# # errorloc6size24L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc6size24-rgapproxloc6size24),1152)
# # fiosave(errorloc6size24L1norm, name="errorloc6size24L1norm")
# # fiosave(errorloc6size24L2norm, name="errorloc6size24L2norm")


# H_size48, correlations_size48 = H_trivial(48,ComplexF64(-1),t_a,t_b)
# fiosave(H_size48, name="hamiltonian_size48")
# fiosave(correlations_size48, name="correlations_size48")

# blockedcorrelationsloc2size48 = startingcorrelations(correlations_size48,2)
# fiosave(blockedcorrelationsloc2size48, name="blockedcorrelationsloc2size48")
# rgresultsloc2size48 = rg(blockedcorrelationsloc2size48,modeselectionstrategy)
# noofrgsteploc2size48 = length(rgresultsloc2size48)
# fiosave(noofrgsteploc2size48, name="noofrgsteploc2size48")
# noofcoremodesloc2size48 = countnoofcoremodes(rgresultsloc2size48)
# fiosave(noofcoremodesloc2size48, name="noofcoremodesloc2size48")
# rgapproxloc2size48 = rgapproximation(rgresultsloc2size48)
# fiosave(rgapproxloc2size48, name="rgapproxloc2size48")

# # errorloc2size48L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc2size48-rgapproxloc2size48),4608)
# # errorloc2size48L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc2size48-rgapproxloc2size48),4608)
# # fiosave(errorloc2size48L1norm, name="errorloc2size48L1norm")
# # fiosave(errorloc2size48L2norm, name="errorloc2size48L2norm")

# blockedcorrelationsloc3size48 = startingcorrelations(correlations_size48,3)
# fiosave(blockedcorrelationsloc3size48, name="blockedcorrelationsloc3size48")
# rgresultsloc3size48 = rg(blockedcorrelationsloc3size48,modeselectionstrategy)
# noofrgsteploc3size48 = length(rgresultsloc3size48)
# fiosave(noofrgsteploc3size48, name="noofrgsteploc3size48")
# noofcoremodesloc3size48 = countnoofcoremodes(rgresultsloc3size48)
# fiosave(noofcoremodesloc3size48, name="noofcoremodesloc3size48")
# rgapproxloc3size48 = rgapproximation(rgresultsloc3size48)
# fiosave(rgapproxloc3size48, name="rgapproxloc3size48")

# # errorloc3size48L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc3size48-rgapproxloc3size48),4608)
# # errorloc3size48L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc3size48-rgapproxloc3size48),4608)
# # fiosave(errorloc3size48L1norm, name="errorloc3size48L1norm")
# # fiosave(errorloc3size48L2norm, name="errorloc3size48L2norm")


# blockedcorrelationsloc4size48 = startingcorrelations(correlations_size48,4)
# fiosave(blockedcorrelationsloc4size48, name="blockedcorrelationsloc4size48")
# rgresultsloc4size48 = rg(blockedcorrelationsloc4size48,modeselectionstrategy)
# noofrgsteploc4size48 = length(rgresultsloc4size48)
# fiosave(noofrgsteploc4size48, name="noofrgsteploc4size48")
# noofcoremodesloc4size48 = countnoofcoremodes(rgresultsloc4size48)
# fiosave(noofcoremodesloc4size48, name="noofcoremodesloc4size48")
# rgapproxloc4size48 = rgapproximation(rgresultsloc4size48)
# fiosave(rgapproxloc4size48, name="rgapproxloc4size48")

# # errorloc4size48L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc4size48-rgapproxloc4size48),4608)
# # errorloc4size48L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc4size48-rgapproxloc4size48),4608)
# # fiosave(errorloc4size48L1norm, name="errorloc4size48L1norm")
# # fiosave(errorloc4size48L2norm, name="errorloc4size48L2norm")

# blockedcorrelationsloc6size48 = startingcorrelations(correlations_size48,6)
# fiosave(blockedcorrelationsloc6size48, name="blockedcorrelationsloc6size48")
# rgresultsloc6size48 = rg(blockedcorrelationsloc6size48,modeselectionstrategy)
# noofrgsteploc6size48 = length(rgresultsloc6size48)
# fiosave(noofrgsteploc6size48, name="noofrgsteploc6size48")
# noofcoremodesloc6size48 = countnoofcoremodes(rgresultsloc6size48)
# fiosave(noofcoremodesloc6size48, name="noofcoremodesloc6size48")
# rgapproxloc6size48 = rgapproximation(rgresultsloc6size48)
# fiosave(rgapproxloc6size48, name="rgapproxloc6size48")

# # errorloc6size48L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc6size48-rgapproxloc6size48),4608)
# # errorloc6size48L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc6size48-rgapproxloc6size48),4608)
# # fiosave(errorloc6size48L1norm, name="errorloc6size48L1norm")
# # fiosave(errorloc6size48L2norm, name="errorloc6size48L2norm")

# blockedcorrelationsloc8size48 = startingcorrelations(correlations_size48,8)
# fiosave(blockedcorrelationsloc8size48, name="blockedcorrelationsloc8size48")
# rgresultsloc8size48 = rg(blockedcorrelationsloc8size48,modeselectionstrategy)
# noofrgsteploc8size48 = length(rgresultsloc8size48)
# fiosave(noofrgsteploc8size48, name="noofrgsteploc8size48")
# noofcoremodesloc8size48 = countnoofcoremodes(rgresultsloc8size48)
# fiosave(noofcoremodesloc8size48, name="noofcoremodesloc8size48")
# rgapproxloc8size48 = rgapproximation(rgresultsloc8size48)
# fiosave(rgapproxloc8size48, name="rgapproxloc8size48")

# # errorloc8size48L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc8size48-rgapproxloc8size48),4608)
# # errorloc8size48L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc8size48-rgapproxloc8size48),4608)
# # fiosave(errorloc8size48L1norm, name="errorloc8size48L1norm")
# # fiosave(errorloc8size48L2norm, name="errorloc8size48L2norm")

# blockedcorrelationsloc12size48 = startingcorrelations(correlations_size48,12)
# fiosave(blockedcorrelationsloc12size48, name="blockedcorrelationsloc12size48")
# rgresultsloc12size48 = rg(blockedcorrelationsloc12size48,modeselectionstrategy)
# noofrgsteploc12size48 = length(rgresultsloc12size48)
# fiosave(noofrgsteploc12size48, name="noofrgsteploc12size48")
# noofcoremodesloc12size48 = countnoofcoremodes(rgresultsloc12size48)
# fiosave(noofcoremodesloc12size48, name="noofcoremodesloc12size48")
# rgapproxloc12size48 = rgapproximation(rgresultsloc12size48)
# fiosave(rgapproxloc12size48, name="rgapproxloc12size48")

# errorloc12size48L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc12size48-rgapproxloc12size48),4608)
# errorloc12size48L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc12size48-rgapproxloc12size48),4608)
# fiosave(errorloc12size48L1norm, name="errorloc12size48L1norm")
# fiosave(errorloc12size48L2norm, name="errorloc12size48L2norm")


H_size96, correlations_size96 = H_trivial(96,ComplexF64(-1),t_a,t_b)
fiosave(H_size96, name="hamiltonian_size96")
fiosave(correlations_size96, name="correlations_size96")

# blockedcorrelationsloc2size96 = startingcorrelations(correlations_size96,2)
# fiosave(blockedcorrelationsloc2size96, name="blockedcorrelationsloc2size96")
# rgresultsloc2size96 = rg(blockedcorrelationsloc2size96,modeselectionstrategy)
# noofrgsteploc2size96 = length(rgresultsloc2size96)
# fiosave(noofrgsteploc2size96, name="noofrgsteploc2size96")
# noofcoremodesloc2size96 = countnoofcoremodes(rgresultsloc2size96)
# fiosave(noofcoremodesloc2size96, name="noofcoremodesloc2size96")
# rgapproxloc2size96 = rgapproximation(rgresultsloc2size96)
# fiosave(rgapproxloc2size96, name="rgapproxloc2size96")

# errorloc2size96L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc2size96-rgapproxloc2size96),4608*4)
# errorloc2size96L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc2size96-rgapproxloc2size96),4608*4)
# fiosave(errorloc2size96L1norm, name="errorloc2size96L1norm")
# fiosave(errorloc2size96L2norm, name="errorloc2size96L2norm")

# @info("system size 96,local size 3")
# blockedcorrelationsloc3size96 = startingcorrelations(correlations_size96,3)
# fiosave(blockedcorrelationsloc3size96, name="blockedcorrelationsloc3size96")
# rgresultsloc3size96 = rg(blockedcorrelationsloc3size96,modeselectionstrategy)
# noofrgsteploc3size96 = length(rgresultsloc3size96)
# fiosave(noofrgsteploc3size96, name="noofrgsteploc3size96")
# noofcoremodesloc3size96 = countnoofcoremodes(rgresultsloc3size96)
# fiosave(noofcoremodesloc3size96, name="noofcoremodesloc3size96")
# rgapproxloc3size96 = rgapproximation(rgresultsloc3size96)
# @info("finish rg approximation")
# fiosave(rgapproxloc3size96, name="rgapproxloc3size96")

# # errorloc3size96L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc3size96-rgapproxloc3size96),4608*4)
# # errorloc3size96L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc3size96-rgapproxloc3size96),4608*4)
# # fiosave(errorloc3size96L1norm, name="errorloc3size96L1norm")
# # fiosave(errorloc3size96L2norm, name="errorloc3size96L2norm")

# @info("system size 96,local size 4")
# blockedcorrelationsloc4size96 = startingcorrelations(correlations_size96,4)
# fiosave(blockedcorrelationsloc4size96, name="blockedcorrelationsloc4size96")
# rgresultsloc4size96 = rg(blockedcorrelationsloc4size96,modeselectionstrategy)
# noofrgsteploc4size96 = length(rgresultsloc4size96)
# fiosave(noofrgsteploc4size96, name="noofrgsteploc4size96")
# noofcoremodesloc4size96 = countnoofcoremodes(rgresultsloc4size96)
# fiosave(noofcoremodesloc4size96, name="noofcoremodesloc4size96")
# rgapproxloc4size96 = rgapproximation(rgresultsloc4size96)
# @info("finish rg approximation")
# fiosave(rgapproxloc4size96, name="rgapproxloc4size96")

# # errorloc4size96L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc4size96-rgapproxloc4size96),4608*4)
# # errorloc4size96L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc4size96-rgapproxloc4size96),4608*4)
# # fiosave(errorloc4size96L1norm, name="errorloc4size96L1norm")
# # fiosave(errorloc4size96L2norm, name="errorloc4size96L2norm")

# @info("system size 96,local size 6")
# blockedcorrelationsloc6size96 = startingcorrelations(correlations_size96,6)
# fiosave(blockedcorrelationsloc6size96, name="blockedcorrelationsloc6size96")
# rgresultsloc6size96 = rg(blockedcorrelationsloc6size96,modeselectionstrategy)
# noofrgsteploc6size96 = length(rgresultsloc6size96)
# fiosave(noofrgsteploc6size96, name="noofrgsteploc6size96")
# noofcoremodesloc6size96 = countnoofcoremodes(rgresultsloc6size96)
# fiosave(noofcoremodesloc6size96, name="noofcoremodesloc6size96")
# rgapproxloc6size96 = rgapproximation(rgresultsloc6size96)
# @info("finish rg approximation")
# fiosave(rgapproxloc6size96, name="rgapproxloc6size96")

# errorloc6size96L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc6size96-rgapproxloc6size96),4608*4)
# errorloc6size96L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc6size96-rgapproxloc6size96),4608*4)
# fiosave(errorloc6size96L1norm, name="errorloc6size96L1norm")
# fiosave(errorloc6size96L2norm, name="errorloc6size96L2norm")

# @info("system size 96,local size 8")
# blockedcorrelationsloc8size96 = startingcorrelations(correlations_size96,8)
# fiosave(blockedcorrelationsloc8size96, name="blockedcorrelationsloc8size96")
# rgresultsloc8size96 = rg(blockedcorrelationsloc8size96,modeselectionstrategy)
# noofrgsteploc8size96 = length(rgresultsloc8size96)
# fiosave(noofrgsteploc8size96, name="noofrgsteploc8size96")
# noofcoremodesloc8size96 = countnoofcoremodes(rgresultsloc8size96)
# fiosave(noofcoremodesloc8size96, name="noofcoremodesloc8size96")
# rgapproxloc8size96 = rgapproximation(rgresultsloc8size96)
# @info("finish rg approximation")
# fiosave(rgapproxloc8size96, name="rgapproxloc8size96")

# # errorloc8size96L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc8size96-rgapproxloc8size96),4608*4)
# # errorloc8size96L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc8size96-rgapproxloc8size96),4608*4)
# # fiosave(errorloc8size96L1norm, name="errorloc8size96L1norm")
# # fiosave(errorloc8size96L2norm, name="errorloc8size96L2norm")

# @info("system size 96,local size 12")
# blockedcorrelationsloc12size96 = startingcorrelations(correlations_size96,12)
# fiosave(blockedcorrelationsloc12size96, name="blockedcorrelationsloc12size96")
# rgresultsloc12size96 = rg(blockedcorrelationsloc12size96,modeselectionstrategy)
# noofrgsteploc12size96 = length(rgresultsloc12size96)
# fiosave(noofrgsteploc12size96, name="noofrgsteploc12size96")
# noofcoremodesloc12size96 = countnoofcoremodes(rgresultsloc12size96)
# fiosave(noofcoremodesloc12size96, name="noofcoremodesloc12size96")
# rgapproxloc12size96 = rgapproximation(rgresultsloc12size96)
# @info("finish rg approximation")
# fiosave(rgapproxloc12size96, name="rgapproxloc12size96")

# # errorloc12size96L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc12size96-rgapproxloc12size96),4608*4)
# # errorloc12size96L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc12size96-rgapproxloc12size96),4608*4)
# # fiosave(errorloc12size96L1norm, name="errorloc12size96L1norm")
# # fiosave(errorloc12size96L2norm, name="errorloc12size96L2norm")

# @info("system size 96,local size 16")
# blockedcorrelationsloc16size96 = startingcorrelations(correlations_size96,16)
# fiosave(blockedcorrelationsloc16size96, name="blockedcorrelationsloc16size96")
# rgresultsloc16size96 = rg(blockedcorrelationsloc16size96,modeselectionstrategy)
# noofrgsteploc16size96 = length(rgresultsloc16size96)
# fiosave(noofrgsteploc16size96, name="noofrgsteploc16size96")
# noofcoremodesloc16size96 = countnoofcoremodes(rgresultsloc16size96)
# fiosave(noofcoremodesloc16size96, name="noofcoremodesloc16size96")
# rgapproxloc16size96 = rgapproximation(rgresultsloc16size96)
# @info("finish rg approximation")
# fiosave(rgapproxloc16size96, name="rgapproxloc16size96")

# # errorloc16size96L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc16size96-rgapproxloc16size96),4608*4)
# # errorloc16size96L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc16size96-rgapproxloc16size96),4608*4)
# # fiosave(errorloc16size96L1norm, name="errorloc16size96L1norm")
# # fiosave(errorloc16size96L2norm, name="errorloc16size96L2norm")

@info("system size 96,local size 24")
blockedcorrelationsloc24size96 = startingcorrelations(correlations_size96,24)
fiosave(blockedcorrelationsloc24size96, name="blockedcorrelationsloc24size96")
rgresultsloc24size96 = rg(blockedcorrelationsloc24size96,modeselectionstrategy)
noofrgsteploc24size96 = length(rgresultsloc24size96)
fiosave(noofrgsteploc24size96 , name="noofrgsteploc24size96 ")
noofcoremodesloc24size96 = countnoofcoremodes(rgresultsloc24size96)
fiosave(noofcoremodesloc24size96, name="noofcoremodesloc24size96")
rgapproxloc24size96 = rgapproximation(rgresultsloc24size96)
@info("finish rg approximation")
fiosave(rgapproxloc24size96, name="rgapproxloc24size96")

# errorloc24size96L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc24size96-rgapproxloc24size96),4608*4)
# errorloc24size96L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc24size96-rgapproxloc24size96),4608*4)
# fiosave(errorloc24size96L1norm, name="errorloc24size96L1norm")
# fiosave(errorloc24size96L2norm, name="errorloc24size96L2norm")

H_size192, correlations_size192 = H_trivial(192,ComplexF64(-1),t_a,t_b)
fiosave(H_size192, name="hamiltonian_size192")
fiosave(correlations_size192, name="correlations_size192")

@info("system size 192,local size 2")
blockedcorrelationsloc2size192 = startingcorrelations(correlations_size192,2)
fiosave(blockedcorrelationsloc2size192, name="blockedcorrelationsloc2size192")
rgresultsloc2size192 = rg(blockedcorrelationsloc2size192,modeselectionstrategy)
noofrgsteploc2size192 = length(rgresultsloc2size192)
fiosave(noofrgsteploc2size192, name="noofrgsteploc2size192")
noofcoremodesloc2size192 = countnoofcoremodes(rgresultsloc2size192)
fiosave(noofcoremodesloc2size192, name="noofcoremodesloc2size192")
rgapproxloc2size192 = rgapproximation(rgresultsloc2size192)
@info("finish rg approximation")
fiosave(rgapproxloc2size192, name="rgapproxloc2size192")

# errorloc2size192L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc2size192-rgapproxloc2size192),4608*4*4)
# errorloc2size192L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc2size192-rgapproxloc2size192),4608*4*4)
# fiosave(errorloc2size192L1norm, name="errorloc2size192L1norm")
# fiosave(errorloc2size192L2norm, name="errorloc2size192L2norm")

# @info("system size 192,local size 3")
# blockedcorrelationsloc3size192 = startingcorrelations(correlations_size192,3)
# fiosave(blockedcorrelationsloc3size192, name="blockedcorrelationsloc3size192")
# rgresultsloc3size192 = rg(blockedcorrelationsloc3size192,modeselectionstrategy)
# noofrgsteploc3size192 = length(rgresultsloc3size192)
# fiosave(noofrgsteploc3size192, name="noofrgsteploc3size192")
# noofcoremodesloc3size192 = countnoofcoremodes(rgresultsloc3size192)
# fiosave(noofcoremodesloc3size192, name="noofcoremodesloc3size192")
# rgapproxloc3size192 = rgapproximation(rgresultsloc3size192)
# @info("finish rg approximation")
# fiosave(rgapproxloc3size192, name="rgapproxloc3size192")

# # errorloc3size192L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc3size192-rgapproxloc3size192),4608*4*4)
# # errorloc3size192L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc3size192-rgapproxloc3size192),4608*4*4)
# # fiosave(errorloc3size192L1norm, name="errorloc3size192L1norm")
# # fiosave(errorloc3size192L2norm, name="errorloc3size192L2norm")

# @info("system size 192,local size 4")
# blockedcorrelationsloc4size192 = startingcorrelations(correlations_size192,4)
# fiosave(blockedcorrelationsloc4size192, name="blockedcorrelationsloc4size192")
# rgresultsloc4size192 = rg(blockedcorrelationsloc4size192,modeselectionstrategy)
# noofrgsteploc4size192 = length(rgresultsloc4size192)
# fiosave(noofrgsteploc4size192, name="noofrgsteploc4size192")
# noofcoremodesloc4size192 = countnoofcoremodes(rgresultsloc4size192)
# fiosave(noofcoremodesloc4size192, name="noofcoremodesloc4size192")
# rgapproxloc4size192 = rgapproximation(rgresultsloc4size192)
# @info("finish rg approximation")
# fiosave(rgapproxloc4size192, name="rgapproxloc4size192")

# # errorloc4size192L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc4size192-rgapproxloc4size192),4608*4*4)
# # errorloc4size192L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc4size192-rgapproxloc4size192),4608*4*4)
# # fiosave(errorloc4size192L1norm, name="errorloc4size192L1norm")
# # fiosave(errorloc4size192L2norm, name="errorloc4size192L2norm")

# @info("system size 192,local size 6")
# blockedcorrelationsloc6size192 = startingcorrelations(correlations_size192,6)
# fiosave(blockedcorrelationsloc6size192 , name="blockedcorrelationsloc6size192")
# rgresultsloc6size192 = rg(blockedcorrelationsloc6size192,modeselectionstrategy)
# noofrgsteploc6size192 = length(rgresultsloc6size192)
# fiosave(noofrgsteploc6size192, name="noofrgsteploc6size192")
# noofcoremodesloc6size192 = countnoofcoremodes(rgresultsloc6size192)
# fiosave(noofcoremodesloc6size192, name="noofcoremodesloc6size192")
# rgapproxloc6size192 = rgapproximation(rgresultsloc6size192)
# @info("finish rg approximation")
# fiosave(rgapproxloc6size192, name="rgapproxloc6size192")

# # errorloc6size192L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc6size192-rgapproxloc6size192),4608*4*4)
# # errorloc6size192L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc6size192-rgapproxloc6size192),4608*4*4)
# # fiosave(errorloc6size192L1norm, name="errorloc6size192L1norm")
# # fiosave(errorloc6size192L2norm, name="errorloc6size192L2norm")

# @info("system size 192,local size 8")
# blockedcorrelationsloc8size192 = startingcorrelations(correlations_size192,8)
# fiosave(blockedcorrelationsloc8size192, name="blockedcorrelationsloc8size192")
# rgresultsloc8size192 = rg(blockedcorrelationsloc8size192,modeselectionstrategy)
# noofrgsteploc8size192 = length(rgresultsloc8size192)
# fiosave(noofrgsteploc8size192, name="noofrgsteploc8size192")
# noofcoremodesloc8size192 = countnoofcoremodes(rgresultsloc8size192)
# fiosave(noofcoremodesloc8size96, name="noofcoremodesloc8size96")
# rgapproxloc8size192 = rgapproximation(rgresultsloc8size192)
# @info("finish rg approximation")
# fiosave(rgapproxloc8size96, name="rgapproxloc8size96")

# # errorloc8size192L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc8size192-rgapproxloc8size192),4608*4*4)
# # errorloc8size192L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc8size192-rgapproxloc8size192),4608*4*4)
# # fiosave(errorloc8size96L1norm, name="errorloc8size96L1norm")
# # fiosave(errorloc8size96L2norm, name="errorloc8size96L2norm")

# @info("system size 192,local size 12")
# blockedcorrelationsloc12size192 = startingcorrelations(correlations_size192,12)
# fiosave(blockedcorrelationsloc12size192, name="blockedcorrelationsloc12size192")
# rgresultsloc12size192 = rg(blockedcorrelationsloc12size192,modeselectionstrategy)
# noofrgsteploc12size192 = length(rgresultsloc12size192)
# fiosave(noofrgsteploc12size192, name="noofrgsteploc12size192")
# noofcoremodesloc12size192 = countnoofcoremodes(rgresultsloc12size192)
# fiosave(noofcoremodesloc12size192, name="noofcoremodesloc12size192")
# rgapproxloc12size192 = rgapproximation(rgresultsloc12size192)
# @info("finish rg approximation")
# fiosave(rgapproxloc12size192, name="rgapproxloc12size192")

# # errorloc12size192L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc12size192-rgapproxloc12size192),4608*4*4)
# # errorloc12size192L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc12size192-rgapproxloc12size192),4608*4*4)
# # fiosave(errorloc12size192L1norm, name="errorloc12size192L1norm")
# # fiosave(errorloc12size192L2norm, name="errorloc12size192L2norm")

# @info("system size 192,local size 16")
# blockedcorrelationsloc16size192 = startingcorrelations(correlations_size192,16)
# fiosave(blockedcorrelationsloc16size192, name="blockedcorrelationsloc16size192")
# rgresultsloc16size192 = rg(blockedcorrelationsloc16size192,modeselectionstrategy)
# noofrgsteploc16size192 = length(rgresultsloc16size192)
# fiosave(noofrgsteploc16size192, name="noofrgsteploc16size192")
# noofcoremodesloc16size192 = countnoofcoremodes(rgresultsloc16size192)
# fiosave(noofcoremodesloc16size192, name="noofcoremodesloc16size192")
# rgapproxloc16size192 = rgapproximation(rgresultsloc16size192)
# @info("finish rg approximation")
# fiosave(rgapproxloc16size192, name="rgapproxloc16size192")

# # errorloc16size192L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc16size192-rgapproxloc16size192),4608*4*4)
# # errorloc16size192L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc16size192-rgapproxloc16size192),4608*4*4)
# # fiosave(errorloc16size192L1norm, name="errorloc16size192L1norm")
# # fiosave(errorloc16size192L2norm, name="errorloc16size192L2norm")

# @info("system size 192,local size 24")
# blockedcorrelationsloc24size192 = startingcorrelations(correlations_size192,24)
# fiosave(blockedcorrelationsloc24size192, name="blockedcorrelationsloc24size192")
# rgresultsloc24size192 = rg(blockedcorrelationsloc24size192,modeselectionstrategy)
# noofrgsteploc24size192 = length(rgresultsloc24size192)
# fiosave(noofrgsteploc24size192, name="noofrgsteploc24size192")
# noofcoremodesloc24size192 = countnoofcoremodes(rgresultsloc24size192)
# fiosave(noofcoremodesloc24size192, name="noofcoremodesloc24size192")
# rgapproxloc24size192 = rgapproximation(rgresultsloc24size192)
# @info("finish rg approximation")
# fiosave(rgapproxloc24size192, name="rgapproxloc24size192")

# # errorloc24size192L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc24size192-rgapproxloc24size192),4608*4*4)
# # errorloc24size192L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc24size192-rgapproxloc24size192),4608*4*4)
# # fiosave(errorloc24size192L1norm, name="errorloc24size192L1norm")
# # fiosave(errorloc24size192L2norm, name="errorloc24size192L2norm")

# @info("system size 192,local size 32")
# blockedcorrelationsloc32size192 = startingcorrelations(correlations_size192,32)
# fiosave(blockedcorrelationsloc32size192, name="blockedcorrelationsloc32size192")
# rgresultsloc32size192 = rg(blockedcorrelationsloc32size192,modeselectionstrategy)
# noofrgsteploc32size192 = length(rgresultsloc32size192)
# fiosave(noofrgsteploc32size192, name="noofrgsteploc32size192")
# noofcoremodesloc32size192 = countnoofcoremodes(rgresultsloc32size192)
# fiosave(noofcoremodesloc32size192, name="noofcoremodesloc32size192")
# rgapproxloc32size192 = rgapproximation(rgresultsloc32size192)
# @info("finish rg approximation")
# fiosave(rgapproxloc32size192, name="rgapproxloc32size192")

# # errorloc32size192L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc32size192-rgapproxloc32size192),4608*4*4)
# # errorloc32size192L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc32size192-rgapproxloc32size192),4608*4*4)
# # fiosave(errorloc32size192L1norm, name="errorloc32size192L1norm")
# # fiosave(errorloc32size192L2norm, name="errorloc32size192L2norm")

# @info("system size 192,local size 48")
# blockedcorrelationsloc48size192 = startingcorrelations(correlations_size192,48)
# fiosave(blockedcorrelationsloc48size192, name="blockedcorrelationsloc48size192")
# rgresultsloc48size192 = rg(blockedcorrelationsloc48size192,modeselectionstrategy)
# noofrgsteploc48size192 = length(rgresultsloc48size192)
# fiosave(noofrgsteploc48size192, name="noofrgsteploc48size192")
# noofcoremodesloc48size192 = countnoofcoremodes(rgresultsloc48size192)
# fiosave(noofcoremodesloc48size192, name="noofcoremodesloc48size192")
# rgapproxloc48size192 = rgapproximation(rgresultsloc48size192)
# @info("finish rg approximation")
# fiosave(rgapproxloc48size192, name="rgapproxloc48size192")

# # errorloc48size192L1norm = focktraceL1norm(FockMap(blockedcorrelationsloc48size192-rgapproxloc48size192),4608*4*4)
# # errorloc48size192L2norm = focktraceL2norm(FockMap(blockedcorrelationsloc48size192-rgapproxloc48size192),4608*4*4)
# # fiosave(errorloc48size192L1norm, name="errorloc48size192L1norm")
# # fiosave(errorloc48size192L2norm, name="errorloc48size192L2norm")