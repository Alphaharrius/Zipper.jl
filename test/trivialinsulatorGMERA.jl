using Plotly, SmithNormalForm, LinearAlgebra, OrderedCollections, SparseArrays, Combinatorics
using Zipper

setmaxthreads(Threads.nthreads())

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

H, correlations = H_trivial(96,ComplexF64(-1),-0.2,-0.8)

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

blockedcorrelations = startingcorrelations(correlations,8)

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
    transistionmapdict = Dict(:orig=>transistionmap)

    firstcenterlist = [[1/2,0] ∈ rgblockedspace,[-1/2,0] ∈ rgblockedspace]
    secondcenterlist = [[0,1/2] ∈ rgblockedspace,[0,-1/2] ∈ rgblockedspace]
    thirdcenterlist = [[1/2,1/2] ∈ rgblockedspace,[-1/2,-1/2] ∈ rgblockedspace, [-1/2,1/2] ∈ rgblockedspace,[1/2,-1/2] ∈ rgblockedspace]
    finalcenterlist = [[0,0] ∈ rgblockedspace]

    @info ("1st gmera step...")
    gmera1 = @time gmerastep(rgblockedcorrelations,rgblockedcorrelations,firstcenterlist, modeselection1stbycountthenbythreshold(1,0.005))
    # gmera1 = @time gmerastep1(rgblockedcorrelations,firstcenterlist)
    @info ("1st gmera approximation to correlations...")
    gmera1approx = transistionmap*gmera1[:emptyisometry]*gmera1[:emptyisometry]'*transistionmap'
    approxdict = Dict(:first=>gmera1approx)
    filledisometrydict = Dict(:first=>gmera1[:filledisometry])
    emptyisometrydict = Dict(:first=>gmera1[:emptyisometry])

    if haskey(gmera1,:correlations)
        transistionmap = transistionmap*gmera1[:courierisometry]
        correlationsdict = Dict(:first=>gmera1[:correlations])
        transistionmapdict[:first] = transistionmap
        @info ("2nd gmera step...")
        gmera2 = @time gmerastep(rgblockedcorrelations,gmera1[:correlations],secondcenterlist,modeselection1stbycountthenbythreshold(1,0.005))
        @info ("2nd gmera approximation to correlations...")
        gmera2approx = transistionmap*gmera2[:emptyisometry]*gmera2[:emptyisometry]'*transistionmap'
        approxdict[:second] = gmera2approx
        filledisometrydict[:second] = gmera2[:filledisometry]
        emptyisometrydict[:second] = gmera2[:emptyisometry]
        if haskey(gmera2,:correlations)
            transistionmap = transistionmap*gmera2[:courierisometry]
            transistionmapdict[:second] = transistionmap
            correlationsdict[:second] = gmera2[:correlations]
            @info ("3rd gmera step...")
            gmera3 = @time gmerastep(rgblockedcorrelations,gmera2[:correlations],thirdcenterlist,modeselection1stbycountthenbythreshold(1,0.005))
            @info ("3rd gmera approximation to correlations...")
            gmera3approx = transistionmap*gmera3[:emptyisometry]*gmera3[:emptyisometry]'*transistionmap'
            approxdict[:third] = gmera3approx
            filledisometrydict[:third] = gmera3[:filledisometry]
            emptyisometrydict[:third] = gmera3[:emptyisometry]
            if haskey(gmera3,:correlations)
                transistionmap = transistionmap*gmera3[:courierisometry]
                transistionmapdict[:third] = transistionmap
                correlationsdict[:third] = gmera3[:correlations]
                @info ("final gmera step...")
                gmera4 = @time gmerastep(rgblockedcorrelations,gmera3[:correlations],finalcenterlist,modeselection1stbycountthenbythreshold(1,0.005))
                @info ("4th gmera approximation to correlations...")
                gmera4approx = transistionmap*gmera4[:emptyisometry]*gmera4[:emptyisometry]'*transistionmap'
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
                        :approxresults => approxdict,
                        :transistionmaps => transistionmapdict,
                        :correlations => gmera4[:correlations])
                else
                    @info("terminate at 4th gmera step")
                    return Dict(
                        :rgblockedmap => rgblocker,
                        :emptyisometryresults => emptyisometrydict,
                        :filledisometryresults => filledisometrydict,
                        :correlationsresults => correlationsdict,
                        :approxresults => approxdict,
                        :transistionmaps => transistionmapdict)
                end
            else
                @info("terminate at 3rd gmera step")
                return Dict(
                    :rgblockedmap => rgblocker,
                    :emptyisometryresults => emptyisometrydict,
                    :filledisometryresults => filledisometrydict,
                    :correlationsresults => correlationsdict,
                    :approxresults => approxdict,
                    :transistionmaps => transistionmapdict)
            end
        else
            @info("terminate at 2nd gmera step")
            return Dict(
                :rgblockedmap => rgblocker,
                :emptyisometryresults => emptyisometrydict,
                :filledisometryresults => filledisometrydict,
                :correlationsresults => correlationsdict,
                :approxresults => approxdict,
                :transistionmaps => transistionmapdict
                )
        end
    else
        @info("terminate at 1st gmera step")
        return Dict(
            :rgblockedmap => rgblocker,
            :emptyisometryresults => emptyisometrydict,
            :filledisometryresults => filledisometrydict,
            :approxresults => approxdict,
            :transistionmaps => transistionmapdict)
    end
end

rg1 = gmera(blockedcorrelations,idmap(blockedcorrelations|>getinspace))
rg2 = gmera(rg1[:correlations],rg1[:transistionmaps][:forth])
# rg3 = gmera(rg2[:correlations],rg2[:transistionmap])

rg1[:approxresults][:first]+rg1[:approxresults][:second]
rg2[:approxresults]

sum([pair[2] for rgresult in [rg1[:approxresults],rg2[:approxresults]] for pair in rgresult ])

a = [1,2,3]
b = [4,5,6]
[(i,j) for i in a for j in b]

FockMap(blockedcorrelations*blockedcorrelations') |> eigspech|>geteigenvalues

