using Zipper
using LinearAlgebra,Plots
plotlyjs()

usecrystaldensemap()

function gmeraterminatedstep(correlations,H,noofflavourpermode::Number)
    crystalfock = correlations|>getoutspace
    crystal::Crystal = crystalfock|>getcrystal
    space::RealSpace = crystal|>getspace
    firstcenter = [0,0] ∈ space
    refrot = inv([2/3 -1/3; -1/3 2/3]')
    hexagonalregion = gethexagonalregion(rot = refrot,crystal=crystal, center=firstcenter, metricspace=space)
    hexagonalregionfock = quantize(hexagonalregion,noofflavourpermode)

    nooffrozenmodes = length(hexagonalregionfock|>orderedmodes)
    nooffilledmodes = div(nooffrozenmodes,2)
    noofemptymodes = nooffilledmodes
    @info("no of filled modes = ",nooffilledmodes)
    @info("no of empty modes = ",noofemptymodes)

    # @info "Computing local correlations..."
    localrestrict = fourier(crystalfock, hexagonalregionfock) / (crystal|>vol|>sqrt)
    localcorrelations = localrestrict'*correlations*localrestrict
    localspectrum = localcorrelations|>eigspech
    localspectrum|>visualize|>display

    localstates = getregionstates(localcorrelations=localcorrelations, grouping=[nooffilledmodes, noofemptymodes])
    localfilledisometry = localstates[1]|>FockMap
    localemptyisometry = localstates[2]|>FockMap

    filledunitcellfock = localfilledisometry|>getinspace|>unitcellfock
    filledunitcellmodes = Subset(m for m in filledunitcellfock)
    emptyunitcellfock = localemptyisometry|>getinspace|>unitcellfock
    emptyunitcellmodes = Subset(m for m in emptyunitcellfock)

        
    newfilledcrystalfock = getcrystalfock(localfilledisometry|>getinspace|>unitcellfock, Crystal(filledunitcellmodes|>basisoffsetofmodes,crystal.sizes))
    newemptycrystalfock = getcrystalfock(localemptyisometry|>getinspace|>unitcellfock, Crystal(emptyunitcellmodes|>basisoffsetofmodes,crystal.sizes))

    newlocalfilledrestrict = fourier(newfilledcrystalfock, localfilledisometry|>getinspace) * (crystal|>vol|>sqrt)
    newlocalemptyrestrict = fourier(newemptycrystalfock, localemptyisometry|>getinspace) * (crystal|>vol|>sqrt)

    globalfilledisometry = broadcast(*,(localrestrict*localfilledisometry), newlocalfilledrestrict')
    globalemptyisometry = broadcast(*,(localrestrict*localemptyisometry), newlocalemptyrestrict')

    filledcorrelations = globalfilledisometry'*correlations*globalfilledisometry
    emptycorrelations = globalemptyisometry'*correlations*globalemptyisometry

    filledH = globalfilledisometry' * H * globalfilledisometry
    emptyH = globalemptyisometry' * H * globalemptyisometry

        return Dict(
        :localcorrelations=>localcorrelations,
        :filledcorrelations=>filledcorrelations,
        :emptycorrelations=>emptycorrelations,
        :globalfilledisometry=>globalfilledisometry,
        :globalemptyisometry=>globalemptyisometry,
        :filledH=>filledH,
        :emptyH=>emptyH)
end

# Dirac
#systemsize24
fiodir("/Users/slwongag/Desktop/data/dirac/systemsize24/localsize2/")
systemsize = 24
blockedcorrelations = fioload("blockedcorrelations")
blockedH = fioload("blockedH")
gmerafirststep(blockedcorrelations,blockedH,0,1)

size24local2gmera1distilled = distillation(fioload("gmera1thirdcouriercorrelations")|>crystalspectrum, :filled=>(v->v < 0.5), :empty=>(v->v > 0.5))
size24local2gmera1emptybands = size24local2gmera1distilled[:empty]
size24local2gmera1emptyisometry = crystalisometry(size24local2gmera1emptybands)|>CrystalDenseMap
size24local2gmera1firstglobalcourierisometry = fioload("gmera1firstglobalcourierisometry")
size24local2gmera1secondglobalcourierisometry = fioload("gmera1secondglobalcourierisometry")
size24local2gmera1thirdglobalcourierisometry = fioload("gmera1thirdglobalcourierisometry")
size24local2gmera1couriercomposemap = size24local2gmera1firstglobalcourierisometry*size24local2gmera1secondglobalcourierisometry*size24local2gmera1thirdglobalcourierisometry
size24local2gmera1terminatedcouriercorrelation = (size24local2gmera1couriercomposemap*size24local2gmera1emptyisometry)*(size24local2gmera1couriercomposemap*size24local2gmera1emptyisometry)'
size24local2gmera1approximatecorrelation = fioload("gmera1approximatecorrelation")
size24local2gmera1terminatedapproximatecorrelation = size24local2gmera1approximatecorrelation+size24local2gmera1terminatedcouriercorrelation

gmera1terminatedalldiff = blockedcorrelations - size24local2gmera1terminatedapproximatecorrelation
# fiosave(gmera1terminatedalldiff, name="gmeraalldiff")

@info("calculating L1 norm of the difference between the blocked correlations and the terminated approximate correlation")
gmera1terminatedalldiffL1norm = focktraceL1norm(gmera1terminatedalldiff,systemsize^2*6)
# fiosave(gmera1terminatedalldiffL1norm, name="gmera1terminatedalldiffL1norm")


size24local2gmera2distilled = distillation(fioload("gmera2thirdcouriercorrelations")|>crystalspectrum, :filled=>(v->v < 0.5), :empty=>(v->v > 0.5))
size24local2gmera2emptybands = size24local2gmera2distilled[:empty]
size24local2gmera2emptyisometry = crystalisometry(size24local2gmera2emptybands)|>CrystalDenseMap
size24local2gmera2couriercomposemap = fioload("gmera2couriercomposemapgmera")
size24local2gmera1terminatedcouriercorrelation = (size24local2gmera2couriercomposemap*size24local2gmera2emptyisometry)*(size24local2gmera2couriercomposemap*size24local2gmera2emptyisometry)'
size24local2gmera2approximatecorrelation = fioload("gmera2approximatecorrelation")
size24local2gmera2terminatedapproximatecorrelation = size24local2gmera2approximatecorrelation+size24local2gmera1terminatedcouriercorrelation

gmera2terminatedalldiff = blockedcorrelations - size24local2gmera2terminatedapproximatecorrelation
fiosave(gmera2terminatedalldiff, name="gmera2terminatedgmeraalldiff")

@info("calculating L1 norm of the difference between the blocked correlations and the terminated approximate correlation")
gmera2terminatedalldiffL1norm = focktraceL1norm(gmera2terminatedalldiff,systemsize^2*6)
fiosave(gmera2terminatedalldiffL1norm, name="gmera2terminatedalldiffL1norm")