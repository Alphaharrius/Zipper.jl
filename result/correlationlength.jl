using LinearAlgebra
using Zipper, Plots
plotlyjs()

setmaxthreads(Threads.nthreads())

power = 3
onsitepotential = 0
nnhopping = 0.2im
systemsize=2^power*3
correlations,H = generatesystem(-onsitepotential+0im, onsitepotential+0im,-1 + 0im,-nnhopping,systemsize)

crystal=correlations|>getinspace|>getcrystal

radius = 2
generatingradius::Integer = ceil(Int, radius * 1.5) # Multiply by 1.5 to ensure all unitcell points fits.
generatinglength::Integer = generatingradius * 2
generatingcrystal::Crystal = Crystal(crystal|>getunitcell, [generatinglength, generatinglength])
crystalregion::Region = generatingcrystal|>sitepoints
crystalregion|>visualize
regionfock = getregionfock(H|>getoutspace, crystalregion)

firstgroupofsites = [getspace(crystal)*(1/3, 0),getspace(crystal)*(2/3, 0),getspace(crystal)*(2/3, 1/3)]
secondgroupofsites = [getspace(crystal)*(0, 1/3),getspace(crystal)*(0, 2/3),getspace(crystal)*(1/3, 2/3)]
sitesinchosendir = Subset(getspace(crystal)*(1/3, 0))

for r in range(0,10)
    if r%2==0
        ref = getspace(crystal)*(floor(Int,r/2),floor(Int,r/2))
        for site in firstgroupofsites
            sitesinchosendir=sitesinchosendir+Subset(ref+site)
        end
    elseif r%2==1
        ref = getspace(crystal)*(floor(Int,r/2)+1,floor(Int,r/2))
        for site in secondgroupofsites
            sitesinchosendir=sitesinchosendir+Subset(ref+site)
        end
    end
end
sitesinchosendir|>visualize
sitesinchosendir
correlationvallist = []

transform = fourier(H|>getoutspace, regionfock)
restricted = transform' * correlations * transform / (crystal|>vol)
ref = getspace(crystal)*(1/3, 0)
visualmode = getregionfock(H|>getoutspace, ref|>Subset)|>first

for offset in sitesinchosendir
    m = getregionfock(H|>getoutspace, offset|>Subset)|>first
    h = (restricted[visualmode, m]|>rep)[1]|>abs
    append!(correlationvallist, h)
end
scatter(log.(abs.((correlationvallist))))

model(t, p) = p[1] * exp.(-p[2] * t)
tdata = range(0,32)
ydata = abs.(correlationvallist)
p0 = [0.5, 0.5]
fit1 = curve_fit(model, tdata, ydata, p0)
param = fit1.param

r0 = getspace(crystal)*(1/3, 0)
m0 = getregionfock(H|>getoutspace, r0|>Subset)|>first
h0 = (restricted[visualmode, m0]|>rep)[1]|>abs
r1 = getspace(crystal)*(2/3, 0)
m1 = getregionfock(H|>getoutspace, r1|>Subset)|>first
h1 = (restricted[visualmode, m1]|>rep)[1]|>abs
r2 = getspace(crystal)*(2/3, 1/3)
m2 = getregionfock(H|>getoutspace, r2|>Subset)|>first
h2 = (restricted[visualmode, m2]|>rep)[1]|>abs
r3 = getspace(crystal)*(1, 1/3)
m3 = getregionfock(H|>getoutspace, r3|>Subset)|>first
h3 = (restricted[visualmode, m3]|>rep)[1]|>abs
r4 = getspace(crystal)*(1, 2/3)
m4 = getregionfock(H|>getoutspace, r4|>Subset)|>first
h4 = (restricted[visualmode, m4]|>rep)[1]|>abs
r5 = getspace(crystal)*(4/3, 2/3)
m5 = getregionfock(H|>getoutspace, r5|>Subset)|>first
h5 = (restricted[visualmode, m5]|>rep)[1]|>abs
r6 = getspace(crystal)*(4/3, 1)
m6 = getregionfock(H|>getoutspace, r6|>Subset)|>first
h6 = (restricted[visualmode, m6]|>rep)[1]|>abs
r7 = getspace(crystal)*(5/3, 1)
m7 = getregionfock(H|>getoutspace, r7|>Subset)|>first
h7 = (restricted[visualmode, m7]|>rep)[1]|>abs

Subset(r1,r2,r3,r4,r5,r6,r7)|>visualize
sitelist = [0,1,2,3,4,5,6,7]
correlationvallist = [h0,h1,h2,h3,h4,h5,h6,h7]
scatter(sitelist,correlationvallist)

Subset(ref,r1,r2,r3,r4)|>visualize
crystalregion|>visualize
sum(p for p in crystalregion)
crystalregion|>getcenter
centeredregion::Region = crystalregion .- (crystalregion|>getcenter)
centeredregion|>visualize

@info "Computing physical data..."
region = getsphericalregion(crystal=crystal, radius=4, metricspace=euclidean(RealSpace, 2))
region|>visualize

@info "Preparing environment..."
triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')
kspace = convert(MomentumSpace, triangular)

pa = [1/3, 2/3] ∈ triangular
pb = [2/3, 1/3] ∈ triangular
pc = (pa + pb) / 2
spatialsnappingcalibration((pa, pb, pc))

c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
c3 = c6^2

unitcell = Subset(pa, pb)
crystal = Crystal(unitcell, [192, 192])
reciprocalhashcalibration(crystal.sizes)

@info "Loading data..."
fiodir("/Users/alphaharrius/ZERData/chern192x192")
hamiltonian = fioload("hamiltonian")
correlations = fioload("correlations")

fiodir("/Users/alphaharrius/ZERData/chern192x192/rg1")
rg1correlations = fioload("couriercorrelations")
rg1block = fioload("block")
rg1zipper = fioload("wanniercourierisometry")' * rg1block

fiodir("/Users/alphaharrius/ZERData/chern192x192/rg2")
rg2correlations = fioload("couriercorrelations")
rg2block = fioload("block")
rg2zipper = fioload("wanniercourierisometry")' * rg2block

fiodir("/Users/alphaharrius/ZERData/chern192x192/rg3")
rg3correlations = fioload("couriercorrelations")
rg3block = fioload("block")
rg3zipper = fioload("wanniercourierisometry")' * rg3block

fiodir("/Users/alphaharrius/ZERData/chern192x192/rg4")
rg4correlations = fioload("couriercorrelations")
rg4block = fioload("block")
rg4zipper = fioload("wanniercourierisometry")' * rg4block

fiodir("/Users/alphaharrius/ZERData/chern192x192/rg5")
rg5correlations = fioload("couriercorrelations")
rg5block = fioload("block")
rg5zipper = fioload("wanniercourierisometry")' * rg5block

fiodir("/Users/alphaharrius/ZERData/chern192x192/rg6")
rg6correlations = fioload("couriercorrelations")
rg6block = fioload("block")
rg6zipper = fioload("wanniercourierisometry")' * rg6block

@info "Computing physical data..."
region = getsphericalregion(crystal=crystal, radius=6, metricspace=euclidean(RealSpace, 2))
ref = getspace(crystal)*(2/3, 1/3)
centerregion = Subset(r for r in region if ((ref-r)|>euclidean|>norm) <= 2.5)
regionfock = getregionfock(hamiltonian|>getoutspace, centerregion)
transform = fourier(hamiltonian|>getoutspace, regionfock)
restricted = transform' * hamiltonian * transform / (crystal|>vol)
visualmode = getregionfock(hamiltonian|>getoutspace, ref|>Subset)|>first
r1 = getspace(crystal)*(1/3, 2/3)
m1 = getregionfock(hamiltonian|>getoutspace, r1|>Subset)|>first
h1 = restricted[visualmode, m1]|>rep
r2 = getspace(crystal)*(2/3, 1 + 1/3)
m2 = getregionfock(hamiltonian|>getoutspace, r2|>Subset)|>first
h2 = restricted[visualmode, m2]|>rep
r3 = getspace(crystal)*(1 + 1/3, 1 + 2/3)
m3 = getregionfock(hamiltonian|>getoutspace, r3|>Subset)|>first
h3 = restricted[visualmode, m3]|>rep
r4 = getspace(crystal)*(1 + 2/3, 2 + 1/3)
m4 = getregionfock(hamiltonian|>getoutspace, r4|>Subset)|>first
h4 = restricted[visualmode, m4]|>rep
r5 = getspace(crystal)*(2 + 1/3, 1 + 2/3)
m5 = getregionfock(hamiltonian|>getoutspace, r5|>Subset)|>first
h5 = restricted[visualmode, m5]|>rep
r6 = getspace(crystal)*(2 + 2/3, 2 + 1/3)
m6 = getregionfock(hamiltonian|>getoutspace, r6|>Subset)|>first
h6 = restricted[visualmode, m6]|>rep
visualize(centerregion, ref|>Subset, r1|>Subset, r2|>Subset, r3|>Subset, r4|>Subset, r5|>Subset, r6|>Subset)

@info "Computing renormalization group data..."
function computergspread(rghamiltonian)
    rgcrystal = rghamiltonian|>getoutspace|>getcrystal
    rgspace = rgcrystal|>getspace
    scale = Scale([2 0; 0 2], rgspace)
    ospace = inv(scale)*rgspace|>orthospace
    region = getsphericalregion(crystal=rgcrystal, radius=12, metricspace=ospace)
    ref = rgspace*(2/3, 1/3)
    centerregion = Subset(r for r in region if (ospace*(ref-r)|>norm) < 7)
    regionfock = getregionfock(rghamiltonian|>getoutspace, centerregion)
    transform = fourier(rghamiltonian|>getoutspace, regionfock) / (rgcrystal|>vol|>sqrt)
    restricted = transform' * rghamiltonian * transform
    visualmode = getregionfock(rghamiltonian|>getoutspace, ref|>Subset)|>first
    r1 = rgspace*(1/3, 2/3)
    m1 = getregionfock(rghamiltonian|>getoutspace, r1|>Subset)|>first
    h1 = restricted[visualmode, m1]|>rep
    r2 = rgspace*(2/3, 1 + 1/3)
    m2 = getregionfock(rghamiltonian|>getoutspace, r2|>Subset)|>first
    h2 = restricted[visualmode, m2]|>rep
    r3 = rgspace*(1 + 1/3, 1 + 2/3)
    m3 = getregionfock(rghamiltonian|>getoutspace, r3|>Subset)|>first
    h3 = restricted[visualmode, m3]|>rep
    r4 = rgspace*(1 + 2/3, 2 + 1/3)
    m4 = getregionfock(rghamiltonian|>getoutspace, r4|>Subset)|>first
    h4 = restricted[visualmode, m4]|>rep
    r5 = rgspace*(2 + 1/3, 1 + 2/3)
    m5 = getregionfock(rghamiltonian|>getoutspace, r5|>Subset)|>first
    h5 = restricted[visualmode, m5]|>rep
    r6 = rgspace*(2 + 2/3, 2 + 1/3)
    m6 = getregionfock(rghamiltonian|>getoutspace, r6|>Subset)|>first
    h6 = restricted[visualmode, m6]|>rep
    return Dict(:h1=>h1[1, 1], :h2=>h2[1, 1], :h3=>h3[1, 1], :h4=>h4[1, 1], :h5=>h5[1, 1], :h6=>h6[1, 1])
end

rg1hamiltonian = rg1zipper * hamiltonian * rg1zipper'
rg1spread = computergspread(rg1hamiltonian)