using Zipper
using LinearAlgebra,Plots
plotlyjs()

# Dirac
#systemsize24
fiodir("/Users/slwongag/Desktop/data/dirac/systemsize24/localsize2/")
size24local2gmeraalldiffL1normdirac = fioload("gmeraalldiffL1norm")
size24local2blockedcorrelations = fioload("blockedcorrelations")
size24local2blockedcrystalfock = size24local2blockedcorrelations|>getoutspace
size24local2blockedcrystal::Crystal = size24local2blockedcrystalfock|>getcrystal
size24local2blockedspace::RealSpace = size24local2blockedcrystal|>getspace

refrot = inv([2/3 -1/3; -1/3 2/3]')
c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
c3 = c6^2
size24local2center = [0,0] ∈ size24local2blockedspace
size24local2hexagonalregion = gethexagonalregion(rot = refrot,crystal=size24local2blockedcrystal, center=size24local2center, metricspace=size24local2blockedspace)
size24local2hexagonalregionfock = quantize(size24local2hexagonalregion,1)

# @info "Computing local correlations..."
size24local2localrestrict = fourier(size24local2blockedcrystalfock, size24local2hexagonalregionfock) / (size24local2blockedcrystal|>vol|>sqrt)
size24local2localcorrelations = size24local2localrestrict'*size24local2blockedcorrelations*size24local2localrestrict

size24local2localspectrum = size24local2localcorrelations|>eigspec
size24local2localspectrum|>visualize|>display

fiodir("/Users/slwongag/Desktop/data/dirac/systemsize24/localsize4/")
size24local4gmeraalldiffL1normdirac = fioload("gmeraalldiffL1norm")
size24local4blockedcorrelations = fioload("blockedcorrelations")
size24local4blockedcrystalfock = size24local4blockedcorrelations|>getoutspace
size24local4blockedcrystal::Crystal = size24local4blockedcrystalfock|>getcrystal
size24local4blockedspace::RealSpace = size24local4blockedcrystal|>getspace

refrot = inv([2/3 -1/3; -1/3 2/3]')
c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
c3 = c6^2
size24local4center = [0,0] ∈ size24local4blockedspace
size24local4hexagonalregion = gethexagonalregion(rot = refrot,crystal=size24local4blockedcrystal, center=size24local4center, metricspace=size24local4blockedspace)
size24local4hexagonalregionfock = quantize(size24local4hexagonalregion,1)

# @info "Computing local correlations..."
size24local4localrestrict = fourier(size24local4blockedcrystalfock, size24local4hexagonalregionfock) / (size24local4blockedcrystal|>vol|>sqrt)
size24local4localcorrelations = size24local4localrestrict'*size24local4blockedcorrelations*size24local4localrestrict

size24local4localspectrum = size24local4localcorrelations|>eigspec
size24local4localspectrum|>visualize|>display

fiodir("/Users/slwongag/Desktop/data/dirac/systemsize24/localsize6/")
size24local6gmeraalldiffL1normdirac = fioload("gmeraalldiffL1norm")
size24local6blockedcorrelations = fioload("blockedcorrelations")
size24local6blockedcrystalfock = size24local6blockedcorrelations|>getoutspace
size24local6blockedcrystal::Crystal = size24local6blockedcrystalfock|>getcrystal
size24local6blockedspace::RealSpace = size24local6blockedcrystal|>getspace

refrot = inv([2/3 -1/3; -1/3 2/3]')
c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
c3 = c6^2
size24local6center = [0,0] ∈ size24local6blockedspace
size24local6hexagonalregion = gethexagonalregion(rot = refrot,crystal=size24local6blockedcrystal, center=size24local6center, metricspace=size24local6blockedspace)
size24local6hexagonalregionfock = quantize(size24local6hexagonalregion,1)

# @info "Computing local correlations..."
size24local6localrestrict = fourier(size24local6blockedcrystalfock, size24local6hexagonalregionfock) / (size24local6blockedcrystal|>vol|>sqrt)
size24local6localcorrelations = size24local6localrestrict'*size24local6blockedcorrelations*size24local6localrestrict
size24local6localspectrum = eigspec(size24local6localcorrelations,groupingthreshold=1e-15)
size24local6localspectrum|>visualize|>display

fiodir("/Users/slwongag/Desktop/data/dirac/systemsize24/localsize8/")
size24local8gmeraalldiffL1normdirac = fioload("gmeraalldiffL1norm")
size24local8blockedcorrelations = fioload("blockedcorrelations")
size24local8blockedcrystalfock = size24local8blockedcorrelations|>getoutspace
size24local8blockedcrystal::Crystal = size24local8blockedcrystalfock|>getcrystal
size24local8blockedspace::RealSpace = size24local8blockedcrystal|>getspace

refrot = inv([2/3 -1/3; -1/3 2/3]')
c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
c3 = c6^2
size24local8center = [0,0] ∈ size24local8blockedspace
size24local8hexagonalregion = gethexagonalregion(rot = refrot,crystal=size24local8blockedcrystal, center=size24local8center, metricspace=size24local8blockedspace)
size24local8hexagonalregionfock = quantize(size24local8hexagonalregion,1)

# @info "Computing local correlations..."
size24local8localrestrict = fourier(size24local8blockedcrystalfock, size24local8hexagonalregionfock) / (size24local8blockedcrystal|>vol|>sqrt)
size24local8localcorrelations = size24local8localrestrict'*size24local8blockedcorrelations*size24local8localrestrict
size24local8localspectrum = eigspec(size24local8localcorrelations,groupingthreshold=1e-10)
size24local8localspectrum|>visualize|>display

sort([abs(eval) for (_,eval) in size24local8localspectrum|>geteigenvalues])[97]
sort([abs(eval) for (_,eval) in size24local8localspectrum|>geteigenvalues])[96]

fiodir("/Users/slwongag/Desktop/data/dirac/systemsize24/localsize12/")
size24local12gmeraalldiffL1normdirac = fioload("gmeraalldiffL1norm")
size24local12blockedcorrelations = fioload("blockedcorrelations")
size24local12blockedcrystalfock = size24local12blockedcorrelations|>getoutspace
size24local12blockedcrystal::Crystal = size24local12blockedcrystalfock|>getcrystal
size24local12blockedspace::RealSpace = size24local12blockedcrystal|>getspace

refrot = inv([2/3 -1/3; -1/3 2/3]')
c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
c3 = c6^2
size24local12center = [0,0] ∈ size24local12blockedspace
size24local12hexagonalregion = gethexagonalregion(rot = refrot,crystal=size24local12blockedcrystal, center=size24local12center, metricspace=size24local12blockedspace)
size24local12hexagonalregionfock = quantize(size24local12hexagonalregion,1)

# @info "Computing local correlations..."
size24local12localrestrict = fourier(size24local12blockedcrystalfock, size24local12hexagonalregionfock) / (size24local12blockedcrystal|>vol|>sqrt)
size24local12localcorrelations = size24local12localrestrict'*size24local12blockedcorrelations*size24local12localrestrict
size24local12localspectrum = eigspec(size24local12localcorrelations,groupingthreshold=1e-17)
size24local12localspectrum|>visualize|>display

sort([abs(eval) for (m,eval) in size24local12localspectrum|>geteigenvalues])[217]
sort([abs(eval) for (m,eval) in size24local12localspectrum|>geteigenvalues])[216]


radius_list = [2,4,6,8,12]
radius_list = [2,4,6,8]
# size24L1normdirac_list = log.([size24local2gmeraalldiffL1normdirac,size24local4gmeraalldiffL1normdirac,size24local6gmeraalldiffL1normdirac,size24local8gmeraalldiffL1normdirac,size24local12gmeraalldiffL1normdirac])
size24L1normdirac_list = log.([size24local2gmeraalldiffL1normdirac,size24local4gmeraalldiffL1normdirac,size24local6gmeraalldiffL1normdirac,size24local8gmeraalldiffL1normdirac])

#systemsize48
fiodir("/Users/slwongag/Desktop/data/dirac/systemsize48/localsize2/")
size48local2gmeraalldiffL1normdirac = fioload("gmeraalldiffL1norm")

# Trivial 
#onsite0.05
fiodir("/Users/slwongag/Desktop/data/trivial/systemsize24/onsitepotential0.05-localsize2")
size24local2gmeraalldiffL1normtrivialonsite005 = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/trivial/systemsize24/onsitepotential0.05-localsize4")
size24local4gmeraalldiffL1normtrivialonsite005 = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/trivial/systemsize24/onsitepotential0.05-localsize6")
size24local6gmeraalldiffL1normtrivialonsite005 = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/trivial/systemsize24/onsitepotential0.05-localsize8")
size24local8gmeraalldiffL1normtrivialonsite005 = fioload("gmeraalldiffL1norm")

# fiodir("/Users/slwongag/Desktop/data/trivial/systemsize24/onsitepotential0.05-localsize12")
# size24local12gmeraalldiffL1normtrivialonsite005 = fioload("gmeraalldiffL1norm")

# size24L1normtrivialonsite005_list = log.([size24local2gmeraalldiffL1normtrivialonsite005,size24local4gmeraalldiffL1normtrivialonsite005,size24local6gmeraalldiffL1normtrivialonsite005,size24gmeraalldiffL1normtriviallocal8onsite005,size24gmeraalldiffL1normtriviallocal12onsite005])
size24L1normtrivialonsite005_list = log.([size24local2gmeraalldiffL1normtrivialonsite005,size24local4gmeraalldiffL1normtrivialonsite005,size24local6gmeraalldiffL1normtrivialonsite005,size24local8gmeraalldiffL1normtrivialonsite005])

#onsite0.1
fiodir("/Users/slwongag/Desktop/data/trivial/systemsize24/onsitepotential0.1-localsize2")
size24local2gmeraalldiffL1normtrivialonsite01 = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/trivial/systemsize24/onsitepotential0.1-localsize4")
size24local4gmeraalldiffL1normtrivialonsite01 = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/trivial/systemsize24/onsitepotential0.1-localsize6")
size24local6gmeraalldiffL1normtrivialonsite01 = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/trivial/systemsize24/onsitepotential0.1-localsize8")
size24local8gmeraalldiffL1normtrivialonsite01 = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/trivial/systemsize24/onsitepotential0.1-localsize12")
size24local12gmeraalldiffL1normtrivialonsite01 = fioload("gmeraalldiffL1norm")

# size24L1normtrivialonsite01_list = log.([size24local2gmeraalldiffL1normtrivialonsite01,size24local4gmeraalldiffL1normtrivialonsite01,size24local6gmeraalldiffL1normtrivialonsite01,size24local8gmeraalldiffL1normtrivialonsite01,size24local12gmeraalldiffL1normtrivialonsite01])
size24L1normtrivialonsite01_list = log.([size24local2gmeraalldiffL1normtrivialonsite01,size24local4gmeraalldiffL1normtrivialonsite01,size24local6gmeraalldiffL1normtrivialonsite01,size24local8gmeraalldiffL1normtrivialonsite01])

#onsite0.15
fiodir("/Users/slwongag/Desktop/data/trivial/systemsize24/onsitepotential0.15-localsize2")
size24local2gmeraalldiffL1normtrivialonsite015 = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/trivial/systemsize24/onsitepotential0.15-localsize4")
size24local4gmeraalldiffL1normtrivialonsite015 = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/trivial/systemsize24/onsitepotential0.15-localsize6")
size24local6gmeraalldiffL1normtrivialonsite015 = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/trivial/systemsize24/onsitepotential0.15-localsize8")
size24local8gmeraalldiffL1normtrivialonsite015 = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/trivial/systemsize24/onsitepotential0.15-localsize12")
size24local12gmeraalldiffL1normtrivialonsite015 = fioload("gmeraalldiffL1norm")

# size24L1normtrivialonsite015_list = log.([size24local2gmeraalldiffL1normtrivialonsite015,size24local4gmeraalldiffL1normtrivialonsite015,size24local6gmeraalldiffL1normtrivialonsite015,size24local8gmeraalldiffL1normtrivialonsite015,size24local12gmeraalldiffL1normtrivialonsite015])
size24L1normtrivialonsite015_list = log.([size24local2gmeraalldiffL1normtrivialonsite015,size24local4gmeraalldiffL1normtrivialonsite015,size24local6gmeraalldiffL1normtrivialonsite015,size24local8gmeraalldiffL1normtrivialonsite015])

#onsite0.2
fiodir("/Users/slwongag/Desktop/data/trivial/systemsize24/onsitepotential0.2-localsize2")
size24local2gmeraalldiffL1normtrivialonsite02 = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/trivial/systemsize24/onsitepotential0.2-localsize4")
size24local4gmeraalldiffL1normtrivialonsite02 = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/trivial/systemsize24/onsitepotential0.2-localsize6")
size24local6gmeraalldiffL1normtrivialonsite02 = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/trivial/systemsize24/onsitepotential0.2-localsize8")
size24local8gmeraalldiffL1normtrivialonsite02 = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/trivial/systemsize24/onsitepotential0.2-localsize12")
size24local12gmeraalldiffL1normtrivialonsite02 = fioload("gmeraalldiffL1norm")

# size24L1normtrivialonsite02_list = log.([size24local2gmeraalldiffL1normtrivialonsite02,size24local4gmeraalldiffL1normtrivialonsite02,size24local6gmeraalldiffL1normtrivialonsite02,size24local8gmeraalldiffL1normtrivialonsite02,size24local12gmeraalldiffL1normtrivialonsite02])
size24L1normtrivialonsite02_list = log.([size24local2gmeraalldiffL1normtrivialonsite02,size24local4gmeraalldiffL1normtrivialonsite02,size24local6gmeraalldiffL1normtrivialonsite02,size24local8gmeraalldiffL1normtrivialonsite02])

#result of chern
#0.1im
fiodir("/Users/slwongag/Desktop/data/chern/systemsize24/nnhopping0.0 + 0.1im-localsize2")
size24local2gmeraalldiffL1normchernnn01im = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/chern/systemsize24/nnhopping0.0 + 0.1im-localsize4")
size24local4gmeraalldiffL1normchernnn01im = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/chern/systemsize24/nnhopping0.0 + 0.1im-localsize6")
size24local6gmeraalldiffL1normchernnn01im = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/chern/systemsize24/nnhopping0.0 + 0.1im-localsize8")
size24local8gmeraalldiffL1normchernnn01im = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/chern/systemsize24/nnhopping0.0 + 0.1im-localsize12")
size24local12gmeraalldiffL1normchernnn01im = fioload("gmeraalldiffL1norm")

#0.2im
fiodir("/Users/slwongag/Desktop/data/chern/systemsize24/nnhopping0.0 + 0.2im-localsize2")
size24local2gmeraalldiffL1normchernnn02im = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/chern/systemsize24/nnhopping0.0 + 0.2im-localsize4")
size24local4gmeraalldiffL1normchernnn02im = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/chern/systemsize24/nnhopping0.0 + 0.2im-localsize6")
size24local6gmeraalldiffL1normchernnn02im = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/chern/systemsize24/nnhopping0.0 + 0.2im-localsize8")
size24local8gmeraalldiffL1normchernnn02im = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/chern/systemsize24/nnhopping0.0 + 0.2im-localsize12")
size24local12gmeraalldiffL1normchernnn02im = fioload("gmeraalldiffL1norm")

#0.3im
fiodir("/Users/slwongag/Desktop/data/chern/systemsize24/nnhopping0.0 + 0.3im-localsize2")
size24local2gmeraalldiffL1normchernnn03im = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/chern/systemsize24/nnhopping0.0 + 0.3im-localsize4")
size24local4gmeraalldiffL1normchernnn03im = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/chern/systemsize24/nnhopping0.0 + 0.3im-localsize6")
size24local6gmeraalldiffL1normchernnn03im = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/chern/systemsize24/nnhopping0.0 + 0.3im-localsize8")
size24local8gmeraalldiffL1normchernnn03im = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/chern/systemsize24/nnhopping0.0 + 0.3im-localsize12")
size24local12gmeraalldiffL1normchernnn03im = fioload("gmeraalldiffL1norm")

#0.4im
fiodir("/Users/slwongag/Desktop/data/chern/systemsize24/nnhopping0.0 + 0.4im-localsize2")
size24local2gmeraalldiffL1normchernnn04im = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/chern/systemsize24/nnhopping0.0 + 0.4im-localsize4")
size24local4gmeraalldiffL1normchernnn04im = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/chern/systemsize24/nnhopping0.0 + 0.4im-localsize6")
size24local6gmeraalldiffL1normchernnn04im = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/chern/systemsize24/nnhopping0.0 + 0.4im-localsize8")
size24local8gmeraalldiffL1normchernnn04im = fioload("gmeraalldiffL1norm")

fiodir("/Users/slwongag/Desktop/data/chern/systemsize24/nnhopping0.0 + 0.4im-localsize12")
size24local12gmeraalldiffL1normchernnn04im = fioload("gmeraalldiffL1norm")


size24L1normchernnn01im_list = log.([size24local2gmeraalldiffL1normchernnn01im,size24local4gmeraalldiffL1normchernnn01im,size24local6gmeraalldiffL1normchernnn01im,size24local8gmeraalldiffL1normchernnn01im])
size24L1normchernnn02im_list = log.([size24local2gmeraalldiffL1normchernnn02im,size24local4gmeraalldiffL1normchernnn02im,size24local6gmeraalldiffL1normchernnn02im,size24local8gmeraalldiffL1normchernnn02im])
size24L1normchernnn03im_list = log.([size24local2gmeraalldiffL1normchernnn03im,size24local4gmeraalldiffL1normchernnn03im,size24local6gmeraalldiffL1normchernnn03im,size24local8gmeraalldiffL1normchernnn03im])
size24L1normchernnn04im_list = log.([size24local2gmeraalldiffL1normchernnn04im,size24local4gmeraalldiffL1normchernnn04im,size24local6gmeraalldiffL1normchernnn04im,size24local8gmeraalldiffL1normchernnn04im])

scatter(radius_list,size24L1normdirac_list,title="log(error) vs radius",label="Dirac: t/V=1/0=inf")
scatter!(radius_list,size24L1normtrivialonsite005_list,label="Trivial: t/V=1/0.05=20")
scatter!(radius_list,size24L1normtrivialonsite01_list,label="Trivial: t/V=1/0.1=10")
scatter!(radius_list,size24L1normtrivialonsite015_list,label="Trivial: t/V=1/0.15=6.67")
scatter!(radius_list,size24L1normtrivialonsite02_list,label="Trivial: t/V=1/0.2=5")
xlabel!("radius")
ylabel!("log of vectorize L1 norm of diff")

scatter!(radius_list,size24L1normchernnn01im_list,label="Chern: t=1,tt=0.1im")
scatter!(radius_list,size24L1normchernnn02im_list,label="Chern: t=1,tt=0.2im")
scatter!(radius_list,size24L1normchernnn03im_list,label="Chern: t=1,tt=0.3im")
scatter!(radius_list,size24L1normchernnn04im_list,label="Chern: t=1,tt=0.4im")
xlabel!("radius")
ylabel!("log of vectorize L1 norm of diff")


#Dirac across RG
fiodir("/Users/slwongag/Desktop/data/dirac/systemsize24/localsize2/")
size24local2gmera1terminatedalldiffL1norm = fioload("gmera1terminatedalldiffL1norm")
size24local2gmeraalldiffL1normdirac = fioload("gmeraalldiffL1norm")