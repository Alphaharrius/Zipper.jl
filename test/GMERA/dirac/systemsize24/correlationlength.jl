using LinearAlgebra
using Zipper, Plots
plotlyjs()

setmaxthreads(Threads.nthreads())

power = 3
onsitepotential = 0
nnhopping = 0
systemsize=2^power*3
correlations,H = generatesystem(-onsitepotential+0im, onsitepotential+0im,-1 + 0im,-nnhopping,systemsize)

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