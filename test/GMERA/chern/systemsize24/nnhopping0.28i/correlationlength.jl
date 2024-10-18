using LinearAlgebra
using Zipper, Plots
plotlyjs()

setmaxthreads(Threads.nthreads())

power = 3
onsitepotential = 0
nnhopping = 0.28im
systemsize=2^power*3
correlations,H = generatesystem(-onsitepotential+0im, onsitepotential+0im,-1 + 0im,-nnhopping,systemsize)
crystal=correlations|>getinspace|>getcrystal

radius = 4
generatingradius::Integer = ceil(Int, radius * 1.5) # Multiply by 1.5 to ensure all unitcell points fits.
generatinglength::Integer = generatingradius * 2
generatingcrystal::Crystal = Crystal(crystal|>getunitcell, [generatinglength, generatinglength])
crystalregion::Region = generatingcrystal|>sitepoints
crystalregion|>visualize
regionfock = getregionfock(H|>getoutspace, crystalregion)
transform = fourier(H|>getoutspace, regionfock)
restricted = transform' * correlations * transform / (crystal|>vol)

# Looking at A sites
# firstgroupofsitesA = [getspace(crystal)*(1/3, 0),getspace(crystal)*(2/3, 1/3)]
# secondgroupofsitesA = [getspace(crystal)*(0, 2/3)]
firstgroupofsitesA = [getspace(crystal)*(1/3, 0),getspace(crystal)*(2/3, 0),getspace(crystal)*(2/3, 1/3)]
secondgroupofsitesA = [getspace(crystal)*(0, 1/3),getspace(crystal)*(0, 2/3),getspace(crystal)*(1/3, 2/3)]
sitesinchosendirA = Subset(getspace(crystal)*(1/3, 0))

for r in range(0,15)
    if r%2==0
        ref = getspace(crystal)*(floor(Int,r/2),floor(Int,r/2))
        for site in firstgroupofsitesA
            sitesinchosendirA=sitesinchosendirA+Subset(ref+site)
        end
    elseif r%2==1
        ref = getspace(crystal)*(floor(Int,r/2)+1,floor(Int,r/2))
        for site in secondgroupofsitesA
            sitesinchosendirA=sitesinchosendirA+Subset(ref+site)
        end
    end
end
sitesinchosendirA|>visualize
correlationvallistA = []
refA = getspace(crystal)*(1/3, 0)
visualmodeA = getregionfock(H|>getoutspace, refA|>Subset)|>first

for offset in sitesinchosendirA
    m = getregionfock(H|>getoutspace, offset|>Subset)|>first
    h = (restricted[visualmodeA, m]|>rep)[1]|>abs
    append!(correlationvallistA, h)
end
scatter(log.(abs.((correlationvallistA[2:end]))))


distancelistA = [r for r in range(1,47)]
model(t, p) = p[1] .+ p[2]*t 
p0 = [1.0, 1.0]
fit = curve_fit(model, distancelistA, log.(abs.((correlationvallistA[2:end]))), p0)
param = fit.param
correlationlength = -1/param[2]

fiodir("/Users/slwongag/Desktop/data/chern/systemsize24")
fiosave(correlationlength, name="correlationlength028im")

# # Looking at B sites
# firstgroupofsitesB = [getspace(crystal)*(2/3, 0)]
# secondgroupofsitesB = [getspace(crystal)*(0, 1/3),getspace(crystal)*(1/3, 2/3)]
# sitesinchosendirB = Subset(getspace(crystal)*(2/3, 0))

# for r in range(0,10)
#     if r%2==0
#         ref = getspace(crystal)*(floor(Int,r/2),floor(Int,r/2))
#         for site in firstgroupofsitesB
#             sitesinchosendirB=sitesinchosendirB+Subset(ref+site)
#         end
#     elseif r%2==1
#         ref = getspace(crystal)*(floor(Int,r/2)+1,floor(Int,r/2))
#         for site in secondgroupofsitesB
#             sitesinchosendirB=sitesinchosendirB+Subset(ref+site)
#         end
#     end
# end
# sitesinchosendirB|>visualize
# correlationvallistB = []
# refB = getspace(crystal)*(2/3, 0)
# visualmodeB = getregionfock(H|>getoutspace, refB|>Subset)|>first

# for offset in sitesinchosendirB
#     m = getregionfock(H|>getoutspace, offset|>Subset)|>first
#     h = (restricted[visualmodeB, m]|>rep)[1]|>abs
#     append!(correlationvallistB, h)
# end
# scatter(log.(abs.((correlationvallistB[2:end]))))