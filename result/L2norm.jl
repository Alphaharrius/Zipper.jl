using Zipper
using LinearAlgebra,Plots
plotlyjs()

usecrystaldensemap()

function getchunkinfo(outcrystal::Crystal, incrystal::Crystal)
    chunkcount::Integer = maximum((size(outcrystal)..., size(incrystal)...))
    datalength::Integer = vol(outcrystal) * vol(incrystal)
    chunksize::Integer = datalength/chunkcount|>round
    return datalength, chunkcount, chunksize
end

function getdenseindex(chunkno::Integer, chunkid::Integer, chunksize::Integer)
    return (chunkno-1)*chunksize + chunkid
end

function getdensemomentums(outcrystal::Crystal, incrystal::Crystal, denseindex::Integer)
    outindex = ceil(denseindex/vol(incrystal))|>Integer
    inindex = (denseindex-1)%vol(incrystal)+1
    return outcrystal[outindex], incrystal[inindex]
end

preparedense(chunkcount::Integer, chunksize::Integer) = [
    spzeros(SparseFockMap, chunksize) for _ in 1:chunkcount]
preparedenselocks(chunkcount::Integer) = [ReentrantLock() for _ in 1:chunkcount]

# Dirac
#systemsize24
fiodir("/Users/slwongag/Desktop/data/dirac/systemsize48/localsize2/")
gmeraalldiff = fioload("gmeraalldiff")
gmeraalldiff 

loutspace = gmeraalldiff|>getoutspace
linspace = gmeraalldiff|>getinspace
routspace = gmeraalldiff|>getoutspace
rinspace = gmeraalldiff|>getinspace
@assert hassamespan(linspace, routspace)
rindextable = outspacesubmaps(gmeraalldiff)
_, chunkcount, chunksize = getchunkinfo(loutspace|>getcrystal, rinspace|>getcrystal)
for (chunkno, chunkid) in gmeraalldiff.nonzeroids
    denseindex = getdenseindex(chunkno, chunkid, gmeraalldiff.chunksize)
    outkl, inkl = getdensemomentums(loutspace|>getcrystal, linspace|>getcrystal, denseindex)
    # rblocks = rindextable[inkl]
    lblock = gmeraalldiff.data[chunkno][chunkid]
end

data = preparedense(chunkcount, chunksize)
    locks = preparedenselocks(chunkcount)
gmeraalldiff.nonzeroids
rindextable = outspacesubmaps(right)


size24local2gmera1firstdiffL2norm = fioload("gmera1firstalldiffL2norm")["re"]
size24local2gmera1seconddiffL2norm = fioload("gmera1secondalldiffL2norm")["re"]
size24local2gmera1thirddiffL2norm = fioload("gmera1thirdalldiffL2norm")["re"]

size24local2gmera2firstdiffL2norm = fioload("gmera2firstalldiffL2norm")["re"]
size24local2gmera2seconddiffL2norm = fioload("gmera2secondalldiffL2norm")["re"]
size24local2gmera2thirddiffL2norm = fioload("gmera2thirdalldiffL2norm")["re"]

size24local2gmera3firstdiffL2norm = fioload("gmera3firstalldiffL2norm")["re"]
size24local2gmera3seconddiffL2norm = fioload("gmera3secondalldiffL2norm")["re"]
size24local2gmera3thirddiffL2norm = fioload("gmera3thirdalldiffL2norm")["re"]

size24local2gmeraalldiffL2normdirac = fioload("gmeraalldiffL2norm")["re"]

size24local2gmeradirac_acrossRGlist = [size24local2gmera1thirddiffL1norm,size24local2gmera2thirddiffL1norm,size24local2gmera3thirddiffL1norm,size24local2gmeraalldiffL1normdirac]
scatter(size24local2gmeradirac_acrossRGlist)
