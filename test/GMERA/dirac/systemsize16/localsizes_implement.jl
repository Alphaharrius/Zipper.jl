using Pkg
Pkg.resolve()
Pkg.Registry.update()
Pkg.instantiate()
Pkg.activate("../../../../")
Pkg.resolve()

using Zipper
using LinearAlgebra,Plots

setmaxthreads(setmaxthreads(Threads.nthreads()))
usecrystaldensemap()
# plotlyjs()

power = 4
onsitepotential = 0
nnhopping = 0
systemsize=2^power
correlations,H = generatesystem(-onsitepotential+0im, onsitepotential+0im,-1 + 0im,-nnhopping,systemsize)
crystalfock = correlations|>getoutspace

scalingswifnoofrgsteps = [(2^i,power-i+1) for i in 1:power-1]

for (scaling,noofrgsteps) in scalingswifnoofrgsteps
    fiodir("../../../../../../../storage/data/slwongag/GMERA/dirac/systemsize$systemsize/localsize$scaling")
    @info("first gmera step")
    rgH,rgcorrelations,couriercomposemap,gmeraapproximatecorrelationsofar,blockedcorrelations,noofflavourpermodeforlaterrg = firstgmerastep(correlations,H,scaling)
    refcorrelations = rgcorrelations
    refH = rgH
    refgmeraapproximatecorrelation = gmeraapproximatecorrelationsofar
    refcouriercomposemap = couriercomposemap
    for rgstep in 2:noofrgsteps
        if rgstep == noofrgsteps
            @info("final gmera step")
            finalgmerastep(refcorrelations,refH,refcouriercomposemap,refgmeraapproximatecorrelation,blockedcorrelations,rgstep,systemsize,noofflavourpermodeforlaterrg)
        else
            @info("intermediate gmera step")
            rgH,rgcorrelations,couriercomposemap,gmeraapproximatecorrelationsofar = intermediategmerastep(refcorrelations,refH,refcouriercomposemap,refgmeraapproximatecorrelation,rgstep,noofflavourpermodeforlaterrg)
        end
        refcorrelations = rgcorrelations
        refH = rgH
        refgmeraapproximatecorrelation = gmeraapproximatecorrelationsofar
        refcouriercomposemap = couriercomposemap
    end
end

