using Zipper
using LinearAlgebra,Plots
plotlyjs()
fiodir("/Users/slwongag/Desktop/data/dirac/systemsize32/debugging/testingdifferentnoofdistillmodes")
gmera1firstlocalcorrelations = fioload("gmera1firstlocalcorrelations")
localcourierisometry = fioload("localcourierisometry")
localfilledisometry = fioload("localfilledisometry")
localemptyisometry = fioload("localemptyisometry")
localwanniercourier = fioload("localwanniercourier")
blockedcorrelations = fioload("blockedcorrelations")
blockedcorrelations|>getinspace|>getcrystal|>getunitcell|>visualize
globalfilledisometry = fioload("globalfilledisometry")
# globalemptyisometry = fioload("globalemptyisometry")
wanniercourierisometry = fioload("wanniercourierisometry")


(localcourierisometry'*gmera1firstlocalcorrelations*localcourierisometry)|>eigspech|>visualize
(localfilledisometry'*gmera1firstlocalcorrelations*localfilledisometry)|>eigspech|>visualize
(localemptyisometry'*gmera1firstlocalcorrelations*localemptyisometry)|>eigspech|>visualize
(localwanniercourier'*gmera1firstlocalcorrelations*localwanniercourier)|>eigspech|>visualize
filledcorrelations
filledcorrelations = globalfilledisometry'*blockedcorrelations*globalfilledisometry
filledidentity = globalfilledisometry'globalfilledisometry
courieridentity = wanniercourierisometry'wanniercourierisometry
filledidentity|>crystalspectrum|>geteigenvalues
[eval for (m,eval) in filledidentity|>crystalspectrum|>geteigenvalues]
couriercorrelations = wanniercourierisometry'*blockedcorrelations*wanniercourierisometry
couriercorrelations|>CrystalFockMap|>crystalspectrum|>visualize

filledcorrelations|>crystalspectrum|>visualize

