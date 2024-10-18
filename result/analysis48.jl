using Zipper
using LinearAlgebra,Plots
plotlyjs()
usecrystaldensemap()

# Dirac
#systemsize48
fiodir("/Users/slwongag/Desktop/data/dirac/systemsize48/localsize2/L2norm/")
size48local2gmera1thirddiffL2norm = fioload("gmera1thirdalldiffL2norm")["re"]
size48local2gmera2thirddiffL2norm = fioload("gmera2thirdalldiffL2norm")["re"]
size48local2gmera3thirddiffL2norm = fioload("gmera3thirdalldiffL2norm")["re"]
size48local2gmera4thirddiffL2norm = fioload("gmera4thirdalldiffL2norm")["re"]
size48local2gmeraalldiffL2normdirac = fioload("gmeraalldiffL2norm")["re"]

fiodir("/Users/slwongag/Desktop/data/dirac/systemsize48/localsize4/L2norm/")
size48local4gmera1thirddiffL2norm = fioload("gmera1thirdalldiffL2norm")["re"]
size48local4gmera2thirddiffL2norm = fioload("gmera2thirdalldiffL2norm")["re"]
size48local4gmera3thirddiffL2norm = fioload("gmera3thirdalldiffL2norm")["re"]
size48local4gmeraalldiffL2normdirac = fioload("gmeraalldiffL2norm")["re"]

fiodir("/Users/slwongag/Desktop/data/dirac/systemsize48/localsize6/L2norm/")
size48local6gmera1thirddiffL2norm = fioload("gmera1thirdalldiffL2norm")["re"]
size48local6gmera2thirddiffL2norm = fioload("gmera2thirdalldiffL2norm")["re"]
size48local6gmera3thirddiffL2norm = fioload("gmera3thirdalldiffL2norm")["re"]
size48local6gmeraalldiffL2normdirac = fioload("gmeraalldiffL2norm")["re"]

fiodir("/Users/slwongag/Desktop/data/dirac/systemsize48/localsize8/L2norm/")
size48local8gmera1thirddiffL2norm = fioload("gmera1thirdalldiffL2norm")["re"]
size48local8gmera2thirddiffL2norm = fioload("gmera2thirdalldiffL2norm")["re"]
size48local8gmeraalldiffL2normdirac = fioload("gmeraalldiffL2norm")["re"]

radius_list = [2,4,6,8]
size48L2normdirac_list = [size48local2gmeraalldiffL2normdirac,size48local4gmeraalldiffL2normdirac,size48local6gmeraalldiffL2normdirac,size48local8gmeraalldiffL2normdirac]
scatter(radius_list,log.(size48L2normdirac_list))

# Trivial 
#onsite0.05
fiodir("/Users/slwongag/Desktop/data/trivial/systemsize48/onsitepotential0.05-localsize2/L2norm/")
size48local2gmera1thirddiffL2normtrivialonsite005  = fioload("gmera1thirdalldiffL2norm")["re"]
size48local2gmera2thirddiffL2normtrivialonsite005 = fioload("gmera2thirdalldiffL2norm")["re"]
size48local2gmera3thirddiffL2normtrivialonsite005 = fioload("gmera3thirdalldiffL2norm")["re"]
size48local2gmera4thirddiffL2normtrivialonsite005 = fioload("gmera4thirdalldiffL2norm")["re"]
size48local2gmeraalldiffL2normtrivialonsite005 = fioload("gmeraalldiffL2norm")["re"]

fiodir("/Users/slwongag/Desktop/data/trivial/systemsize48/onsitepotential0.05-localsize4/L2norm/")
size48local4gmera1thirddiffL2normtrivialonsite005  = fioload("gmera1thirdalldiffL2norm")["re"]
size48local4gmera2thirddiffL2normtrivialonsite005 = fioload("gmera2thirdalldiffL2norm")["re"]
size48local4gmera3thirddiffL2normtrivialonsite005 = fioload("gmera3thirdalldiffL2norm")["re"]
size48local4gmeraalldiffL2normtrivialonsite005 = fioload("gmeraalldiffL2norm")["re"]

fiodir("/Users/slwongag/Desktop/data/trivial/systemsize48/onsitepotential0.05-localsize6/L2norm/")
size48local6gmera1thirddiffL2normtrivialonsite005  = fioload("gmera1thirdalldiffL2norm")["re"]
size48local6gmera2thirddiffL2normtrivialonsite005 = fioload("gmera2thirdalldiffL2norm")["re"]
size48local6gmera3thirddiffL2normtrivialonsite005 = fioload("gmera3thirdalldiffL2norm")["re"]
size48local6gmeraalldiffL2normtrivialonsite005 = fioload("gmeraalldiffL2norm")["re"]

fiodir("/Users/slwongag/Desktop/data/trivial/systemsize48/onsitepotential0.05-localsize8/L2norm/")
size48local8gmera1thirddiffL2normtrivialonsite005  = fioload("gmera1thirdalldiffL2norm")["re"]
size48local8gmera2thirddiffL2normtrivialonsite005 = fioload("gmera2thirdalldiffL2norm")["re"]
size48local8gmeraalldiffL2normtrivialonsite005 = fioload("gmeraalldiffL2norm")["re"]

size48L2normtrivialonsite005_list = [size48local2gmeraalldiffL2normtrivialonsite005,size48local4gmeraalldiffL2normtrivialonsite005,size48local6gmeraalldiffL2normtrivialonsite005,size48local8gmeraalldiffL2normtrivialonsite005]
scatter!(radius_list,log.(size48L2normtrivialonsite005_list))

#onsite0.1
fiodir("/Users/slwongag/Desktop/data/trivial/systemsize48/onsitepotential0.1-localsize2/L2norm/")
size48local2gmera1thirddiffL2normtrivialonsite01  = fioload("gmera1thirdalldiffL2norm")["re"]
size48local2gmera2thirddiffL2normtrivialonsite01 = fioload("gmera2thirdalldiffL2norm")["re"]
size48local2gmera3thirddiffL2normtrivialonsite01 = fioload("gmera3thirdalldiffL2norm")["re"]
size48local2gmera4thirddiffL2normtrivialonsite01 = fioload("gmera4thirdalldiffL2norm")["re"]
size48local2gmeraalldiffL2normtrivialonsite01 = fioload("gmeraalldiffL2norm")["re"]

fiodir("/Users/slwongag/Desktop/data/trivial/systemsize48/onsitepotential0.1-localsize4/L2norm/")
size48local4gmera1thirddiffL2normtrivialonsite01  = fioload("gmera1thirdalldiffL2norm")["re"]
size48local4gmera2thirddiffL2normtrivialonsite01 = fioload("gmera2thirdalldiffL2norm")["re"]
size48local4gmera3thirddiffL2normtrivialonsite01 = fioload("gmera3thirdalldiffL2norm")["re"]
size48local4gmeraalldiffL2normtrivialonsite01 = fioload("gmeraalldiffL2norm")["re"]

fiodir("/Users/slwongag/Desktop/data/trivial/systemsize48/onsitepotential0.1-localsize6/L2norm/")
size48local6gmera1thirddiffL2normtrivialonsite01  = fioload("gmera1thirdalldiffL2norm")["re"]
size48local6gmera2thirddiffL2normtrivialonsite01 = fioload("gmera2thirdalldiffL2norm")["re"]
size48local6gmera3thirddiffL2normtrivialonsite01 = fioload("gmera3thirdalldiffL2norm")["re"]
size48local6gmeraalldiffL2normtrivialonsite01 = fioload("gmeraalldiffL2norm")["re"]

fiodir("/Users/slwongag/Desktop/data/trivial/systemsize48/onsitepotential0.1-localsize8/L2norm/")
size48local8gmera1thirddiffL2normtrivialonsite01  = fioload("gmera1thirdalldiffL2norm")["re"]
size48local8gmera2thirddiffL2normtrivialonsite01 = fioload("gmera2thirdalldiffL2norm")["re"]
size48local8gmeraalldiffL2normtrivialonsite01 = fioload("gmeraalldiffL2norm")["re"]

size48L2normtrivialonsite01_list = [size48local2gmeraalldiffL2normtrivialonsite01,size48local4gmeraalldiffL2normtrivialonsite01,size48local6gmeraalldiffL2normtrivialonsite01,size48local8gmeraalldiffL2normtrivialonsite01]
scatter!(radius_list,log.(size48L2normtrivialonsite01_list))

#onsite0.15
fiodir("/Users/slwongag/Desktop/data/trivial/systemsize48/onsitepotential0.15-localsize2/L2norm/")
size48local2gmera1thirddiffL2normtrivialonsite015  = fioload("gmera1thirdalldiffL2norm")["re"]
size48local2gmera2thirddiffL2normtrivialonsite015 = fioload("gmera2thirdalldiffL2norm")["re"]
size48local2gmera3thirddiffL2normtrivialonsite015 = fioload("gmera3thirdalldiffL2norm")["re"]
size48local2gmera4thirddiffL2normtrivialonsite015 = fioload("gmera4thirdalldiffL2norm")["re"]
size48local2gmeraalldiffL2normtrivialonsite015 = fioload("gmeraalldiffL2norm")["re"]

fiodir("/Users/slwongag/Desktop/data/trivial/systemsize48/onsitepotential0.15-localsize4/L2norm/")
size48local4gmera1thirddiffL2normtrivialonsite015  = fioload("gmera1thirdalldiffL2norm")["re"]
size48local4gmera2thirddiffL2normtrivialonsite015 = fioload("gmera2thirdalldiffL2norm")["re"]
size48local4gmera3thirddiffL2normtrivialonsite015 = fioload("gmera3thirdalldiffL2norm")["re"]
size48local4gmeraalldiffL2normtrivialonsite015 = fioload("gmeraalldiffL2norm")["re"]

fiodir("/Users/slwongag/Desktop/data/trivial/systemsize48/onsitepotential0.15-localsize6/L2norm/")
size48local6gmera1thirddiffL2normtrivialonsite015  = fioload("gmera1thirdalldiffL2norm")["re"]
size48local6gmera2thirddiffL2normtrivialonsite015 = fioload("gmera2thirdalldiffL2norm")["re"]
size48local6gmera3thirddiffL2normtrivialonsite015 = fioload("gmera3thirdalldiffL2norm")["re"]
size48local6gmeraalldiffL2normtrivialonsite015 = fioload("gmeraalldiffL2norm")["re"]

fiodir("/Users/slwongag/Desktop/data/trivial/systemsize48/onsitepotential0.15-localsize8/L2norm/")
size48local8gmera1thirddiffL2normtrivialonsite015  = fioload("gmera1thirdalldiffL2norm")["re"]
size48local8gmera2thirddiffL2normtrivialonsite015 = fioload("gmera2thirdalldiffL2norm")["re"]
size48local8gmeraalldiffL2normtrivialonsite015 = fioload("gmeraalldiffL2norm")["re"]

size48L2normtrivialonsite015_list = [size48local2gmeraalldiffL2normtrivialonsite015,size48local4gmeraalldiffL2normtrivialonsite015,size48local6gmeraalldiffL2normtrivialonsite015,size48local8gmeraalldiffL2normtrivialonsite015]
scatter!(radius_list,log.(size48L2normtrivialonsite015_list))

#onsite0.2
fiodir("/Users/slwongag/Desktop/data/trivial/systemsize48/onsitepotential0.2-localsize2/L2norm/")
size48local2gmera1thirddiffL2normtrivialonsite02  = fioload("gmera1thirdalldiffL2norm")["re"]
size48local2gmera2thirddiffL2normtrivialonsite02 = fioload("gmera2thirdalldiffL2norm")["re"]
size48local2gmera3thirddiffL2normtrivialonsite02 = fioload("gmera3thirdalldiffL2norm")["re"]
size48local2gmera4thirddiffL2normtrivialonsite02 = fioload("gmera4thirdalldiffL2norm")["re"]
size48local2gmeraalldiffL2normtrivialonsite02 = fioload("gmeraalldiffL2norm")["re"]

fiodir("/Users/slwongag/Desktop/data/trivial/systemsize48/onsitepotential0.2-localsize4/L2norm/")
size48local4gmera1thirddiffL2normtrivialonsite02  = fioload("gmera1thirdalldiffL2norm")["re"]
size48local4gmera2thirddiffL2normtrivialonsite02 = fioload("gmera2thirdalldiffL2norm")["re"]
size48local4gmera3thirddiffL2normtrivialonsite02 = fioload("gmera3thirdalldiffL2norm")["re"]
size48local4gmeraalldiffL2normtrivialonsite02 = fioload("gmeraalldiffL2norm")["re"]

fiodir("/Users/slwongag/Desktop/data/trivial/systemsize48/onsitepotential0.2-localsize6/L2norm/")
size48local6gmera1thirddiffL2normtrivialonsite02  = fioload("gmera1thirdalldiffL2norm")["re"]
size48local6gmera2thirddiffL2normtrivialonsite02 = fioload("gmera2thirdalldiffL2norm")["re"]
size48local6gmera3thirddiffL2normtrivialonsite02 = fioload("gmera3thirdalldiffL2norm")["re"]
size48local6gmeraalldiffL2normtrivialonsite02 = fioload("gmeraalldiffL2norm")["re"]

fiodir("/Users/slwongag/Desktop/data/trivial/systemsize48/onsitepotential0.2-localsize8/L2norm/")
size48local8gmera1thirddiffL2normtrivialonsite02  = fioload("gmera1thirdalldiffL2norm")["re"]
size48local8gmera2thirddiffL2normtrivialonsite02 = fioload("gmera2thirdalldiffL2norm")["re"]
size48local8gmeraalldiffL2normtrivialonsite02 = fioload("gmeraalldiffL2norm")["re"]

size48L2normtrivialonsite02_list = [size48local2gmeraalldiffL2normtrivialonsite02,size48local4gmeraalldiffL2normtrivialonsite02,size48local6gmeraalldiffL2normtrivialonsite02,size48local8gmeraalldiffL2normtrivialonsite02]
scatter!(radius_list,log.(size48L2normtrivialonsite02_list))

# Chern
#nn 0.1im
fiodir("/Users/slwongag/Desktop/data/chern/systemsize48/nnhopping0.1im-localsize2/L2norm/")
size48local2gmera1thirddiffL2normchern01im = fioload("gmera1thirdalldiffL2norm")["re"]
size48local2gmera2thirddiffL2normchern01im = fioload("gmera2thirdalldiffL2norm")["re"]
size48local2gmera3thirddiffL2normchern01im = fioload("gmera3thirdalldiffL2norm")["re"]
size48local2gmera4thirddiffL2normchern01im = fioload("gmera4thirdalldiffL2norm")["re"]
size48local2gmeraalldiffL2normchernnn01im = fioload("gmeraalldiffL2norm")["re"]

fiodir("/Users/slwongag/Desktop/data/chern/systemsize48/nnhopping0.1im-localsize4/L2norm")
size48local4gmera1thirddiffL2normchern01im = fioload("gmera1thirdalldiffL2norm")["re"]
size48local4gmera2thirddiffL2normchern01im = fioload("gmera2thirdalldiffL2norm")["re"]
size48local4gmera3thirddiffL2normchern01im = fioload("gmera3thirdalldiffL2norm")["re"]
size48local4gmeraalldiffL2normchernnn01im = fioload("gmeraalldiffL2norm")["re"]

fiodir("/Users/slwongag/Desktop/data/chern/systemsize48/nnhopping0.1im-localsize6/L2norm")
size48local6gmera1thirddiffL2normchern01im = fioload("gmera1thirdalldiffL2norm")["re"]
size48local6gmera2thirddiffL2normchern01im = fioload("gmera2thirdalldiffL2norm")["re"]
size48local6gmera3thirddiffL2normchern01im = fioload("gmera3thirdalldiffL2norm")["re"]
size48local6gmeraalldiffL2normchernnn01im = fioload("gmeraalldiffL2norm")["re"]

fiodir("/Users/slwongag/Desktop/data/chern/systemsize48/nnhopping0.1im-localsize8/L2norm")
size48local8gmera1thirddiffL2normchern01im = fioload("gmera1thirdalldiffL2norm")["re"]
size48local8gmera1thirddiffL2normchern01im = fioload("gmera2thirdalldiffL2norm")["re"]
size48local8gmeraalldiffL2normchernnn01im = fioload("gmeraalldiffL2norm")["re"]

size48L2normchernnn01im_list = [size48local2gmeraalldiffL2normchernnn01im,size48local4gmeraalldiffL2normchernnn01im,size48local6gmeraalldiffL2normchernnn01im,size48local8gmeraalldiffL2normchernnn01im]
scatter!(radius_list,log.(size48L2normchernnn01im_list))

#nn 0.2im
fiodir("/Users/slwongag/Desktop/data/chern/systemsize48/nnhopping0.2im-localsize2/L2norm/")
size48local2gmera1thirddiffL2normchern02im = fioload("gmera1thirdalldiffL2norm")["re"]
size48local2gmera2thirddiffL2normchern02im = fioload("gmera2thirdalldiffL2norm")["re"]
size48local2gmera3thirddiffL2normchern02im = fioload("gmera3thirdalldiffL2norm")["re"]
size48local2gmera4thirddiffL2normchern02im = fioload("gmera4thirdalldiffL2norm")["re"]
size48local2gmeraalldiffL2normchernnn02im = fioload("gmeraalldiffL2norm")["re"]

fiodir("/Users/slwongag/Desktop/data/chern/systemsize48/nnhopping0.2im-localsize4/L2norm")
size48local4gmera1thirddiffL2normchern02im = fioload("gmera1thirdalldiffL2norm")["re"]
size48local4gmera2thirddiffL2normchern02im = fioload("gmera2thirdalldiffL2norm")["re"]
size48local4gmera3thirddiffL2normchern02im = fioload("gmera3thirdalldiffL2norm")["re"]
size48local4gmeraalldiffL2normchernnn02im = fioload("gmeraalldiffL2norm")["re"]

fiodir("/Users/slwongag/Desktop/data/chern/systemsize48/nnhopping0.2im-localsize6/L2norm")
size48local6gmera1thirddiffL2normchern02im = fioload("gmera1thirdalldiffL2norm")["re"]
size48local6gmera2thirddiffL2normchern02im = fioload("gmera2thirdalldiffL2norm")["re"]
size48local6gmera3thirddiffL2normchern02im = fioload("gmera3thirdalldiffL2norm")["re"]
size48local6gmeraalldiffL2normchernnn02im = fioload("gmeraalldiffL2norm")["re"]

fiodir("/Users/slwongag/Desktop/data/chern/systemsize48/nnhopping0.2im-localsize8/L2norm")
size48local8gmera1thirddiffL2normchern02im = fioload("gmera1thirdalldiffL2norm")["re"]
size48local8gmera1thirddiffL2normchern02im = fioload("gmera2thirdalldiffL2norm")["re"]
size48local8gmeraalldiffL2normchernnn02im = fioload("gmeraalldiffL2norm")["re"]

size48L2normchernnn02im_list = [size48local2gmeraalldiffL2normchernnn02im,size48local4gmeraalldiffL2normchernnn02im,size48local6gmeraalldiffL2normchernnn02im,size48local8gmeraalldiffL2normchernnn02im]
scatter!(radius_list,log.(size48L2normchernnn02im_list))

#nn 0.3im
fiodir("/Users/slwongag/Desktop/data/chern/systemsize48/nnhopping0.3im-localsize2/L2norm/")
size48local2gmera1thirddiffL2normchern03im = fioload("gmera1thirdalldiffL2norm")["re"]
size48local2gmera2thirddiffL2normchern03im = fioload("gmera2thirdalldiffL2norm")["re"]
size48local2gmera3thirddiffL2normchern03im = fioload("gmera3thirdalldiffL2norm")["re"]
size48local2gmera4thirddiffL2normchern03im = fioload("gmera4thirdalldiffL2norm")["re"]
size48local2gmeraalldiffL2normchernnn03im = fioload("gmeraalldiffL2norm")["re"]

fiodir("/Users/slwongag/Desktop/data/chern/systemsize48/nnhopping0.3im-localsize4/L2norm")
size48local4gmera1thirddiffL2normchern03im = fioload("gmera1thirdalldiffL2norm")["re"]
size48local4gmera2thirddiffL2normchern03im = fioload("gmera2thirdalldiffL2norm")["re"]
size48local4gmera3thirddiffL2normchern03im = fioload("gmera3thirdalldiffL2norm")["re"]
size48local4gmeraalldiffL2normchernnn03im = fioload("gmeraalldiffL2norm")["re"]

fiodir("/Users/slwongag/Desktop/data/chern/systemsize48/nnhopping0.3im-localsize6/L2norm")
size48local6gmera1thirddiffL2normchern03im = fioload("gmera1thirdalldiffL2norm")["re"]
size48local6gmera2thirddiffL2normchern03im = fioload("gmera2thirdalldiffL2norm")["re"]
size48local6gmera3thirddiffL2normchern03im = fioload("gmera3thirdalldiffL2norm")["re"]
size48local6gmeraalldiffL2normchernnn03im = fioload("gmeraalldiffL2norm")["re"]

fiodir("/Users/slwongag/Desktop/data/chern/systemsize48/nnhopping0.3im-localsize8/L2norm")
size48local8gmera1thirddiffL2normchern03im = fioload("gmera1thirdalldiffL2norm")["re"]
size48local8gmera1thirddiffL2normchern03im = fioload("gmera2thirdalldiffL2norm")["re"]
size48local8gmeraalldiffL2normchernnn03im = fioload("gmeraalldiffL2norm")["re"]

size48L2normchernnn03im_list = [size48local2gmeraalldiffL2normchernnn03im,size48local4gmeraalldiffL2normchernnn03im,size48local6gmeraalldiffL2normchernnn03im,size48local8gmeraalldiffL2normchernnn03im]
scatter!(radius_list,log.(size48L2normchernnn03im_list))

#nn 0.4im
fiodir("/Users/slwongag/Desktop/data/chern/systemsize48/nnhopping0.4im-localsize2/L2norm/")
size48local2gmera1thirddiffL2normchern04im = fioload("gmera1thirdalldiffL2norm")["re"]
size48local2gmera2thirddiffL2normchern04im = fioload("gmera2thirdalldiffL2norm")["re"]
size48local2gmera3thirddiffL2normchern04im = fioload("gmera3thirdalldiffL2norm")["re"]
size48local2gmera4thirddiffL2normchern04im = fioload("gmera4thirdalldiffL2norm")["re"]
size48local2gmeraalldiffL2normchernnn04im = fioload("gmeraalldiffL2norm")["re"]

fiodir("/Users/slwongag/Desktop/data/chern/systemsize48/nnhopping0.4im-localsize4/L2norm")
size48local4gmera1thirddiffL2normchern04im = fioload("gmera1thirdalldiffL2norm")["re"]
size48local4gmera2thirddiffL2normchern04im = fioload("gmera2thirdalldiffL2norm")["re"]
size48local4gmera3thirddiffL2normchern04im = fioload("gmera3thirdalldiffL2norm")["re"]
size48local4gmeraalldiffL2normchernnn04im = fioload("gmeraalldiffL2norm")["re"]

fiodir("/Users/slwongag/Desktop/data/chern/systemsize48/nnhopping0.4im-localsize6/L2norm")
size48local6gmera1thirddiffL2normchern04im = fioload("gmera1thirdalldiffL2norm")["re"]
size48local6gmera2thirddiffL2normchern04im = fioload("gmera2thirdalldiffL2norm")["re"]
size48local6gmera3thirddiffL2normchern04im = fioload("gmera3thirdalldiffL2norm")["re"]
size48local6gmeraalldiffL2normchernnn04im = fioload("gmeraalldiffL2norm")["re"]

fiodir("/Users/slwongag/Desktop/data/chern/systemsize48/nnhopping0.4im-localsize8/L2norm")
size48local8gmera1thirddiffL2normchern04im = fioload("gmera1thirdalldiffL2norm")["re"]
size48local8gmera1thirddiffL2normchern04im = fioload("gmera2thirdalldiffL2norm")["re"]
size48local8gmeraalldiffL2normchernnn04im = fioload("gmeraalldiffL2norm")["re"]

size48L2normchernnn04im_list = [size48local2gmeraalldiffL2normchernnn04im,size48local4gmeraalldiffL2normchernnn04im,size48local6gmeraalldiffL2normchernnn04im,size48local8gmeraalldiffL2normchernnn04im]
scatter!(radius_list,log.(size48L2normchernnn04im_list))

