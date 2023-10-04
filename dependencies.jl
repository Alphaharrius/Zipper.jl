using Pkg

Pkg.add("LinearAlgebra")
Pkg.add("Combinatorics")
Pkg.add("Statistics")
Pkg.add("Plotly")
Pkg.add("OrderedCollections")
Pkg.add("ColorTypes")
Pkg.add("Compat")
Pkg.Registry.add(RegistrySpec(url="https://github.com/wildart/BoffinStuff.git"))
Pkg.add("SmithNormalForm")
Pkg.add("SparseArrays")

Pkg.add("Revise")