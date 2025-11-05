module TriSpinD

using LinearAlgebra

include("Types.jl")
include("TriangularLattice.jl")

export plot_triangular_cluster
function plot_triangular_cluster(args...; kwarg...)
    _needs_glmakie(:plot_triangular_cluster)
end

_needs_glmakie(fname::Symbol) = error(string(
    fname, " requires GLMakie. Please `using GLMakie` first."
))

end
