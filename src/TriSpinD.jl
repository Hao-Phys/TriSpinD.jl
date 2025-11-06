module TriSpinD

using LinearAlgebra
using ITensors
using ITensorMPS
using FFTW

include("Types.jl")
export TriangularLatticeModel, CorrelationData, IntensityData

include("TriangularLattice.jl")

export plot_triangular_cluster, plot_Sz_expectations

function plot_triangular_cluster(args...; kwarg...)
    _needs_glmakie(:plot_triangular_cluster)
end

function plot_Sz_expectations(args...; kwarg...)
    _needs_glmakie(:plot_Sz_expectations)
end

_needs_glmakie(fname::Symbol) = error(string(
    fname, " requires GLMakie. Please `using GLMakie` first."
))

include("ITensorTools.jl")
export construct_mpo, initialize_UUD, initialize_UUUD, dmrg_gs, correlation_function, Sz_expectations

include("PostProcessing.jl")

include("DSSF.jl")
export dssf

end
