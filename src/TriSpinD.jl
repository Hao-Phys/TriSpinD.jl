module TriSpinD

using LinearAlgebra
using ITensors
using ITensorMPS
using FFTW

include("Types.jl")
export TriangularLatticeModel, CorrelationData, IntensityData

include("TriangularLattice.jl")

export plot_triangular_cluster, plot_Sz_expectations, plot_dssf_heatmap, plot_available_q_points

function plot_triangular_cluster(args...; kwarg...)
    _needs_glmakie(:plot_triangular_cluster)
end

function plot_Sz_expectations(args...; kwarg...)
    _needs_glmakie(:plot_Sz_expectations)
end

function plot_dssf_heatmap(args...; kwarg...)
    _needs_glmakie(:plot_dssf_heatmap)
end

function plot_available_q_points(args...)
    _needs_glmakie(:available_q_points)
end

_needs_glmakie(fname::Symbol) = error(string(
    fname, " requires GLMakie. Please `using GLMakie` first."
))

include("ITensorTools.jl")
export construct_mpo, initialize_UUD, initialize_UUUD, initialize_polarized, dmrg_gs, correlation_function, Sz_expectations

include("FourierTransforms.jl")

include("DSSF.jl")
export dssf

end
