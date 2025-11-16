struct Bond
    i :: Int
    j :: Int
    r_i :: Tuple{Float64, Float64} # Cartesian coordinates of site i
    r_j :: Tuple{Float64, Float64} # Cartesian coordinates of site j (No periodic wrapping)
end

struct TriangularLatticeModel
    Lx :: Int
    Ly :: Int
    order :: Int
    J₁ :: Float64
    J₂ :: Float64
    Δ  :: Float64
    h  :: Float64
end

TriangularLatticeModel(Lx::Int, Ly::Int, J₁::Float64, J₂::Float64, Δ::Float64, h::Float64; order::Int=1) = TriangularLatticeModel(Lx, Ly, order, J₁, J₂, Δ, h)

struct CorrelationData
    correlations :: Array{ComplexF64, 4}
    r0s :: Tuple{Vararg{Int}}
    combiners :: Tuple{Vararg{Tuple{Int, Int}}}
    dt :: Float64
    tf :: Float64
    Lx :: Int
    Ly :: Int
end

struct IntensityData
    intensities :: Array{Float64, 2}
    Sxx :: Array{Float64, 2}
    Szz :: Array{Float64, 2}
    q_points :: Matrix{Vector{Float64}}
    energies_full :: Vector{Float64}
    Lx :: Int
    Ly :: Int
    measure :: Symbol
end