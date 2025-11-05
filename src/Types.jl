struct Bond
    i :: Int
    j :: Int
    r_i :: Tuple{Float64, Float64} # Cartesian coordinates of site i
    r_j :: Tuple{Float64, Float64} # Cartesian coordinates of site j (No periodic wrapping)
end