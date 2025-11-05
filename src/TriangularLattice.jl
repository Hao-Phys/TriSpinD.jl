# Here, we implement the XC boundary condition for triangular lattice
"""
    siteindex(Lx, Ly, x, y; order::Int=1)

Map the 2D coordinate `(x, y)` on an Lx×Ly triangular XC lattice
to a 1D site index according to the winding order.

- `order = 1`: odd-y sites first (1,3,5,...) then even-y (2,4,6,...)
  within each column, matching Hao’s "order = 1" figure.
- `order = 2`: simple column-major order (bottom→top).
"""
function siteindex(Lx::Int, Ly::Int, x::Int, y::Int; order::Int = 1)
    @boundscheck 1 ≤ x ≤ Lx || error("x out of range 1:$Lx")
    @boundscheck 1 ≤ y ≤ Ly || error("y out of range 1:$Ly")

    if order == 1
        n_odd = (Ly + 1) ÷ 2
        if isodd(y)
            site = (x - 1) * Ly + (y + 1) ÷ 2
        else
            site = (x - 1) * Ly + n_odd + y ÷ 2
        end
    elseif order == 2
        site = (x - 1) * Ly + y
    else
        error("Invalid order = $order. Use 1 or 2.")
    end

    return site
end

"""
    inverse_siteindex(Lx, Ly, site; order::Int=1)

Inverse of `siteindex`: return `(x, y)` for a given site index.
"""
function inverse_siteindex(Lx::Int, Ly::Int, site::Int; order::Int = 1)
    N = Lx * Ly
    @boundscheck 1 ≤ site ≤ N || error("site out of range 1:$N")

    x = div(site - 1, Ly) + 1
    pos = site - (x - 1) * Ly

    if order == 1
        n_odd = (Ly + 1) ÷ 2
        if pos ≤ n_odd
            y = 2 * pos - 1
        else
            y = 2 * (pos - n_odd)
        end
    elseif order == 2
        y = pos
    else
        error("Invalid order = $order. Use 1 or 2.")
    end

    return x, y
end

function cartesian_coordinates(x::Int, y::Int)
    r_x = x - 1
    if iseven(y)
        r_x += 0.5
    end
    r_y = (y - 1) * √3 / 2
    return (r_x, r_y)
end

function triangular_lattice_bonds(Lx::Int, Ly::Int; order::Int=1)
    bonds_1nn = Bond[]
    bonds_2nn = Bond[]
    @assert order == 1 || order == 2 "Invalid order = $order. Use 1 or 2."
    for x in 1:Lx, y in 1:Ly
        i = siteindex(Lx, Ly, x, y; order=order)
        r_i = cartesian_coordinates(x, y)

        # 1nn
        if x < Lx
            j1 = siteindex(Lx, Ly, x+1, y; order=order)
            r_j1 = cartesian_coordinates(x+1, y)
            push!(bonds_1nn, Bond(i, j1, r_i, r_j1))
        end
        if x > 1 || iseven(y)
            j2 = isodd(y) ? siteindex(Lx, Ly, x-1, mod1(y+1, Ly); order=order) : siteindex(Lx, Ly, x, mod1(y+1, Ly); order=order)
            r_j2 = cartesian_coordinates(isodd(y) ? x-1 : x, y+1)
            push!(bonds_1nn, Bond(i, j2, r_i, r_j2))
        end
        if x < Lx || isodd(y)
            j3 = iseven(y) ? siteindex(Lx, Ly, x+1, mod1(y+1, Ly); order=order) : siteindex(Lx, Ly, x, mod1(y+1, Ly); order=order)
            r_j3 = cartesian_coordinates(iseven(y) ? x+1 : x, y+1)
            push!(bonds_1nn, Bond(i, j3, r_i, r_j3))
        end

        # 2nn
        if x < Lx - 1
            j21 = isodd(y) ? siteindex(Lx, Ly, x+1, mod1(y+1, Ly); order=order) : siteindex(Lx, Ly, x+2, mod1(y+1, Ly); order=order)
            r_j21 = cartesian_coordinates(isodd(y) ? x+1 : x+2, y+1)
            push!(bonds_2nn, Bond(i, j21, r_i, r_j21))
        elseif x == Lx - 1
            if isodd(y)
                j21 = siteindex(Lx, Ly, x+1, mod1(y+1, Ly); order=order)
                r_j21 = cartesian_coordinates(x+1, y+1)
                push!(bonds_2nn, Bond(i, j21, r_i, r_j21))
            end
        end

        j22 = siteindex(Lx, Ly, x, mod1(y+2, Ly); order=order)
        r_j22 = cartesian_coordinates(x, y+2)
        push!(bonds_2nn, Bond(i, j22, r_i, r_j22))
        if x > 2
            j23 = iseven(y) ? siteindex(Lx, Ly, x-1, mod1(y+1, Ly); order=order) : siteindex(Lx, Ly, x-2, mod1(y+1, Ly); order=order)
            r_j23 = cartesian_coordinates(iseven(y) ? x-1 : x-2, y+1)
            push!(bonds_2nn, Bond(i, j23, r_i, r_j23))
        elseif x == 2
            if iseven(y)
                j23 = siteindex(Lx, Ly, x-1, mod1(y+1, Ly); order=order)
                r_j23 = cartesian_coordinates(x-1, y+1)
                push!(bonds_2nn, Bond(i, j23, r_i, r_j23))
            end
        end
    end

    return bonds_1nn, bonds_2nn
end
