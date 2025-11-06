module TriSpinDPlotExt

using TriSpinD
import TriSpinD: siteindex, inverse_siteindex, triangular_lattice_bonds, cartesian_coordinates
using GLMakie

function TriSpinD.plot_triangular_cluster(Lx::Int, Ly::Int; order::Int=1)
    bonds_1nn, bonds_2nn = triangular_lattice_bonds(Lx, Ly; order=order)

    with_theme(theme_latexfonts()) do
        fig = Figure()
        ax = Axis(fig[1,1], title="Triangular Lattice Cluster", titlesize=30, aspect=DataAspect())
        hidedecorations!(ax)
        hidexdecorations!(ax)
        hideydecorations!(ax)
        for x in 1:Lx, y in 1:Ly
            site = siteindex(Lx, Ly, x, y; order=order)
            r = cartesian_coordinates(x, y)
            scatter!(ax, r[1], r[2], markersize=10, color=:blue)
            text!(ax, r[1]+0.1, r[2]+0.02; text="$site", color=:black, fontsize=25)
        end

        for bond in bonds_1nn
            (; r_i, r_j) = bond
            lines!(ax, [r_i[1], r_j[1]], [r_i[2], r_j[2]], color=:red)
        end

        for bond in bonds_2nn
            (; r_i, r_j) = bond
            lines!(ax, [r_i[1], r_j[1]], [r_i[2], r_j[2]], color=:blue, linestyle=:dash)
        end
        fig
    end
end

function TriSpinD.plot_Sz_expectations(Lx::Int, Ly::Int, Szs::Vector{Float64}; order::Int=1)
    with_theme(theme_latexfonts()) do
        fig = Figure()
        ax = Axis(fig[1,1], title="Sᶻ Expectation Values", titlesize=30, aspect=DataAspect())
        hidedecorations!(ax)
        hidexdecorations!(ax)
        hideydecorations!(ax)
        for x in 1:Lx, y in 1:Ly
            site = siteindex(Lx, Ly, x, y; order=order)
            r = cartesian_coordinates(x, y)
            scatter!(ax, r[1], r[2], markersize=20, color=Szs[site], colorrange=(-0.5, 0.5), colormap=:bwr)
            text!(ax, r[1]+0.1, r[2]+0.02; text="$site", color=:black, fontsize=25)
        end
        xlims!(ax, -0.1, Lx+0.1)
        ylims!(ax, -0.1, √3/2*Ly+0.1)
        Colorbar(fig[1,2]; colorrange=(-0.5, 0.5), colormap=:bwr)
        fig
    end
end

end