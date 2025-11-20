module TriSpinDPlotExt

using TriSpinD
import TriSpinD: siteindex, inverse_siteindex, triangular_lattice_bonds, cartesian_coordinates, available_q_points
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
        ax = Axis(fig[1,1], title=L"\langle S^z \rangle", titlesize=30, aspect=DataAspect())
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

function TriSpinD.plot_dssf_heatmap(intensity_data::IntensityData, q_indices::AbstractVector{Int}; xlabel, xticks, energy_min=0.0, energy_max=30.0, colorrange=nothing, colormap=:viridis)
    (; intensities, energies_full, Lx, Ly, measure) = intensity_data

    Nsites = Lx * Ly
    all(1 ≤ x ≤ Nsites for x in q_indices) || error("q_indices must be between 1 and $Nsites.")
    num_qs = length(q_indices)
    inten_cuts = zeros(length(q_indices), length(energies_full))
    for (i, q_index) in enumerate(q_indices)
        inten_cuts[i, :] = intensities[:, q_index]
    end
    isnothing(colorrange) && (colorrange = (0.0, maximum(inten_cuts)))

    title_str = "DSSF-" * string(measure) * " Heatmap"

    with_theme(theme_latexfonts()) do
        fig = Figure()
        ax  = Axis(fig[1,1], title=title_str, xlabel=xlabel, ylabel=L"E/J", titlesize=30, xlabelsize=30, ylabelsize=30, xticks=xticks, xticklabelsize=25, yticklabelsize=25)
        cm = heatmap!(ax, 1:num_qs, energies_full, inten_cuts, colormap=colormap, colorrange=colorrange)
        xlims!(ax, 1, num_qs)
        ylims!(ax, energy_min, energy_max)
        Colorbar(fig[1,2], cm)
        fig
    end

end

function TriSpinD.plot_available_q_points(Lx::Int, Ly::Int)
    b₁ = 2π * [1, 1/√3]
    b₂ = 2π * [0, 2/√3]

    BZ_vertices = [2b₁/3-b₂/3, b₁/3+b₂/3, -b₁/3+2b₂/3, -2b₁/3+b₂/3, -b₁/3-b₂/3, b₁/3-2b₂/3, 2b₁/3-b₂/3]

    q_points = available_q_points(Lx, Ly)
    with_theme(theme_latexfonts()) do
        fig = Figure()
        ax = Axis(fig[1,1], title="Available q-points", xlabel=L"q_x", ylabel=L"q_y", titlesize=30, xlabelsize=30, ylabelsize=30, aspect=DataAspect())
        for (iq, q) in enumerate(q_points)
            qx, qy = q
            scatter!(ax, qx, qy, markersize=10, color=:blue)
            text!(ax, qx+0.1, qy+0.02; text="$(iq)", color=:black, fontsize=20)
        end
        lines!(ax, [v[1] for v in BZ_vertices], [v[2] for v in BZ_vertices], color=:black)
        text!(ax, 0, 0; text=L"\Gamma", color=:purple, fontsize=20, align=(:right, :top))
        text!(ax, 4π/3, 0; text=L"K", color=:purple, fontsize=20, align=(:left, :top))
        text!(ax, 2π/3, 2π/√3; text=L"K", color=:purple, fontsize=20, align=(:right, :top))
        text!(ax, π, π/√3; text=L"M", color=:purple, fontsize=20, align=(:right, :top))
        text!(ax, 0, 2π/√3; text=L"M", color=:purple, fontsize=20, align=(:right, :top))
        xlims!(ax, -2π, 2π)
        ylims!(ax, -2π/√3-0.6, 2π/√3+0.6)
        fig
    end
end


end