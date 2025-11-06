function construct_mpo(tlm::TriangularLatticeModel, sites)
    (; Lx, Ly, order, J₁, J₂, Δ, h) = tlm

    Ns = Lx * Ly
    bonds_1nn, bonds_2nn = triangular_lattice_bonds(Lx, Ly; order=order)
    os = AutoMPO()

    for bond in bonds_1nn
        (; i, j) = bond
        os += J₁, "Sx", i, "Sx", j
        os += J₁, "Sy", i, "Sy", j
        os += J₁ * Δ, "Sz", i, "Sz", j
    end

    for bond in bonds_2nn
        (; i, j) = bond
        os += J₂, "Sx", i, "Sx", j
        os += J₂, "Sy", i, "Sy", j
        os += J₂ * Δ, "Sz", i, "Sz", j
    end

    for i in 1:Ns
        os += h, "Sz", i
    end

    H_mpo = MPO(os, sites)

    return H_mpo
end

function initialize_UUD(tlm::TriangularLatticeModel)
    (; Lx, Ly, order) = tlm
    Ns = Lx * Ly
    sites = siteinds("S=1/2", Ns; conserve_qns=true)
    ψ0 = MPS(sites)

    for x in 1:Lx, y in 1:Ly
        i = siteindex(Lx, Ly, x, y; order)
        s = sites[i]
        if isodd(y)
            if mod(x, 3) == 0
                ψ0[i] = ITensors.state(s, "Dn")
            else
                ψ0[i] = ITensors.state(s, "Up")
            end
        else
            if mod(x, 3) == 1
                ψ0[i] = ITensors.state(s, "Dn")
            else
                ψ0[i] = ITensors.state(s, "Up")
            end
        end
    end

    return ψ0, sites
end

function initialize_UUUD(tlm::TriangularLatticeModel)
    (; Lx, Ly, order) = tlm
    Ns = Lx * Ly
    sites = siteinds("S=1/2", Ns; conserve_qns=true)
    ψ0 = MPS(sites)

    for x in 1:Lx, y in 1:Ly
        i = siteindex(Lx, Ly, x, y; order)
        s = sites[i]
        if isodd(y)
            if mod(x, 2) == 1
                ψ0[i] = ITensors.state(s, "Up")
            else
                ψ0[i] = ITensors.state(s, "Dn")
            end
        else
            ψ0[i] = ITensors.state(s, "Up")
        end
    end

    return ψ0, sites
end

"""
    dmrg_gs(tlm::TriangularLatticeModel, ψ0, sites; kwargs...)
"""
function dmrg_gs(tlm::TriangularLatticeModel, ψ0, sites; kwargs...)
    H_mpo = construct_mpo(tlm, sites)
    Ns = length(sites)

    E_gs, ψ_gs = dmrg(H_mpo, ψ0; kwargs...)

    @show "Ground state energy per site: ", E_gs / Ns

    Sxs = [expect(ψ_gs, "Sx"; sites=i) for i in 1:Ns]
    @show "Sx expectation values: ", Sxs
    println("Sx expectation values: ", Sxs)
    Sys = [expect(ψ_gs, "Sy"; sites=i) for i in 1:Ns]
    @show "Sy expectation values: ", Sys
    Szs = [expect(ψ_gs, "Sz"; sites=i) for i in 1:Ns]
    @show "Sz expectation values: ", Szs

    return ψ_gs, E_gs, sites, H_mpo
end

function correlation_function_aux(j, E_gs, ψ_gs, sites, H_mpo, dt, tf, α, β; tdvp_kwargs...)
    Ns = length(sites)
    spinops = ("Sx", "Sy", "Sz")
    Siαs = [op(spinops[α], sites, i) for i in 1:Ns]
    Sjβ  = op(spinops[β], sites, j)

    nt = Int(round(tf / dt)) + 1
    corrs = zeros(ComplexF64, nt, Ns)

    Si0_exps = zeros(ComplexF64, Ns)
    Sj0_exp = expect(ψ_gs, spinops[β]; sites=j)

    ψt = apply(Sjβ, ψ_gs)

    for i in 1:Ns
        corrs[1, i] += inner(ψ_gs, apply(Siαs[i], ψt))
        Si0_exps[i] = expect(ψ_gs, spinops[α]; sites=i)
        corrs[1, i] -= Si0_exps[i] * Sj0_exp
    end

    for it in 2:nt
        @show it, α, β
        ψt = tdvp(H_mpo, -1im*dt, ψt; nsteps=1, tdvp_kwargs...)
        ψt *= exp(1im * dt * E_gs)
        for i in 1:Ns
            corrs[it, i] += inner(ψ_gs, apply(Siαs[i], ψt))
            corrs[it, i] -= Si0_exps[i] * Sj0_exp
        end
    end

    return corrs
end

"""
    correlation_function(tlm::TriangularLatticeModel, dt, tf, combiners::NTuple{Tuple{Int, Int}}; kwargs...)

Compute the time-dependent spin correlation functions ⟨Sᵅ(r, t) Sᵝ(r₀, 0)⟩ - ⟨Sᵅ(r)⟩⟨Sᵝ(r₀)⟩ for all sites r in the triangular lattice model `tlm`, where `r₀` is chosen to be the center site of sublattice A and B. Here The time evolution is performed up to time `tf` with time step `dt`. The spin component pairs (α, β) are specified in `combiners`, which is a tuple of tuples, e.g., `((1, 1), (2, 2), (3, 3))` for Sx-Sx, Sy-Sy, and Sz-Sz correlations.

Tips on the choice of `dt` and `tf`: Given a `max_energy` scale of interest, dt ~ (1/max_energy) * π; given a `energy_resolution` scale of interest, tf ~ (1/energy_resolution) * π.
"""
function correlation_function(tlm::TriangularLatticeModel, ψ0, sites, dt, tf, r0s::Tuple{Vararg{Int}}, combiners::Tuple{Vararg{Tuple{Int, Int}}}, dmrg_kwargs::NamedTuple; tdvp_kwargs::NamedTuple=(;))
    ψ_gs, E_gs, sites, H_mpo = dmrg_gs(tlm, ψ0, sites; dmrg_kwargs...)
    @show "Maximum bond dimension in ground state MPS: " maxlinkdim(ψ_gs)

    num_combiners = length(combiners)
    iters = [(combiner_index, r0_index) for combiner_index in 1:num_combiners for r0_index in 1:2]

    Ns = length(sites)
    nt = Int(round(tf / dt)) + 1
    correlations = zeros(ComplexF64, nt, Ns, length(combiners), length(r0s))

    Threads.@threads for i in eachindex(iters)
        combiner_index, r0_index = iters[i]
        α, β = combiners[combiner_index]
        r0 = r0s[r0_index]
        corrs = correlation_function_aux(r0, E_gs, ψ_gs, sites, H_mpo, dt, tf, α, β; tdvp_kwargs...)
        correlations[:, :, combiner_index, r0_index] = corrs
    end

    return CorrelationData(correlations, r0s, combiners, dt, tf, L1, L2)
end