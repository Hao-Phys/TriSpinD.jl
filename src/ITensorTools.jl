function construct_mpo(tlm::TriangularLatticeModel, sites)
    (; Lx, Ly, order, J₁, J₂, Δ, h) = tlm

    Ns = Lx * Ly
    bonds_1nn, bonds_2nn = triangular_lattice_bonds(Lx, Ly; order=order)
    os = AutoMPO()

    for bond in bonds_1nn
        (; i, j) = bond
        os += 0.5 * J₁, "S+", i, "S-", j
        os += 0.5 * J₁, "S-", i, "S+", j
        os += J₁ * Δ, "Sz", i, "Sz", j
    end

    for bond in bonds_2nn
        (; i, j) = bond
        os += 0.5 * J₂, "S+", i, "S-", j
        os += 0.5 * J₂, "S-", i, "S+", j
        os += J₂ * Δ, "Sz", i, "Sz", j
    end

    for i in 1:Ns
        os += -h, "Sz", i
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

function initialize_polarized(tlm::TriangularLatticeModel)
    (; Lx, Ly, order) = tlm
    Ns = Lx * Ly
    sites = siteinds("S=1/2", Ns; conserve_qns=true)
    ψ0 = MPS(sites)

    for i in 1:Ns
        ψ0[i] = ITensors.state(sites[i], "Up")
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

    return ψ_gs, E_gs, sites, H_mpo
end

function Sz_expectations(ψ, sites)
    Ns = length(sites)
    Szs = [expect(ψ, "Sz"; sites=i) for i in 1:Ns]
    return Szs
end

function correlation_function_aux(j, E_gs, ψ_gs, sites, H_mpo, dt, tf, α, β; tdvp_kwargs...)
    Ns = length(sites)
    spinops = ("S+", "S-", "Sz")
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
    correlation_function(tlm::TriangularLatticeModel, ψ0, sites, dt, tf, r0s::Tuple{Vararg{Int}}, combiner::Tuple{Int, Int}, dmrg_kwargs::NamedTuple; tdvp_kwargs::NamedTuple=(;))

Compute connected time-dependent spin correlations ⟨Sᵅ(r, t) Sᵝ(r₀, 0)⟩ - ⟨Sᵅ(r)⟩⟨Sᵝ(r₀)⟩ for every lattice site `r` relative to each reference site in `r0s`. The ground state is obtained via `dmrg_gs`, using `ψ0` and `sites` as the initial MPS and site set, and the excited state is evolved by TDVP from `t = 0` to `tf` in steps of `dt`. The pair `combiner = (α, β)` selects the spin operators `("S+", "S-", "Sz")[α]` and `("S+", "S-", "Sz")[β]`. Extra keyword arguments are forwarded to `dmrg` through `dmrg_kwargs` and to the TDVP evolution through `tdvp_kwargs`. The result is returned as a `CorrelationData` object storing an `nt × Ns × length(r0s)` array together with the lattice metadata.

Tips on the choice of `dt` and `tf`: Given a `max_energy` scale of interest, dt ≈ π / max_energy; given an `energy_resolution` scale of interest, tf ≈ π / energy_resolution.
"""
function correlation_function(tlm::TriangularLatticeModel, ψ0, sites, dt, tf, r0s::Tuple{Vararg{Int}}, combiner::Tuple{Int, Int}, dmrg_kwargs::NamedTuple; tdvp_kwargs::NamedTuple=(;))
    @assert all(1 ≤ r0 ≤ length(sites) for r0 in r0s) "All reference sites in r0s must be between 1 and $(length(sites))."
    @assert combiner in ((1,2), (2,1), (3,3)) "Combiner must be one of (1,2), (2,1), or (3,3)."

    (; Lx, Ly) = tlm

    ψ_gs, E_gs, sites, H_mpo = dmrg_gs(tlm, ψ0, sites; dmrg_kwargs...)
    @show "Maximum bond dimension in ground state MPS: " maxlinkdim(ψ_gs)

    Ns = length(sites)
    nt = Int(round(tf / dt)) + 1
    correlations = zeros(ComplexF64, nt, Ns, length(r0s))

    Threads.@threads for r0_index in eachindex(r0s)
        α, β = combiner
        r0 = r0s[r0_index]
        corrs = correlation_function_aux(r0, E_gs, ψ_gs, sites, H_mpo, dt, tf, α, β; tdvp_kwargs...)
        correlations[:, :, r0_index] = corrs
    end

    return CorrelationData(correlations, r0s, combiner, dt, tf, Lx, Ly)
end