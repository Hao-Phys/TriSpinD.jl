"""
    dssf(corr_data::CorrelationData, σt, σr; measure::Symbol=:sperp, order::Int=1)
Compute the dynamical structure factor from time-dependent correlation data. The time smoothing parameter is `σt`, and the spatial smoothing parameter is `σr`. The `measure` argument specifies the type of intensity measure to compute, which can be `:sperp`, `:trace`, or `:component`.
"""
function dssf(corr_data::CorrelationData, σt, measure::Symbol; order::Int=1)
    (; r0s, correlations, dt, tf, Lx, Ly) = corr_data

    ts = range(0, tf; length=Int(round(tf / dt)) + 1)
    ts_full = [reverse(-ts[2:end]); ts]
    energies_full = fftshift(fftfreq(length(ts_full))) * 2π / dt

    Nsites = Lx * Ly
    intensities = zeros(Float64, length(energies_full), Nsites)

    q_points = available_q_points(Lx, Ly)

    for r0_index in eachindex(r0s)
        intensities += smooth_and_fourier_corr(correlations[:, :, r0_index], r0s[r0_index], Lx, Ly, dt, tf, σt; order)
    end

    return IntensityData(intensities, q_points, energies_full, Lx, Ly, measure)
end

function intensities_Sxx(Spm::IntensityData)
    (; q_points, energies_full, Lx, Ly, intensities) = Spm
    Sxx_intensities = 0.25 * intensities
    return IntensityData(Sxx_intensities, q_points, energies_full, Lx, Ly, :Sxx)
end

function intensities_Sxx(Spm::IntensityData, Smp::IntensityData)
    Spm_intensities = Spm.intensities
    Smp_intensities = Smp.intensities
    (; q_points, energies_full, Lx, Ly) = Spm
    Sxx_intensities = 0.25 * (Spm_intensities + Smp_intensities)
    return IntensityData(Sxx_intensities, q_points, energies_full, Lx, Ly, :Sxx)
end

function intensities(Spm::IntensityData, Szz::IntensityData, measure::Symbol)
    @assert measure in (:sperp, :trace) "Unsupported measure: $measure. Supported measures are :sperp, :trace"
    Sxx_intensities = intensities_Sxx(Spm).intensities
    Szz_intensities = Szz.intensities

    (; q_points, energies_full, Lx, Ly) = Spm
    if measure == :sperp
        Nsites = Lx * Ly
        intensities = zeros(Float64, length(energies_full), Nsites)
        for (iq, q_point) in enumerate(q_points)
            q_norm = norm(q_point)
            if q_norm ≈ 0
                intensities[:, iq] += Sxx_intensities[:, iq] + Szz_intensities[:, iq]
            else
                intensities[:, iq] += (2 - q_point[1]^2/q_norm^2 - q_point[2]^2/q_norm^2) * Sxx_intensities[:, iq] + Szz_intensities[:, iq]
            end
        end
        return IntensityData(intensities, q_points, energies_full, Lx, Ly, measure)
    elseif measure == :trace
        intensities = Sxx_intensities + Szz_intensities
        return IntensityData(intensities, q_points, energies_full, Lx, Ly, measure)
    end
end

function intensities(Spm::IntensityData, Smp::IntensityData, Szz::IntensityData, measure::Symbol)
    @assert measure in (:sperp, :trace) "Unsupported measure: $measure. Supported measures are :sperp, :trace"
    Sxx_intensities = intensities_Sxx(Spm, Smp).intensities
    Szz_intensities = Szz.intensities

    (; q_points, energies_full, Lx, Ly) = Spm
    if measure == :sperp
        Nsites = Lx * Ly
        intensities = zeros(Float64, length(energies_full), Nsites)
        for (iq, q_point) in enumerate(q_points)
            q_norm = norm(q_point)
            if q_norm ≈ 0
                intensities[:, iq] += Sxx_intensities[:, iq] + Szz_intensities[:, iq]
            else
                intensities[:, iq] += (2 - q_point[1]^2/q_norm^2 - q_point[2]^2/q_norm^2) * Sxx_intensities[:, iq] + Szz_intensities[:, iq]
            end
        end
        return IntensityData(intensities, q_points, energies_full, Lx, Ly, measure)
    elseif measure == :trace
        intensities = Sxx_intensities + Szz_intensities
        return IntensityData(intensities, q_points, energies_full, Lx, Ly, measure)
    end
end