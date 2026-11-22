#!/usr/bin/env julia

# ============================================================
# Temporal Variance Spectrum (TVS)
# Fullerton et al. (1996)
#
# INPUT:
#   wave        :: Vector{Float64}
#   fluxes      :: Matrix{Float64}  (Nspectra × Nwave) (NaN values are not allowed!)
#   wave_range  :: Tuple(minλ, maxλ)
#
# OUTPUT:
#   ASCII .dat file:
#       wavelength   TVS   sqrtTVS
#
# ============================================================

using NaNStatistics
using DelimitedFiles

# ------------------------------------------------------------
# Estimate continuum noise σ_ic
# ------------------------------------------------------------
function estimate_continuum_noise(
    flux::Vector{Float64},
    continuum_mask=nothing
)

    if continuum_mask === nothing
        continuum_flux = flux
    else
        continuum_flux = flux[continuum_mask]
    end

    σ = nanstd(continuum_flux)

    # avoid division by zero
    return max(σ, 1e-12)
end


# ------------------------------------------------------------
# Compute TVS
# ------------------------------------------------------------
function compute_tvs(
    wave::Vector{Float64},
    fluxes::Matrix{Float64};
    wave_range::Tuple{Float64,Float64},
    continuum_mask=nothing,
    output_file::String="TVS.dat"
)

    if any(isnan, fluxes)
        error("NaN value found in input spectrum")
    end

    nspec, nwave = size(fluxes)

    if length(wave) != nwave
        error("wave and flux dimensions do not match")
    end

    # --------------------------------------------------------
    # Select wavelength range
    # --------------------------------------------------------
    wmin, wmax = wave_range

    mask = (wave .>= wmin) .& (wave .<= wmax)

    wave_sel = wave[mask]
    flux_sel = fluxes[:, mask]

    nsel = length(wave_sel)

    # --------------------------------------------------------
    # Estimate σ_ic for each spectrum
    # --------------------------------------------------------
    σ_ic = zeros(Float64, nspec)

    for i in 1:nspec
        σ_ic[i] = estimate_continuum_noise(
            vec(fluxes[i, :]),
            continuum_mask
        )
    end

    # --------------------------------------------------------
    # Fullerton weights
    #
    # σ0² = [ (1/N) Σ (1/σ_ic²) ]⁻¹
    # w_i = (σ0 / σ_ic)²
    # --------------------------------------------------------
    σ0² = 1.0 / nanmean(1.0 ./ (σ_ic .^ 2))
    σ0 = sqrt(σ0²)

    weights = (σ0 ./ σ_ic) .^ 2

    # --------------------------------------------------------
    # Weighted mean spectrum
    # --------------------------------------------------------
    mean_spec = zeros(Float64, nsel)

    weight_sum = sum(weights)

    for j in 1:nsel
        mean_spec[j] = sum(weights .* flux_sel[:, j]) / weight_sum
    end

    # --------------------------------------------------------
    # Compute TVS
    #
    # High-S/N approximation:
    # α_ij ≈ S_ij
    # --------------------------------------------------------
    tvs = zeros(Float64, nsel)

    for i in 1:nspec

        Sij = flux_sel[i, :]

        αij = max.(Sij, 1e-12)

        tvs .+= weights[i] .* ((Sij .- mean_spec) .^ 2) ./ αij
    end

    tvs ./= (nspec - 1)

    sqrt_tvs = sqrt.(tvs)

    # --------------------------------------------------------
    # Save output
    # --------------------------------------------------------
    output = hcat(wave_sel, tvs, sqrt_tvs)

    header = """
    # TVS computed following Fullerton et al. (1996)
    # Columns:
    # wavelength   TVS   sqrt(TVS)
    """
    
    if output_file != ""
    
        open(output_file, "w") do io

            write(io, header)

            writedlm(io, output)
        end
        println("TVS saved to: $output_file")
    
    end


    return wave_sel, tvs, sqrt_tvs
end


# ============================================================
# Example usage
# ============================================================
if abspath(PROGRAM_FILE) == @__FILE__

    # --------------------------------------------------------
    # Synthetic example
    # --------------------------------------------------------

    nspec = 10
    nwave = 5000

    wave = collect(range(6500.0, 6600.0, length=nwave))

    spectra = zeros(Float64, nspec, nwave)

    for i in 1:nspec

        flux = ones(Float64, nwave)

        # Variable Gaussian absorption line
        amp = 0.2 + 0.02 * randn()

        line = 0.2 .* exp.(-(wave .- 6562.8).^2 ./ (2 * 0.5^2))

        flux .-= line

        # Add noise
        flux .+= 0.005 .* randn(nwave)

        spectra[i, :] .= flux
    end

    # Example continuum regions
    continuum_mask =
        ((wave .> 6505.0) .& (wave .< 6520.0)) .|
        ((wave .> 6580.0) .& (wave .< 6595.0))

    # --------------------------------------------------------
    # Compute TVS
    # --------------------------------------------------------
    compute_tvs(
        wave,
        spectra;
        wave_range=(6540.0, 6585.0),
        continuum_mask=continuum_mask,
        output_file="Halpha_TVS.dat"
    )

print(a)

end