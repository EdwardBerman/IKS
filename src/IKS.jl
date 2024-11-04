module IterativeKaisserSquires

using LinearAlgebra
using WCS
using FFTW
using Flux
# using DeconvOptim
using SpecialFunctions

#=
αJ2000 , δJ2000 are astrometric coordinates right ascension and declination
x, y are pixel coordinates from source extractor tangent plane projection (with potential corrections)
g1, g2 are shear
k1, k2 are wave numbers
κ_E and κ_B are the convergence fields
M is the Mask
W is the wavelet transform
note also that α = ϕk in later contexts
=#

function γ_to_κ(γ1::AbstractMatrix{<:Real}, γ2::AbstractMatrix{<:Real})::AbstractMatrix{<:Complex}
    @assert size(γ1) == size(γ2) "γ1 and γ2 must have the same size"
    k1 = fftfreq(size(γ1, 1))
    k2 = fftfreq(size(γ1, 2))
    temp = zeros(size(k2, 1), size(k1, 1))
    k1 = zeros(size(k2, 1), size(k1, 1)) .+ k1'
    k1 = k1'
    for col in eachcol(temp)
        col .= col .+ k2  
    end
    k2 = temp 
    k2 = k2'
    k_squared = k1.^2 + k2.^2
    ϵ = 1e-10
    denominator = k_squared .+ ϵ
    numerator = (k1 .- k2*im) .^2
    κ = ifft((numerator ./ denominator) .* fft(γ1 .+ γ2*im))
    return κ
end

function λ_max(γ1::AbstractMatrix{<:Real}, γ2::AbstractMatrix{<:Real})::Float64
    κ = γ_to_κ(γ1, γ2)
    κ_E = real(κ)
    α = dct(κ_E) # TO-DO: Check the normalization on this. Note α = ϕ^T k
    λ_max = maximum(α)
    return λ_max
end

function add_zero_padding(image::AbstractMatrix{<:Real}, padding::Int64=2)::AbstractMatrix{<:Real}
    @assert padding >= 0 "padding must be a non-negative integer"
    padded_image = zeros(Int(size(image, 1)*padding), Int(size(image, 2)*padding))
    
    start_row = floor(Int, (padded_rows - original_rows) / 2) + 1
    start_col = floor(Int, (padded_cols - original_cols) / 2) + 1
    end_row = start_row + original_rows - 1
    end_col = start_col + original_cols - 1
    
    @assert start_row > 0 && start_col > 0 "Padding too small to center the image."
    @assert end_row <= padded_rows && end_col <= padded_cols "Padding too small to center the image."
    
    padded_image[start_row:end_row, start_col:end_col] .= image
    return padded_image
end

function average_shear_binning(g1::AbstractVector{<:Real}, 
        g2::AbstractVector{<:Real}, 
        x::AbstractVector{<:Real},
        y::AbstractVector{<:Real},
        resolution::Float64)::Tuple{AbstractMatrix{<:Real}, AbstractMatrix{<:Real}, AbstractMatrix{<:Real}}

    @assert size(g1) == size(g2) "g1 and g2 must have the same size"
    @assert size(x) == size(y) "x and y must have the same size"

    x_min = minimum(x)
    x_max = maximum(x)
    x_bins = Int(round((x_max - x_min) / resolution))

    y_min = minimum(y)
    y_max = maximum(y)
    y_bins = Int(round((y_max - y_min) / resolution))

    x_bin_indices = floor.((x .- x_min) / resolution) .+ 1
    y_bin_indices = floor.((y .- y_min) / resolution) .+ 1

    sum_g1 = zeros(x_bins, y_bins)
    sum_g2 = zeros(x_bins, y_bins)
    counts = zeros(x_bins, y_bins)
    M = zeros(x_bins, y_bins)

    for k in eachindex(x)
        i = x_bin_indices[k]
        j = y_bin_indices[k]
        sum_g1[i, j] += g1[k]
        sum_g2[i, j] += g2[k]
        counts[i, j] += 1
    end

    M = counts .> 0

    ϵ = 1e-10
    shear_field_g1 = sum_g1 ./ (counts .+ ϵ)
    shear_field_g2 = sum_g2 ./ (counts .+ ϵ)

    return shear_field_g1, shear_field_g2, M

end

function b3spline_smoothing(step::Int64=1)::AbstractMatrix{<:Complex}
    step_hole = Int(step)
    c1 = 1.0 / 16.0
    c2 = 1.0 / 4.0
    c3 = 3.0 / 8.0
    length = Int(4 * step_hole + 1)
    kernel1d = zeros(Float64, 1, length)
    kernel1d[1, 1] = c1
    kernel1d[1, end] = c1
    kernel1d[1, step_hole + 1] = c2
    kernel1d[1, end - step_hole] = c2
    kernel1d[1, 2 * step_hole + 1] = c3
    kernel2d = kernel1d' * kernel1d
    return kernel2d
end

function W(κ::AbstractMatrix{<:Complex}, scales::Int64)::AbstractMatrix{<:Real}
    wavelet_coefficients = zeros(scales, size(κ)...)
    image_in = κ
    image_out = zeros(size(κ))
    image_auxiliary = zeros(size(κ))

    normalization_image_input = zeros(size(κ))
    normalization_image_input[size(κ, 1) ÷ 2, size(κ, 2) ÷ 2] = 1.0
    normalization_coefficients = zeros(scales, size(κ)...)
    normalization_image_output = zeros(size(κ))

    kernel_size = size(b3spline_smoothing(step=1))  
    padding = ((kernel_size[1] - 1) ÷ 2, (kernel_size[2] - 1) ÷ 2)
    conv_layer = Conv(kernel_size, 1 => 1, stride=(1, 1), pad=padding)
    conv_layer.weight .= reshape(Float32.(kernel), kernel_size[1], kernel_size[2], 1, 1)
    
    for i in 1:(scales - 1)
        kernel = b3spline_smoothing(step=2^i)

        input_data = reshape(Float32.(image_in), size(image_in)..., 1, 1)
        image_out = reshape(conv_layer(input_data), size(κ))
        image_auxiliary = reshape(conv_layer(image_out), size(κ))
        wavelet_coefficients[i, :, :] .= (image_in .- image_auxiliary)
        image_in = image_auxiliary

        normalization_image_output = conv_layer(reshape(Float32.(normalization_image_input), size(κ)..., 1, 1))
        normalization_image_input = reshape(conv_layer(normalization_image_output), size(κ))
        squared_norm = normalization_image_output .^ 2
        normalization_coefficients[i, :, :] .= squared_norm
    end
    
    wavelet_coefficients[end, :, :] .= image_out
    normalization_coefficients[end, :, :] .= normalization_image_output .^ 2
    tabNs = zeros(scales, 1)

    for i in 1:scales
        tabNs[i] = sqrt(sum(sum(normalization_coefficients[i, :, :], dims=1), dims=2)[1])
    end
    
    normalized_wavelength_coefficients = wavelet_coefficients ./ tabNs
    
    return normalized_wavelength_coefficients
end

function WT(wavelet_coefficients::AbstractMatrix{<:Real})::AbstractMatrix{<:Real}

end

function normalize_wavelets(wavelet_coefficients::AbstractMatrix{<:Real}, M::AbstractMatrix{<:Real})::AbstractMatrix{<:Real}
    normalized_wavelets = zeros(size(wavelet_coefficients))
    for i in 1:size(wavelet_coefficients, 1)
        in_mask = M .> 0
        σ_in = std(wavelet_coefficients[i, :, :][Bool.(in_mask)])

        out_mask = M .<= 0
        σ_out = std(wavelet_coefficients[i, :, :][Bool.(out_mask)])

        normalized_wavelets[i, :, :] = wavelet_coefficients[i, :, :] .* (σ_out / σ_in)
    end
    return normalized_wavelets
end

function IterativeKaisserSquires(g1::AbstractVector{<:Real}, 
        g2::AbstractVector{<:Real}, 
        x::AbstractVector{<:Real},
        y::AbstractVector{<:Real},
        resolution::Float64,
        λmin::Float64=0.0,
        max_iters_inner::Int64=3,
        max_iters_outer::Int64=100)::Tuple{AbstractMatrix{<:Real}, AbstractMatrix{<:Real}}

    @assert size(g1) == size(g2) "g1 and g2 must have the same size"
    @assert size(x) == size(y) "x and y must have the same size"

    # to start, we take γ ≈ g and use iterative methods to improve the estimate
    γ1, γ2, M = average_shear_binning(g1, g2, x, y, resolution)
    γ1, γ2, M = add_zero_padding(γ1), add_zero_padding(γ2), add_zero_padding(M)

    λmin = 0.0
    λmax = λ_max(γ1, γ2)

    κ_E = zeros(size(γ1))
    γ_init = γ1 + γ2*im

    wavelet_scales = Int(round(log(minimum(size(κ_E)))))
    wavelet_coefficients = zeros(wavelet_scales, size(κ_E)...)

    for k in 1:max_iters_outer
        γ_k = γ_init .* (1 .- κ_E)
        κ_k = γ_to_κ(real(γ_k), imag(γ_k))
        κ_i = κ_k
        λ_i = λmax

        for i in 1:max_iters_inner
            α = dct(κ_i)  # α = ϕ^T k^i
            α_tilde = [abs(α[i]) > λ_i ? α[i] : 0 for i in 1:length(α)]
            κ_i = idct(α_tilde)
            wavelet_coefficients = normalize_wavelets(W(κ_i, wavelet_scales), M) 


            # TO-DO: Fill in this steps, the Starlet wave transform in particular, also double check when to use DCT and IDCT and how to normalize them

            λ_i = λmin + (λmax - λmin) * (1 - erf(2.8 * i / max_iters_inner))
        end
        κ_E = real(κ_i) 
    end
    return κ_E
end

function IterativeKaisserSquires(g1::AbstractVector{<:Real}, 
        g2::AbstractVector{<:Real}, 
        worldcoords::AbstractMatrix{<:Real},
        wcs::WCSTransform, 
        resolution::Float64,
        λmin::Float64=0.0,
        max_iters_inner::Int64=3,
        max_iters_outer::Int64=1000)::Tuple{AbstractMatrix{<:Real}, AbstractMatrix{<:Real}}

    @assert size(g1) == size(g2) "g1 and g2 must have the same size"
    pixcoords = world_to_pix(wcs, worldcoords)
    x = pixcoords[1,:]
    y = pixcoords[2,:]

end

end 
