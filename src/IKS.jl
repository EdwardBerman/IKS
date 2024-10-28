module IterativeKaisserSquires

using LinearAlgebra
using WCS

#=
α , δ are astrometric coordinates right ascension and declination
x, y are pixel coordinates from source extractor tangent plane projection (with potential corrections)
g1, g2 are shear
M is the Mask
=#

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
        resolution::Float64)::Tuple{AbstractMatrix{<:Real}, AbstractMatrix{<:Real}}

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

function IterativeKaisserSquires(g1::AbstractVector{<:Real}, 
        g2::AbstractVector{<:Real}, 
        x::AbstractVector{<:Real},
        y::AbstractVector{<:Real},
        resolution::Float64,
        sharpness::Float64,
        max_iters::Int64=100)::Tuple{AbstractMatrix{<:Real}, AbstractMatrix{<:Real}}

    @assert size(g1) == size(g2) "g1 and g2 must have the same size"
    @assert size(x) == size(y) "x and y must have the same size"

    # to start, we take γ ≈ g and use iterative methods to improve the estimate
    γ1, γ2, M = average_shear_binning(g1, g2, x, y, resolution)
    γ1, γ2, M = add_zero_padding(γ1), add_zero_padding(γ2), add_zero_padding(M)

    γ = γ1 + γ2*im

    λ_min = 0.0
    k = 0
end

function IterativeKaisserSquires(g1::AbstractVector{<:Real}, 
        g2::AbstractVector{<:Real}, 
        worldcoords::AbstractMatrix{<:Real},
        wcs::WCSTransform, 
        resolution::Float64,
        sharpness::Float64,
        max_iters::Int64=1000)::Tuple{AbstractMatrix{<:Real}, AbstractMatrix{<:Real}}

    @assert size(g1) == size(g2) "g1 and g2 must have the same size"
    pixcoords = world_to_pix(wcs, worldcoords)
    x = pixcoords[1,:]
    y = pixcoords[2,:]

end

end 
