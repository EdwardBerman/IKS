using FITSIO

include("../src/IKS.jl")

using .IKS

f = FITS("../data/forecast_lum_annular.fits")
x = read(f[2], "X_IMAGE_se")
y = read(f[2], "Y_IMAGE_se")
g1 = read(f[2], "g1_Rinv")
g2 = read(f[2], "g2_Rinv")


κ_E = IterativeKaisserSquires(g1, g2, x, y, 0.1)
