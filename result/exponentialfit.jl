using Plots
using LsqFit

model(t, p) = p[1] .+ p[2]*t 

tdata = range(0,10,20)
ydata = model(tdata, [1.0 2.0]) + 0.01*randn(length(tdata))
p0 = [0.5, 0.5]
fits = curve_fit(model, tdata, ydata, p0)
param = fits.param
fits
param = fits.param
scatter(tdata, ydata, label="Data")

plot!(tdata, fits(tdata), label="Exponential Fit")

# Define the exponential model
function exponential_model(x, a, b)
    return a * exp.(b * x)
end

# Generate some sample data
x = 1:10
y = [2.718, 7.389, 20.086, 54.598, 148.413, 403.429, 1096.633, 2980.958, 8103.084, 22026.466]

# Fit the curve to the data
fit = curve_fit(exponential_model, x, y)

# Plot the data points and the fitted curve
scatter(x, y, label="Data")
plot!(x, fit(x), label="Exponential Fit")