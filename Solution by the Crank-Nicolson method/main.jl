u = incompressibleCouetteFlow(1)
y = linspace(0, 1, length(u))
plt = plot(u,
        y,
        color_palette = :Greens,
        legend = false,
        yticks = 0:0.2:1,
        linewidth = 2,
        title = "Velocity profiles for different number of timesteps",
        xlabel = L"$u$",
        ylabel = L"$\frac{y}{D}$")
display(plt)
for i in 1:5:100
    u = incompressibleCouetteFlow(i)
    plot!(u, y, legend = false)
    display(plt)
end
# saves the figure as a pdf
# savefig("testBilde.pdf")
