using ODE # Installed from https://github.com/jgidi/ODE.jl
using Plots, ProgressMeter

# Define the space
Nx = 256
xmin, xmax = 0.0, 2pi
x = range(xmin, step=(xmax-xmin)/Nx, length=Nx)

# Temporal discretization
Nt = 1000
dt = (x[2]-x[1])/4

# Initial conditions
u = @. exp(-100(x-1)^2)

# Differential equation
c(x) = @. 0.2 + sin(x -1)^2
dudt(u, t) = - c(x) .* derivative(u, x)

# Integrate initial conditions with ODE.integrate_ode
time, ut = integrate_ode(dudt, u, Nt, dt, integrator=rungekutta4)

# Make animation with only some instants
instants_to_plot = 1:10:Nt

progressbar = Progress(length(instants_to_plot), desc="Animating...")
anim = @animate for i in instants_to_plot
    p = plot()

    # Plot initial conditions
    plot!(p, x, ut[:, 1], l = (:gray, :dash), label="\$u(x, 0)\$")
    # Plot ut at current time
    plot!(p, x, ut[:, i], l = (:black), label="\$u(x, t)\$")
    # Add title with current time
    plot!(p, title = "t = $(round(time[i], digits=3))")
end

# Save animation as mp4
mp4(anim, "trefethen.mp4")
