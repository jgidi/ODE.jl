using ODE # Installed from https://github.com/jgidi/ODE.jl
using Plots, ProgressMeter

# Define the space
Nx = 256
xmin, xmax = 0.0, 2pi
x = range(xmin, step=(xmax-xmin)/Nx, length=Nx)

# Temporal discretization
dt = 3e-3
Nt = round(Int, 1.5/dt)

# Initial conditions
u = sin.(x)

# Differential equation
dudt(u, t) = -0.5derivative(u.^2, x) + 1e-3*derivative(u, x, order=2)

# Integrate initial conditions with ODE.integrate_ode
time, ut = integrate_ode(dudt, u, Nt, dt, integrator=rungekutta4)

# Make animation with only some instants
instants_to_plot = 1:Nt

progressbar = Progress(length(instants_to_plot), desc="Animating...")
anim = @animate for i in instants_to_plot
    p = plot()

    # Plot initial conditions
    plot!(p, x, ut[:, 1], l = (:gray, :dash), label="\$u(x, 0)\$")
    # Plot ut at current time
    plot!(p, x, ut[:, i], l = (:black), label="\$u(x, t)\$")
    # Add other labels
    plot!(p,
          xlabel = "\$ x\$",
          ylabel = "\$ u\$",
          title = "t = $(round(time[i], digits=2))",
          )
end

# Save animation as mp4
mp4(anim, "burgers.mp4")
