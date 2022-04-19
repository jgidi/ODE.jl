using ODE
using Plots, ProgressMeter

# Definimos variable espacial
Nx = 256
xmin, xmax = 0.0, 2pi
x = range(xmin, step=(xmax-xmin)/Nx, length=Nx)

# Condiciones iniciales
u = sin.(x)

# Discretizcion temporal
Nt = 1000
dt = 1e-3

# Ecuación diferencial
v(x) = 1.0
dudt(u, t) = - v(x) * derivative(u, x)

time, ut = integrate_ode(dudt, u, Nt, dt)

progressbar = Progress(Nt, desc="Animating...")
anim = @animate for i in 1:Nt
    p = plot()

    plot!(p, x, ut[:, 1], l = (:gray, :dash), label="\$u(x, 0)\$")
    plot!(p, x, ut[:, i], l = (:black), label="\$u(x, t)\$")
    plot!(p, title = "t = $(round(time[i], digits=3))")
end

mp4(anim, "traslacion.mp4")
