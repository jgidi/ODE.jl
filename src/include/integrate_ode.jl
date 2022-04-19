"""
    integrate_ode(dydt, y0, Nt, dt; integrator=euler)

Integrates the ODE `dydt`, starting from the array with initial conditions 'y0'
during `Nt` iterations with a time step `dt`, using the integrator `integrator`.
Returns a tuple `(time, trajectory)`, where `time` is the time axis ranging from `dt` to `dt*Nt` with steps `dt`,
and `trajectory` is an accumulator. `trajectory[i, j]` is the `i-th` point of `y` at the `j-th` time.

Notes
====

* By default, uses the simplest integrator, `euler`.
"""
function integrate_ode(dydt, y0, Nt, dt; integrator=euler)

    # Time axis
    time = range(dt, step=dt, length=Nt)

    # Define an accumulator for the i-th value at the j-th time
    trajectory = Array{Float64}(undef, length(y0), Nt)

    y = copy(y0)
    progressbar = Progress(Nt)
    for (i, t) in enumerate(time)

        # Obtain y(t+dt) from y(t) and dydt(y, t)
        y = integrator(dydt, y, t, dt)

        # Save the i-th value of y, y(t=i*dt)
        trajectory[:, i] .= y

        # Update progressbar
        next!(progressbar)
    end

    return time, trajectory
end
