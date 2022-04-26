# Integrators should take the ODE 'dydt',
# an array 'y' representing the function at time 't',
# the time 't', and the time step 'dt'
# They will return an array of the function at time (t+dt).

"""
	rungekutta4(dydt, y, t, dt)

Standard Fourth order Runge-Kutta integrator.
"""
function rungekutta4(dydt, y, t, dt)
    k1 = dt .* dydt(y, t)
    k2 = dt .* dydt(y .+ k1./2, t + dt/2)
    k3 = dt .* dydt(y .+ k2./2, t + dt/2)
    k4 = dt .* dydt(y .+ k3   , t + dt)

    return @. y + (k1 + 2k2 + 2k3 + k4) / 6 # y(t + dt)
end

"""
	rungekutta10_4(dydt, y, t, dt)

Fourth order Runge-Kutta with 10 steps from [1].

[1] : Highly efficient strong stability-preserving runge–kutta methods with low-storage
implementations - David I. Ketcheson (2008)
"""
function rungekutta10_4(dydt, y, t, dt)

    u1 = copy(y)
    u2 = copy(y)
    for _ in 1:5
        F = dydt(u1, t)
        u1 = @. u1 + F * dt/6
    end
    u2 = @. (u2 + 9u1)/25
    u1 = @. 15u2 - 5u1
    for _ in 6:9
        F = dydt(u1, t)
        u1 = @. u1 + F * dt/6
    end
    F = dydt(u1, t)
    u1 = @. u2 + 3u1/5 + F * dt/10

    return u1
end


"""
	euler(dydt, y, t, dt)

First order Euler integrator
"""
function euler(dydt, y, t, dt)
	return y .+ dt .* dydt(y, t) # y(t + dt)
end
