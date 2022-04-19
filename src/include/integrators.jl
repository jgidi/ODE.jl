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
    k2 = dt .* dydt(@. y + k1/2, t + dt/2)
    k3 = dt .* dydt(@. y + k2/2, t + dt/2)
    k4 = dt .* dydt(@. y + k3  , t + dt)

    return @. y + (k1 + 2k2 + 2k3 + k4) / 6 # y(t + dt)
end

"""
	euler(dydt, y, t, dt)

First order Euler integrator
"""
function euler(dydt, y, t, dt)
	return y .+ dt .* dydt(y, t) # y(t + dt)
end
