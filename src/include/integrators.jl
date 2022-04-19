# Definimos un integrador. Runge-Kutta de 4to orden en 'dt'
# Input: dydt(y, t), y(t), t, dt
# Output: y(t + dt)
function rungekutta4(dydt, y, t, dt)
    # El anteponer un punto a un operador le indica a Julia que
    # esa operacion debe ser hecha componente a componente.
    # Esto es equivalente a decir que "vectorializa" la operación
    k1 = dt .* dydt(y, t)
    k2 = dt .* dydt(y .+ k1./2, t + dt/2)
    k3 = dt .* dydt(y .+ k2./2, t + dt/2)
    k4 = dt .* dydt(y .+ k3, t + dt)

    # El macro @. vectorializa todas las operaciones que vengan despues de el
    return @. y + (k1 + 2k2 + 2k3 + k4) / 6 # y(t + dt)
end

# Definimos otro integrador sólo para mostrar que podemos proveerlo a 'integrate_edo'
# Integrador de 1er orden en 'dt'
function euler(dydt, y, t, dt)
    return y .+ dt .* dydt(y, t) # y(t + dt)
end
