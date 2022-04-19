"""
    derivative(fx, x; order=1)

Derivada espectral de orden 'order' entero no negativo.

Notes
=====

* Con 'orden=N' nos referimos a la N-esima derivada, no a la presición del método.
* Por defecto, order = 1
"""
function derivative(fx, x; order=1)
    Nx = length(x)
    dx = x[2]-x[1]

    # Eje de frecuencias multiplicado por i = √(-1)
    ik = im .* rfftfreq(Nx, 2pi/dx)

    return irfft(rfft(fx) .* ik.^order, Nx)
end
