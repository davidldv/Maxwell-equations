# Simulación 1D FDTD de las ecuaciones de Maxwell (Ez, Hy) — ejecución visible
# Produzco una animación en tiempo real de la propagación de un pulso gaussiano.
# Requisitos: numpy, matplotlib (ya disponibles en el entorno).
# No usé seaborn ni colores explícitos (regla del entorno).

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Constantes físicas
c0 = 299792458.0        # velocidad de la luz en el vacío (m/s)
mu0 = 4e-7 * np.pi      # permeabilidad del vacío
eps0 = 1.0 / (mu0 * c0**2)  # permitividad del vacío

# Parámetros numéricos de la malla
nx = 400               # puntos espaciales
dx = 1e-3              # separación espacial (m)
S = 0.99               # número de Courant (debe ser <= 1 en 1D para estabilidad)
dt = S * dx / c0       # paso temporal (s)
nsteps = 1200          # pasos temporales totales
snap_every = 4         # dibujar cada N pasos (para acelerar animación)

# Campos (1D, Yee grid: Ez on integer points, Hy staggered)
Ez = np.zeros(nx)
Hy = np.zeros(nx-1)    # Hy definido entre puntos de Ez

# Parámetros de la fuente (pulso gaussiano en el centro)
t0 = 40.0
spread = 12.0
source_pos = nx // 4   # colocamos la fuente no en el centro para ver reflexión

# Máscara absorbente simple (región de "amortiguamiento" en los bordes)
n_absorb = 40
absorb_coeff = np.ones(nx)
# aplicamos una ventana gaussiana en los extremos para amortiguar
for i in range(n_absorb):
    coeff = np.exp(-((n_absorb - i) / (n_absorb/3))**2)
    absorb_coeff[i] = 1.0 - 0.95 * coeff    # izquierda
    absorb_coeff[-1 - i] = 1.0 - 0.95 * coeff  # derecha

# Para comparar velocidad: calcular distancia en m entre puntos
x = np.arange(nx) * dx

# Listas para la animación
frames = []

# Preparar figura
fig, ax = plt.subplots(figsize=(10,4))
line, = ax.plot(x, Ez)
ax.set_ylim(-1.2, 1.2)
ax.set_xlim(0, x[-1])
ax.set_xlabel("x (m)")
ax.set_ylabel("Ez (arb. units)")
ax.set_title("Simulación 1D FDTD — Pulso gaussiano")

# Bucle de tiempo (actualizaciones FDTD)
def fdtd_step(n):
    global Ez, Hy
    # Actualizar campo magnético Hy (entre i y i+1)
    # Hy[i]^{n+1/2} = Hy[i]^{n-1/2} + (dt/(mu*dx)) * (Ez[i+1]^n - Ez[i]^n)
    Hy += (dt / (mu0 * dx)) * (Ez[1:] - Ez[:-1])

    # Actualizar campo eléctrico Ez
    # Ez[i]^{n+1} = Ez[i]^n + (dt/(eps*dx)) * (Hy[i] - Hy[i-1])
    Ez[1:-1] += (dt / (eps0 * dx)) * (Hy[1:] - Hy[:-1])

    # Fuente: pulso gaussiano en el tiempo (forzado en Ez[source_pos])
    pulse = np.exp(-0.5 * ((t0 - n) / spread)**2)
    Ez[source_pos] += pulse

    # Aplicar máscara absorbente en los bordes multiplicando por coeficiente < 1
    Ez *= absorb_coeff

# Función de actualización de la animación (se llama cada frame)
step_counter = 0
def update(frame):
    global step_counter, Ez, Hy
    # Realizar varias sub-iteraciones para avanzar más entre frames
    for _ in range(snap_every):
        fdtd_step(step_counter)
        step_counter += 1
    line.set_ydata(Ez)
    ax.set_title(f"Simulación 1D FDTD — paso temporal {step_counter}/{nsteps}")
    return (line,)

# Crear animación
nframes = nsteps // snap_every
ani = animation.FuncAnimation(fig, update, frames=nframes, blit=True, interval=10)

# Mostrar animación en salida
plt.close(fig)  # evitar doble impresión estática
ani

