# Simulación 1D FDTD de las ecuaciones de Maxwell (Ez, Hy) con PML
# Implementa Perfectly Matched Layers para absorción realista

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

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

# Parámetros PML
npml = 20              # grosor de las capas PML (puntos)
sigma_max = 0.8        # conductividad máxima normalizada (reducida para estabilidad)
m = 3.0                # orden del perfil de conductividad (típicamente 2-4)
R0 = 1e-6              # reflexión deseada en la interfaz PML

# Calcular perfil de conductividad PML
def calculate_pml_profile(npml, sigma_max, m):
    """Calcula el perfil de conductividad PML con gradiente suave"""
    profile = np.zeros(npml)
    for i in range(npml):
        # Perfil polinomial: σ = σ_max * (d/δ)^m
        # donde d es la distancia desde la interfaz
        d_norm = (npml - i - 1) / npml  # normalizado de 0 a 1
        profile[i] = sigma_max * (d_norm ** m)
    return profile

# Crear perfiles de conductividad para ambos extremos
sigma_e = np.zeros(nx)  # conductividad eléctrica
sigma_h = np.zeros(nx-1)  # conductividad magnética

# PML izquierda
pml_left = calculate_pml_profile(npml, sigma_max, m)
sigma_e[:npml] = pml_left
sigma_h[:npml] = pml_left[:npml] if npml <= nx-1 else pml_left[:nx-1]

# PML derecha  
pml_right = calculate_pml_profile(npml, sigma_max, m)
sigma_e[-npml:] = pml_right[::-1]  # invertir para que aumente hacia el borde
if npml <= nx-1:
    sigma_h[-npml:] = pml_right[::-1]
else:
    sigma_h[-(nx-1):] = pml_right[:nx-1][::-1]

# Coeficientes PML para la actualización FDTD
# Basado en la formulación de Taflove & Hagness
da_e = np.exp(-sigma_e * dt / eps0)
db_e = np.zeros_like(sigma_e)
mask_e = sigma_e > 1e-12  # evitar división por cero
db_e[mask_e] = (1 - da_e[mask_e]) / (sigma_e[mask_e] * dx / eps0)
db_e[~mask_e] = dt / (eps0 * dx)  # valor estándar en región sin PML

da_h = np.exp(-sigma_h * dt / mu0)  
db_h = np.zeros_like(sigma_h)
mask_h = sigma_h > 1e-12  # evitar división por cero
db_h[mask_h] = (1 - da_h[mask_h]) / (sigma_h[mask_h] * dx / mu0)
db_h[~mask_h] = dt / (mu0 * dx)  # valor estándar en región sin PML

# Campos (1D, Yee grid: Ez on integer points, Hy staggered)
Ez = np.zeros(nx)
Hy = np.zeros(nx-1)

# Variables auxiliares para PML (almacenan historiales de campo)
Psi_Ez = np.zeros(nx)  # variable auxiliar para Ez
Psi_Hy = np.zeros(nx-1)  # variable auxiliar para Hy

# Parámetros de la fuente (pulso gaussiano)
t0 = 40.0
spread = 12.0
source_pos = nx // 4   # colocamos la fuente no en el centro para ver reflexión

# Para comparar velocidad: calcular distancia en m entre puntos
x = np.arange(nx) * dx

# Preparar figura con visualización mejorada para PML
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))

# Gráfica principal: campos electromagnéticos
line_ez, = ax1.plot(x, Ez, 'b-', label='Ez (Campo Eléctrico)', linewidth=1.5)
line_hy, = ax1.plot(x[:-1] + dx/2, Hy*377, 'r-', label='Hy × 377 (Campo Magnético)', linewidth=1.5)
ax1.set_ylim(-1.5, 1.5)
ax1.set_xlim(0, x[-1])
ax1.set_ylabel("Campo (arb. units)")
ax1.set_title("Simulación 1D FDTD con PML — Propagación de Pulso Gaussiano")
ax1.legend()
ax1.grid(True, alpha=0.3)

# Gráfica secundaria: perfil de conductividad PML
ax2.plot(x, sigma_e, 'g-', linewidth=2, label='Conductividad PML')
ax2.fill_between(x, 0, sigma_e, alpha=0.3, color='green')
ax2.set_xlim(0, x[-1])
ax2.set_xlabel("x (m)")
ax2.set_ylabel("σ (normalizada)")
ax2.set_title("Perfil de Conductividad PML")
ax2.legend()
ax2.grid(True, alpha=0.3)

# Marcar regiones PML
ax1.axvspan(0, npml*dx, alpha=0.2, color='red', label='PML Izquierda')
ax1.axvspan((nx-npml)*dx, x[-1], alpha=0.2, color='red', label='PML Derecha')

# Bucle de tiempo FDTD con PML
def fdtd_step_pml(n):
    global Ez, Hy, Psi_Ez, Psi_Hy
    
    # 1. Actualizar variable auxiliar Psi_Hy y campo magnético Hy
    curl_e = Ez[1:] - Ez[:-1]  # curl de Ez
    Psi_Hy = da_h * Psi_Hy + db_h * curl_e
    Hy = Hy + dt / (mu0 * dx) * curl_e - dt / mu0 * Psi_Hy

    # 2. Actualizar variable auxiliar Psi_Ez y campo eléctrico Ez
    curl_h = np.zeros(nx)
    curl_h[1:-1] = Hy[1:] - Hy[:-1]  # curl de Hy (interior)
    # Condiciones de frontera: Ez[0] y Ez[-1] se quedan como están
    
    Psi_Ez = da_e * Psi_Ez + db_e * curl_h
    Ez = Ez + dt / (eps0 * dx) * curl_h - dt / eps0 * Psi_Ez

    # 3. Fuente: pulso gaussiano en el tiempo
    pulse = np.exp(-0.5 * ((t0 - n) / spread)**2)
    Ez[source_pos] += pulse

# Función de actualización de la animación
step_counter = 0
def update(frame):
    global step_counter, Ez, Hy
    
    # Realizar varias iteraciones FDTD entre frames
    for _ in range(snap_every):
        fdtd_step_pml(step_counter)
        step_counter += 1
    
    # Actualizar gráficas
    line_ez.set_ydata(Ez)
    line_hy.set_ydata(Hy * 377)  # Factor 377 para escalar impedancia
    
    ax1.set_title(f"FDTD con PML — Paso {step_counter}/{nsteps} — "
                  f"Reflexión teórica: {R0:.2e}")
    
    return (line_ez, line_hy)

# Crear animación
nframes = nsteps // snap_every
ani = animation.FuncAnimation(fig, update, frames=nframes, blit=True, interval=50)

print(f"🔬 Simulación FDTD con PML iniciada:")
print(f"   • Grosor PML: {npml} puntos ({npml*dx*1000:.1f} mm)")
print(f"   • Conductividad máxima: {sigma_max:.2f}")
print(f"   • Reflexión esperada: {R0:.2e}")
print(f"   • Frames totales: {nframes}")

plt.tight_layout()
plt.show()

