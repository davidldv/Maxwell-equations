# Simulación 2D FDTD de las ecuaciones de Maxwell (modo TM: Ez, Hx, Hy)
# Implementa propagación electromagnética en 2D con PML y visualización avanzada

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LinearSegmentedColormap

# Constantes físicas
c0 = 299792458.0        # velocidad de la luz en el vacío (m/s)
mu0 = 4e-7 * np.pi      # permeabilidad del vacío
eps0 = 1.0 / (mu0 * c0**2)  # permitividad del vacío

# Parámetros numéricos de la malla 2D
nx, ny = 200, 150       # puntos espaciales (x, y)
dx = dy = 1e-3          # separación espacial (m)
S = 0.5                 # número de Courant 2D REDUCIDO para estabilidad
dt = S * dx / (c0 * np.sqrt(2))  # paso temporal (s)
nsteps = 1200           # MÁS pasos temporales para ver evolución completa
snap_every = 4          # dibujar cada N pasos

# Parámetros PML 2D
npml = 15               # grosor de las capas PML
sigma_max = 0.4         # conductividad máxima REDUCIDA para estabilidad
m = 3.0                 # orden del perfil de conductividad

def calculate_pml_profile_2d(nx, ny, npml, sigma_max, m):
    """Calcula perfiles de conductividad PML en 2D"""
    sigma_x = np.zeros((nx, ny))
    sigma_y = np.zeros((nx, ny))
    
    # PML en fronteras x (izquierda y derecha)
    for i in range(npml):
        # Izquierda
        d_norm = (npml - i) / npml
        sigma_val = sigma_max * (d_norm ** m)
        sigma_x[i, :] = sigma_val
        
        # Derecha
        sigma_x[nx-1-i, :] = sigma_val
    
    # PML en fronteras y (abajo y arriba)
    for j in range(npml):
        # Abajo
        d_norm = (npml - j) / npml
        sigma_val = sigma_max * (d_norm ** m)
        sigma_y[:, j] = sigma_val
        
        # Arriba
        sigma_y[:, ny-1-j] = sigma_val
    
    return sigma_x, sigma_y

# Crear perfiles de conductividad PML
sigma_x, sigma_y = calculate_pml_profile_2d(nx, ny, npml, sigma_max, m)

# Coeficientes PML 2D
def calculate_pml_coefficients(sigma_x, sigma_y, dt, eps0, mu0, dx, dy):
    """Calcula coeficientes PML para actualización FDTD 2D"""
    # Para campo eléctrico Ez
    da_x = np.exp(-sigma_x * dt / eps0)
    db_x = np.zeros_like(sigma_x)
    mask_x = sigma_x > 1e-12
    db_x[mask_x] = (1 - da_x[mask_x]) / (sigma_x[mask_x] * dx / eps0)
    db_x[~mask_x] = dt / (eps0 * dx)
    
    da_y = np.exp(-sigma_y * dt / eps0)
    db_y = np.zeros_like(sigma_y)
    mask_y = sigma_y > 1e-12
    db_y[mask_y] = (1 - da_y[mask_y]) / (sigma_y[mask_y] * dy / eps0)
    db_y[~mask_y] = dt / (eps0 * dy)
    
    # Para campos magnéticos Hx, Hy
    da_x_h = np.exp(-sigma_x * dt / mu0)
    db_x_h = np.zeros_like(sigma_x)
    db_x_h[mask_x] = (1 - da_x_h[mask_x]) / (sigma_x[mask_x] * dx / mu0)
    db_x_h[~mask_x] = dt / (mu0 * dx)
    
    da_y_h = np.exp(-sigma_y * dt / mu0)
    db_y_h = np.zeros_like(sigma_y)
    db_y_h[mask_y] = (1 - da_y_h[mask_y]) / (sigma_y[mask_y] * dy / mu0)
    db_y_h[~mask_y] = dt / (mu0 * dy)
    
    return (da_x, db_x, da_y, db_y, da_x_h, db_x_h, da_y_h, db_y_h)

# Calcular coeficientes PML
pml_coeffs = calculate_pml_coefficients(sigma_x, sigma_y, dt, eps0, mu0, dx, dy)
da_x, db_x, da_y, db_y, da_x_h, db_x_h, da_y_h, db_y_h = pml_coeffs

# Campos electromagnéticos 2D (modo TM)
Ez = np.zeros((nx, ny))     # Campo eléctrico Ez (perpendicular al plano xy)
Hx = np.zeros((nx, ny-1))   # Campo magnético Hx (en plano xy)
Hy = np.zeros((nx-1, ny))   # Campo magnético Hy (en plano xy)

# Variables auxiliares PML 2D
Psi_Ez_x = np.zeros((nx, ny))    # ∂Hy/∂x para Ez
Psi_Ez_y = np.zeros((nx, ny))    # ∂Hx/∂y para Ez
Psi_Hx_y = np.zeros((nx, ny-1))  # ∂Ez/∂y para Hx
Psi_Hy_x = np.zeros((nx-1, ny))  # ∂Ez/∂x para Hy

# Parámetros de la fuente 2D
source_x, source_y = nx//3, ny//2    # posición de la fuente
t0 = 30.0
spread = 8.0
freq = 2e9              # frecuencia 2 GHz

# Crear grillas de coordenadas para visualización
x = np.arange(nx) * dx * 1000  # convertir a mm
y = np.arange(ny) * dy * 1000  # convertir a mm
X, Y = np.meshgrid(x, y, indexing='ij')

# Configurar visualización
fig = plt.figure(figsize=(14, 10))
gs = fig.add_gridspec(2, 2, height_ratios=[3, 1], width_ratios=[3, 1])

# Gráfica principal: campo Ez en 2D
ax_main = fig.add_subplot(gs[0, 0])
im = ax_main.imshow(Ez.T, extent=[0, x[-1], 0, y[-1]], 
                   vmin=-1, vmax=1, cmap='RdBu_r', origin='lower')
ax_main.set_xlabel('x (mm)')
ax_main.set_ylabel('y (mm)')
ax_main.set_title('Campo Eléctrico Ez - Simulación 2D FDTD')

# Marcar regiones PML
pml_x = npml * dx * 1000
pml_y = npml * dy * 1000
ax_main.axvspan(0, pml_x, alpha=0.2, color='red', label='PML')
ax_main.axvspan(x[-1]-pml_x, x[-1], alpha=0.2, color='red')
ax_main.axhspan(0, pml_y, alpha=0.2, color='red')
ax_main.axhspan(y[-1]-pml_y, y[-1], alpha=0.2, color='red')

# Marcar fuente
ax_main.plot(source_x*dx*1000, source_y*dy*1000, 'wo', 
            markersize=8, markeredgecolor='black', linewidth=2, label='Fuente')

plt.colorbar(im, ax=ax_main, label='Ez (arb. units)')

# Gráfica de perfil PML en X
ax_pml_x = fig.add_subplot(gs[0, 1])
ax_pml_x.plot(sigma_x[:, ny//2], x, 'g-', linewidth=2)
ax_pml_x.fill_betweenx(x, 0, sigma_x[:, ny//2], alpha=0.3, color='green')
ax_pml_x.set_ylabel('x (mm)')
ax_pml_x.set_xlabel('σ_x')
ax_pml_x.set_title('PML X')
ax_pml_x.grid(True, alpha=0.3)

# Gráfica de perfil PML en Y
ax_pml_y = fig.add_subplot(gs[1, 0])
ax_pml_y.plot(y, sigma_y[nx//2, :], 'g-', linewidth=2)
ax_pml_y.fill_between(y, 0, sigma_y[nx//2, :], alpha=0.3, color='green')
ax_pml_y.set_xlabel('y (mm)')
ax_pml_y.set_ylabel('σ_y')
ax_pml_y.set_title('PML Y')
ax_pml_y.grid(True, alpha=0.3)

# Información de la simulación
ax_info = fig.add_subplot(gs[1, 1])
ax_info.axis('off')
info_text = f"""📊 Parámetros 2D FDTD:
• Malla: {nx}×{ny} puntos
• Resolución: {dx*1000:.1f} mm
• PML: {npml} capas
• Courant: S = {S:.2f}
• Frecuencia: {freq/1e9:.1f} GHz
• Paso temporal: {dt*1e12:.1f} ps"""

ax_info.text(0.05, 0.95, info_text, transform=ax_info.transAxes, 
            fontsize=10, verticalalignment='top', fontfamily='monospace')

def fdtd_step_2d(n):
    """Un paso de actualización FDTD 2D con PML"""
    global Ez, Hx, Hy, Psi_Ez_x, Psi_Ez_y, Psi_Hx_y, Psi_Hy_x
    
    # 1. Actualizar campos magnéticos Hx, Hy
    # Hx: -∂Ez/∂y
    curl_ez_y = np.zeros((nx, ny-1))
    curl_ez_y[:, :] = -(Ez[:, 1:] - Ez[:, :-1]) / dy
    
    Psi_Hx_y = da_y[:, :-1] * Psi_Hx_y + db_y_h[:, :-1] * curl_ez_y
    Hx += dt/mu0 * curl_ez_y + dt/mu0 * Psi_Hx_y
    
    # Hy: ∂Ez/∂x
    curl_ez_x = np.zeros((nx-1, ny))
    curl_ez_x[:, :] = (Ez[1:, :] - Ez[:-1, :]) / dx
    
    Psi_Hy_x = da_x[:-1, :] * Psi_Hy_x + db_x_h[:-1, :] * curl_ez_x
    Hy += dt/mu0 * curl_ez_x + dt/mu0 * Psi_Hy_x
    
    # 2. Actualizar campo eléctrico Ez
    # Ez: ∂Hy/∂x - ∂Hx/∂y
    curl_h_x = np.zeros((nx, ny))
    curl_h_y = np.zeros((nx, ny))
    
    curl_h_x[1:-1, :] = (Hy[1:, :] - Hy[:-1, :]) / dx
    curl_h_y[:, 1:-1] = -(Hx[:, 1:] - Hx[:, :-1]) / dy
    
    Psi_Ez_x = da_x * Psi_Ez_x + db_x * curl_h_x
    Psi_Ez_y = da_y * Psi_Ez_y + db_y * curl_h_y
    
    Ez += dt/eps0 * (curl_h_x + curl_h_y) + dt/eps0 * (Psi_Ez_x + Psi_Ez_y)
    
    # 3. Fuente: pulso gaussiano modulado
    pulse = np.exp(-0.5 * ((t0 - n) / spread)**2) * np.sin(2*np.pi*freq*n*dt)
    Ez[source_x, source_y] += pulse

# Función de animación
step_counter = 0
def update_2d(frame):
    global step_counter, Ez
    
    # Ejecutar varios pasos FDTD
    for _ in range(snap_every):
        fdtd_step_2d(step_counter)
        step_counter += 1
    
    # Actualizar visualización
    im.set_array(Ez.T)
    ax_main.set_title(f'Campo Ez - Paso {step_counter}/{nsteps} - t = {step_counter*dt*1e9:.1f} ns')
    
    return [im]

# Crear y ejecutar animación
nframes = nsteps // snap_every
ani = animation.FuncAnimation(fig, update_2d, frames=nframes, blit=False, interval=50)

print(f"🚀 Simulación 2D FDTD iniciada:")
print(f"   • Resolución: {nx}×{ny} = {nx*ny:,} puntos")
print(f"   • Memoria: ~{3*nx*ny*8/1e6:.1f} MB")
print(f"   • PML: {npml} capas ({npml*dx*1000:.1f} mm)")
print(f"   • Frecuencia: {freq/1e9:.1f} GHz")
print(f"   • Frames: {nframes}")

plt.tight_layout()
plt.show()