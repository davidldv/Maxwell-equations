# Simulación 2D FDTD estabilizada - Maxwell en 2D con PML
# Versión optimizada para evitar overflow numérico

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Constantes físicas
c0 = 2.998e8            # velocidad de la luz (m/s)
mu0 = 4e-7 * np.pi      # permeabilidad del vacío
eps0 = 1.0 / (mu0 * c0**2)  # permitividad del vacío

# Parámetros numéricos optimizados
nx, ny = 120, 80        # malla más pequeña para estabilidad
dx = dy = 2e-3          # resolución espacial (2 mm)
S = 0.5                 # Courant más conservador
dt = S * dx / (c0 * np.sqrt(2))  # paso temporal estable
nsteps = 400            # pasos temporales
snap_every = 2          # frames

# Parámetros PML conservadores
npml = 10               # PML más delgada
sigma_max = 0.3         # conductividad reducida
m = 2.0                 # perfil cuadrático

print(f"📊 Parámetros estabilizados:")
print(f"   • Malla: {nx}×{ny} puntos")
print(f"   • dt = {dt*1e12:.2f} ps")
print(f"   • Courant S = {S:.2f}")

def create_pml_2d(nx, ny, npml, sigma_max, m):
    """Crear PML 2D con perfil suave"""
    sigma_x = np.zeros((nx, ny))
    sigma_y = np.zeros((nx, ny))
    
    # PML gradual en X
    for i in range(npml):
        profile = sigma_max * ((npml - i - 1) / npml) ** m
        sigma_x[i, :] = profile          # izquierda
        sigma_x[nx-1-i, :] = profile     # derecha
    
    # PML gradual en Y
    for j in range(npml):
        profile = sigma_max * ((npml - j - 1) / npml) ** m
        sigma_y[:, j] = profile          # abajo
        sigma_y[:, ny-1-j] = profile     # arriba
    
    return sigma_x, sigma_y

# Crear conductividades PML
sigma_x, sigma_y = create_pml_2d(nx, ny, npml, sigma_max, m)

# Coeficientes PML simplificados
ca = np.exp(-sigma_x * dt / eps0)
cb = (1 - ca) / (sigma_x + 1e-15)  # evitar división por cero
ca_y = np.exp(-sigma_y * dt / eps0)
cb_y = (1 - ca_y) / (sigma_y + 1e-15)

# Campos electromagnéticos
Ez = np.zeros((nx, ny), dtype=np.float64)
Hx = np.zeros((nx, ny-1), dtype=np.float64)
Hy = np.zeros((nx-1, ny), dtype=np.float64)

# Variables PML auxiliares
Psi_Ez_x = np.zeros((nx, ny), dtype=np.float64)
Psi_Ez_y = np.zeros((nx, ny), dtype=np.float64)

# Fuente puntual
src_x, src_y = nx//3, ny//2
freq = 1e9              # 1 GHz
t0 = 20.0
spread = 6.0

# Configurar visualización
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 8))

# Mapa de campo Ez
x_mm = np.arange(nx) * dx * 1000
y_mm = np.arange(ny) * dy * 1000
im1 = ax1.imshow(Ez.T, extent=[0, x_mm[-1], 0, y_mm[-1]], 
                vmin=-0.5, vmax=0.5, cmap='RdBu_r', origin='lower')
ax1.set_title('Campo Eléctrico Ez')
ax1.set_xlabel('x (mm)')
ax1.set_ylabel('y (mm)')

# Marcar fuente
ax1.plot(src_x*dx*1000, src_y*dy*1000, 'yo', markersize=6, markeredgecolor='k')

# PML en X
ax2.plot(sigma_x[:, ny//2], x_mm, 'g-', linewidth=2)
ax2.fill_betweenx(x_mm, 0, sigma_x[:, ny//2], alpha=0.3, color='green')
ax2.set_title('PML Perfil X')
ax2.set_xlabel('σ_x')
ax2.set_ylabel('x (mm)')
ax2.grid(True, alpha=0.3)

# PML en Y
ax3.plot(y_mm, sigma_y[nx//2, :], 'g-', linewidth=2)
ax3.fill_between(y_mm, 0, sigma_y[nx//2, :], alpha=0.3, color='green')
ax3.set_title('PML Perfil Y')
ax3.set_xlabel('y (mm)')
ax3.set_ylabel('σ_y')
ax3.grid(True, alpha=0.3)

# Información
ax4.axis('off')
info = f"""Simulación 2D FDTD
Malla: {nx} × {ny}
Resolución: {dx*1000:.1f} mm
PML: {npml} capas
Frecuencia: {freq/1e9:.1f} GHz
Courant: {S:.2f}
"""
ax4.text(0.1, 0.7, info, fontsize=10, verticalalignment='top', fontfamily='monospace')

plt.colorbar(im1, ax=ax1, shrink=0.8)

def fdtd_update(n):
    """Actualización FDTD 2D estabilizada"""
    global Ez, Hx, Hy, Psi_Ez_x, Psi_Ez_y
    
    # Actualizar Hx: -∂Ez/∂y
    dEz_dy = np.diff(Ez, axis=1) / dy
    Hx += (dt / mu0) * (-dEz_dy)
    
    # Actualizar Hy: ∂Ez/∂x
    dEz_dx = np.diff(Ez, axis=0) / dx
    Hy += (dt / mu0) * dEz_dx
    
    # Actualizar Ez con PML
    dHy_dx = np.zeros((nx, ny))
    dHx_dy = np.zeros((nx, ny))
    
    dHy_dx[1:-1, :] = np.diff(Hy, axis=0) / dx
    dHx_dy[:, 1:-1] = -np.diff(Hx, axis=1) / dy
    
    curl_h = dHy_dx + dHx_dy
    
    # Aplicar PML de forma estabilizada
    Psi_Ez_x = ca * Psi_Ez_x + cb * dHy_dx * (dt / eps0)
    Psi_Ez_y = ca_y * Psi_Ez_y + cb_y * dHx_dy * (dt / eps0)
    
    Ez += (dt / eps0) * curl_h + Psi_Ez_x + Psi_Ez_y
    
    # Fuente: pulso sinusoidal gaussiano
    pulse = np.exp(-0.5 * ((t0 - n) / spread)**2) * np.sin(2*np.pi*freq*n*dt)
    if 0 <= src_x < nx and 0 <= src_y < ny:
        Ez[src_x, src_y] += 0.1 * pulse  # amplitud reducida

# Función de animación
step = 0
def animate(frame):
    global step
    
    # Ejecutar pasos FDTD
    for _ in range(snap_every):
        fdtd_update(step)
        step += 1
    
    # Actualizar visualización
    im1.set_array(Ez.T)
    ax1.set_title(f'Ez - Paso {step}/{nsteps} - t = {step*dt*1e9:.1f} ns')
    
    # Estadísticas
    max_field = np.max(np.abs(Ez))
    ax4.clear()
    ax4.axis('off')
    stats = f"""Estado actual:
Paso: {step}/{nsteps}
Tiempo: {step*dt*1e9:.2f} ns
|Ez|_max: {max_field:.3f}
Memoria: {Ez.nbytes/1024:.1f} KB

PML σ_max: {sigma_max:.2f}
Grosor: {npml} puntos
"""
    ax4.text(0.05, 0.95, stats, fontsize=9, transform=ax4.transAxes,
            verticalalignment='top', fontfamily='monospace')
    
    return [im1]

# Ejecutar animación
nframes = nsteps // snap_every
ani = animation.FuncAnimation(fig, animate, frames=nframes, 
                            interval=80, blit=False, repeat=True)

print(f"🚀 Iniciando simulación 2D...")
print(f"   Memoria total: ~{(Ez.nbytes + Hx.nbytes + Hy.nbytes)/1024:.1f} KB")

plt.tight_layout()
plt.show()