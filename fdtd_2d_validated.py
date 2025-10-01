# Simulación FDTD 2D con validación física - versión estable
# Incluye: medición de velocidad, materiales simples, monitoreo de energía

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time

# Constantes físicas
c0 = 2.998e8
mu0 = 4e-7 * np.pi
eps0 = 1.0 / (mu0 * c0**2)

# Parámetros optimizados para estabilidad
nx, ny = 100, 60
dx = dy = 3e-3          # 3 mm de resolución
S = 0.3                 # Courant muy conservador
dt = S * dx / (c0 * np.sqrt(2))
nsteps = 300
snap_every = 2

print(f"📊 Simulación FDTD con validación física:")
print(f"   • Malla: {nx}×{ny}")
print(f"   • dt = {dt*1e12:.1f} ps")
print(f"   • c₀ teórica = {c0/1e8:.3f} × 10⁸ m/s")

# Crear materiales simples
eps_r = np.ones((nx, ny))
eps_r[nx//2-8:nx//2+8, ny//2-6:ny//2+6] = 2.25  # Vidrio: εᵣ = 2.25

# PML simple
npml = 8
sigma_max = 0.15
sigma = np.zeros((nx, ny))

for i in range(npml):
    val = sigma_max * (i / npml) ** 2
    sigma[i, :] = val
    sigma[nx-1-i, :] = val
    sigma[:, i] = np.maximum(sigma[:, i], val)
    sigma[:, ny-1-i] = np.maximum(sigma[:, ny-1-i], val)

# Campos
Ez = np.zeros((nx, ny))
Hx = np.zeros((nx, ny-1))
Hy = np.zeros((nx-1, ny))

# Fuente y detectores
src_x, src_y = 20, ny//2
freq = 8e8  # 800 MHz (más baja para estabilidad)
t0 = 15
spread = 5

# Detectores para medir velocidad
detectors = [(35, ny//2), (50, ny//2), (65, ny//2)]
detector_data = [[] for _ in detectors]
detector_times = []

# Arrays para energía
energy_total = []
energy_times = []

# Configuración de visualización
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
ax_main, ax_energy, ax_detectors, ax_info = axes.flat

# Panel principal
x_mm = np.arange(nx) * dx * 1000
y_mm = np.arange(ny) * dy * 1000

# Mostrar materiales como fondo
im_bg = ax_main.imshow(eps_r.T, extent=[0, x_mm[-1], 0, y_mm[-1]], 
                      cmap='Greys', alpha=0.4, origin='lower', vmin=1, vmax=2.5)

# Campo Ez
im_ez = ax_main.imshow(Ez.T, extent=[0, x_mm[-1], 0, y_mm[-1]], 
                      vmin=-0.1, vmax=0.1, cmap='RdBu_r', origin='lower', alpha=0.8)

ax_main.set_xlabel('x (mm)')
ax_main.set_ylabel('y (mm)')
ax_main.set_title('Campo Ez + Materiales')

# Marcar fuente y detectores
ax_main.plot(src_x*dx*1000, src_y*dy*1000, 'yo', markersize=8, label='Fuente')
for i, (dx_pos, dy_pos) in enumerate(detectors):
    ax_main.plot(dx_pos*dx*1000, dy_pos*dy*1000, 'rs', markersize=5, label=f'Det{i+1}')
ax_main.legend()

plt.colorbar(im_ez, ax=ax_main, label='Ez')

# Configurar otros paneles
ax_energy.set_title('Energía Total vs Tiempo')
ax_energy.set_xlabel('Tiempo (ns)')
ax_energy.set_ylabel('Energía')
ax_energy.grid(True, alpha=0.3)

ax_detectors.set_title('Señales de Detectores')
ax_detectors.set_xlabel('Tiempo (ns)')
ax_detectors.set_ylabel('Ez')
ax_detectors.grid(True, alpha=0.3)

ax_info.axis('off')

def calculate_simple_energy(Ez, Hx, Hy):
    """Calcular energía de forma estable"""
    energy_e = 0.5 * eps0 * np.sum(Ez**2) * dx * dy
    
    # Interpolar H campos de forma simple
    H_total = 0
    if Hx.size > 0:
        H_total += np.sum(Hx**2)
    if Hy.size > 0:
        H_total += np.sum(Hy**2)
    
    energy_h = 0.5 * mu0 * H_total * dx * dy
    return energy_e + energy_h

def fdtd_step_stable(n):
    """Paso FDTD estable"""
    global Ez, Hx, Hy
    
    # Actualizar H (magnético)
    if Hx.shape[1] > 0:
        dEz_dy = (Ez[:, 1:] - Ez[:, :-1]) / dy
        Hx -= (dt / mu0) * dEz_dy
    
    if Hy.shape[0] > 0:
        dEz_dx = (Ez[1:, :] - Ez[:-1, :]) / dx
        Hy += (dt / mu0) * dEz_dx
    
    # Actualizar E (eléctrico)
    # Interior
    for i in range(1, nx-1):
        for j in range(1, ny-1):
            curl_h = 0
            if i < Hy.shape[0]:
                curl_h += (Hy[i, j] - Hy[i-1, j]) / dx
            if j < Hx.shape[1]:
                curl_h -= (Hx[i, j] - Hx[i, j-1]) / dy
            
            # Actualización con material
            eps_local = eps_r[i, j] * eps0
            Ez[i, j] += (dt / eps_local) * curl_h
    
    # Aplicar PML simple
    Ez *= (1 - sigma * dt / (2 * eps0))
    
    # Fuente
    pulse = np.exp(-0.5 * ((t0 - n) / spread)**2) * np.sin(2*np.pi*freq*n*dt)
    Ez[src_x, src_y] += 0.01 * pulse
    
    # Recopilar datos
    current_time = n * dt * 1e9
    detector_times.append(current_time)
    
    for i, (det_x, det_y) in enumerate(detectors):
        detector_data[i].append(Ez[det_x, det_y])
    
    # Energía
    total_energy = calculate_simple_energy(Ez, Hx, Hy)
    energy_total.append(total_energy)
    energy_times.append(current_time)

def estimate_wave_speed():
    """Estimar velocidad usando tiempo de picos"""
    if len(detector_times) < 50:
        return None
    
    # Encontrar máximos en cada detector
    peak_times = []
    for i, data in enumerate(detector_data):
        if len(data) > 20:
            # Buscar el pico principal
            abs_data = np.abs(data)
            if np.max(abs_data) > 0.001:  # Umbral mínimo
                peak_idx = np.argmax(abs_data[10:]) + 10  # Evitar ruido inicial
                peak_time = detector_times[peak_idx]
                peak_times.append(peak_time)
    
    if len(peak_times) >= 2:
        # Calcular velocidad entre primer y último detector
        distance = (detectors[-1][0] - detectors[0][0]) * dx
        time_diff = (peak_times[-1] - peak_times[0]) * 1e-9  # convertir a segundos
        if time_diff > 0:
            velocity = distance / time_diff
            return velocity
    
    return None

# Variables de animación
step = 0
measured_velocity = None

def animate(frame):
    global step, measured_velocity
    
    # Ejecutar pasos FDTD
    for _ in range(snap_every):
        fdtd_step_stable(step)
        step += 1
    
    # Actualizar campo
    im_ez.set_array(Ez.T)
    ax_main.set_title(f'Ez + Materiales - t = {step*dt*1e9:.2f} ns')
    
    # Actualizar energía
    if len(energy_times) > 1:
        ax_energy.clear()
        ax_energy.plot(energy_times, energy_total, 'b-', linewidth=2)
        ax_energy.set_title('Energía Total vs Tiempo')
        ax_energy.set_xlabel('Tiempo (ns)')
        ax_energy.set_ylabel('Energía')
        ax_energy.grid(True, alpha=0.3)
    
    # Actualizar detectores
    if len(detector_times) > 1:
        ax_detectors.clear()
        colors = ['red', 'green', 'blue']
        for i, (color, data) in enumerate(zip(colors, detector_data)):
            ax_detectors.plot(detector_times, data, color=color, 
                            linewidth=2, label=f'Detector {i+1}')
        ax_detectors.set_title('Señales de Detectores')
        ax_detectors.set_xlabel('Tiempo (ns)')
        ax_detectors.set_ylabel('Ez')
        ax_detectors.legend()
        ax_detectors.grid(True, alpha=0.3)
    
    # Medir velocidad
    if step % 30 == 0:
        measured_velocity = estimate_wave_speed()
    
    # Información
    ax_info.clear()
    ax_info.axis('off')
    
    speed_text = f"{measured_velocity/1e8:.3f}" if measured_velocity else "Calculando..."
    error_text = f"{abs(measured_velocity - c0)/c0*100:.1f}" if measured_velocity else "N/A"
    current_energy = energy_total[-1] if energy_total else 0
    max_energy = max(energy_total) if energy_total else 1
    energy_loss = (1 - current_energy/max_energy)*100 if max_energy > 0 else 0
    
    info_text = f"""VALIDACION FISICA FDTD 2D

VELOCIDAD DE ONDA:
  Teorica: {c0/1e8:.3f} × 10⁸ m/s
  Medida:  {speed_text} × 10⁸ m/s
  Error:   {error_text}%

MATERIALES:
  Vacio:  εᵣ = 1.0 (blanco)
  Vidrio: εᵣ = 2.25 (gris)
  c_vidrio = {c0/1.5/1e8:.3f} × 10⁸ m/s

ENERGIA:
  Actual: {current_energy:.2e} J/m
  Maxima: {max_energy:.2e} J/m
  Perdida: {energy_loss:.1f}% (PML)

ESTABILIDAD:
  Paso: {step}/{nsteps}
  Courant: S = {S:.2f} ✓
  Frecuencia: {freq/1e9:.1f} GHz
"""
    
    ax_info.text(0.05, 0.95, info_text, transform=ax_info.transAxes,
                fontsize=9, verticalalignment='top', fontfamily='monospace')
    
    return [im_ez]

# Ejecutar simulación
nframes = nsteps // snap_every
ani = animation.FuncAnimation(fig, animate, frames=nframes, 
                            interval=100, blit=False, repeat=False)

plt.tight_layout()
plt.show()

print(f"Simulación completada. Velocity final: {measured_velocity/1e8 if measured_velocity else 'N/A'} × 10⁸ m/s")