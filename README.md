# Simulaciones FDTD de las Ecuaciones de Maxwell

[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://python.org)
[![NumPy](https://img.shields.io/badge/NumPy-1.20+-orange.svg)](https://numpy.org)
[![Matplotlib](https://img.shields.io/badge/Matplotlib-3.0+-green.svg)](https://matplotlib.org)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## 📋 Descripción

Este proyecto implementa simulaciones **FDTD (Finite-Difference Time-Domain)** para resolver las ecuaciones de Maxwell en 1D, 2D y esquemas 3D. Incluye técnicas avanzadas como **PML (Perfectly Matched Layers)** para absorción de fronteras, materiales heterogéneos, y validación física cuantitativa.

### ⚡ Características principales

- **Simulaciones 1D, 2D y 3D** de propagación electromagnética
- **PML (Perfectly Matched Layers)** para absorción sin reflexiones
- **Materiales heterogéneos** con diferentes permitividades
- **Validación física** con medición de velocidad de onda
- **Monitoreo de energía** para verificar estabilidad numérica
- **Visualización en tiempo real** con animaciones interactivas
- **Análisis cuantitativo** de efectividad de técnicas numéricas

## 🚀 Instalación

### Requisitos del sistema
- Python 3.8 o superior
- Windows, macOS o Linux

### Instalación rápida

```bash
# Clonar el repositorio
git clone https://github.com/davidldv/Maxwell-equations.git
cd Maxwell-equations

# Crear entorno virtual
python -m venv .venv

# Activar entorno virtual
# Windows:
.venv\Scripts\activate
# macOS/Linux:
source .venv/bin/activate

# Instalar dependencias
pip install numpy matplotlib
```

## 📁 Estructura del proyecto

```
ecuaciones_maxwell/
├── README.md                     # Este archivo
├── .gitignore                    # Archivos a ignorar en Git
├── main.py                       # Simulación FDTD 1D con PML
├── fdtd_2d_stable.py             # Simulación FDTD 2D estabilizada
├── fdtd_2d_validated.py          # FDTD 2D con validación física completa
├── fdtd_2d_advanced.py           # FDTD 2D con materiales y monitoreo
├── fdtd_3d_skeleton.py           # Esqueleto para FDTD 3D
├── pml_validation_final.py       # Validación cuantitativa de PML
└── pml_analysis_results.txt      # Resultados de análisis PML
```

## 🎯 Guía de uso

### 1. Simulación básica 1D

```bash
python main.py
```

**Características:**
- Propagación de pulso gaussiano en 1D
- PML (Perfectly Matched Layers) para absorción
- Visualización de campo Ez vs tiempo
- Perfil de conductividad PML

### 2. Simulación 2D estabilizada

```bash
python fdtd_2d_stable.py
```

**Características:**
- Ondas electromagnéticas circulares en 2D
- Modo TM (Ez, Hx, Hy)
- Parámetros optimizados para estabilidad
- Visualización de campo 2D con código de colores

### 3. Simulación 2D con validación física

```bash
python fdtd_2d_validated.py
```

**Características:**
- **Medición de velocidad de onda** con 3 detectores
- **Materiales heterogéneos** (vacío + vidrio)
- **Monitoreo de energía** electromagnética total
- **Validación en tiempo real** de parámetros físicos

### 4. Análisis de efectividad PML

```bash
python pml_validation_final.py
```

**Características:**
- **Comparación cuantitativa** CON vs SIN PML
- **Métricas objetivas** de coeficiente de reflexión
- **Análisis de absorción** de energía
- **Factor de mejora** calculado

## 🔬 Fundamentos científicos

### Ecuaciones de Maxwell implementadas

El proyecto resuelve las ecuaciones de Maxwell en forma diferencial:

**Modo TM (2D):**
```
∂Hx/∂t = -(1/μ₀) ∂Ez/∂y
∂Hy/∂t = (1/μ₀) ∂Ez/∂x  
∂Ez/∂t = (1/ε₀εᵣ)(∂Hy/∂x - ∂Hx/∂y)
```

**Criterio de estabilidad de Courant:**
```
S = c₀ × dt / dx ≤ 1/√(dimensiones)
```

### PML (Perfectly Matched Layers)

Las PML implementan absorción perfecta mediante:

```python
σ(d) = σ_max × (d/δ)^m
```

Donde:
- `σ_max`: conductividad máxima
- `d`: distancia desde la interfaz PML
- `δ`: grosor total de la PML
- `m`: orden del perfil (típicamente 2-4)

## 📊 Ejemplos de resultados

### Propagación 1D
- **Velocidad medida:** c₀ = 2.998 × 10⁸ m/s
- **Error típico:** < 1%
- **Absorción PML:** > 99%

### Simulación 2D
- **Patrones de onda:** Circulares desde fuente puntual
- **Materiales:** Reflexión/refracción en interfaces
- **Estabilidad:** Energía conservada en región libre

### Validación PML
- **Coeficiente de reflexión:** R < 0.01
- **Factor de mejora:** 100x vs fronteras rígidas
- **Absorción de energía:** 95%+ efectiva

## ⚙️ Parámetros de configuración

### Parámetros físicos principales

```python
# Constantes fundamentales
c0 = 2.998e8        # Velocidad de la luz (m/s)
mu0 = 4π × 10⁻⁷     # Permeabilidad del vacío
eps0 = 8.854e-12    # Permitividad del vacío

# Parámetros de malla
nx, ny = 120, 80    # Puntos espaciales
dx = dy = 2e-3      # Resolución espacial (mm)
dt = 2.36e-12       # Paso temporal (ps)
```

### Configuración PML

```python
npml = 10           # Grosor de capas PML
sigma_max = 0.3     # Conductividad máxima
m = 2.0            # Orden del perfil
```

### Materiales

```python
eps_r = 1.0         # Vacío
eps_r = 2.25        # Vidrio (n = 1.5)
eps_r = 4.0         # Dieléctrico típico
```

## 🎨 Visualización

### Mapas de campo 2D
- **Código de colores:** Rojo/Azul para Ez positivo/negativo
- **Escala de grises:** Materiales con εᵣ > 1
- **Regiones PML:** Marcadas en rojo transparente

### Gráficas temporales
- **Detectores múltiples:** Para medición de velocidad
- **Monitoreo de energía:** Conservación y disipación
- **Análisis espectral:** FFT de señales temporales

### Métricas en tiempo real
- **Velocidad de onda:** Comparación con c₀
- **Error porcentual:** Precisión numérica
- **Estabilidad:** Factores de Courant y energía

## 🧪 Validación y pruebas

### Tests físicos implementados

1. **Velocidad de propagación**
   ```python
   v_medida = distancia / tiempo_vuelo
   error = |v_medida - c₀| / c₀ × 100%
   ```

2. **Conservación de energía**
   ```python
   E_total = ½ε|E|² + ½μ|H|²
   ```

3. **Efectividad PML**
   ```python
   R = |E_reflejado| / |E_incidente|
   ```

### Criterios de aceptación
- **Error de velocidad:** < 5%
- **Estabilidad numérica:** Sin overflow
- **Absorción PML:** R < 0.1
- **Conservación energía:** Pérdida < 1% en región libre

## 🔧 Personalización avanzada

### Añadir nuevos materiales

```python
def create_custom_material(nx, ny):
    eps_r = np.ones((nx, ny))
    # Región con εᵣ personalizada
    eps_r[x1:x2, y1:y2] = epsilon_value
    return eps_r
```

### Fuentes personalizadas

```python
def custom_source(n, freq, amplitude):
    return amplitude * np.sin(2*np.pi*freq*n*dt) * envelope(n)
```

### PML con perfiles avanzados

```python
def advanced_pml_profile(d, delta):
    return sigma_max * np.exp(-(delta-d)/decay_length)
```

## 📚 Referencias científicas

1. **Taflove, A. & Hagness, S.C.** - *Computational Electrodynamics: The Finite-Difference Time-Domain Method* (3rd Ed.)
2. **Berenger, J.P.** - *A perfectly matched layer for the absorption of electromagnetic waves* (1994)
3. **Yee, K.S.** - *Numerical solution of initial boundary value problems* (1966)
4. **Gedney, S.D.** - *An anisotropic PML absorbing media for FDTD simulation* (1996)

## 🤝 Contribuciones

¡Las contribuciones son bienvenidas! Por favor:

1. Fork el proyecto
2. Crea una rama para tu feature (`git checkout -b feature/AmazingFeature`)
3. Commit tus cambios (`git commit -m 'Add AmazingFeature'`)
4. Push a la rama (`git push origin feature/AmazingFeature`)
5. Abre un Pull Request

### Áreas de contribución prioritarias
- **FDTD 3D completo** con visualización volumétrica
- **Materiales dispersivos** (Debye, Drude, Lorentz)
- **Fuentes complejas** (antenas, guías de onda)
- **Optimización GPU** con CUDA/OpenCL
- **Interfaz gráfica** para configuración interactiva

## 📄 Licencia

Este proyecto está bajo la Licencia MIT. Ver `LICENSE` para más detalles.

## 👨‍💻 Autor

**David** - [davidldv](https://github.com/davidldv)

## 🙏 Agradecimientos

- Comunidad científica de electromagnetismo computacional
- Desarrolladores de NumPy y Matplotlib
- Referencias bibliográficas de FDTD clásico
- Algoritmos PML de Berenger y Gedney

---

## 📞 Soporte

¿Preguntas o problemas? 

- 📧 Email: [crear issue en GitHub](../../issues)
- 📖 Wiki: [Documentación extendida](../../wiki)
- 💬 Discusiones: [GitHub Discussions](../../discussions)

---

*Simulando el universo electromagnético, una ecuación de Maxwell a la vez* ⚡📐