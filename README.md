# Guiado magnético autónomo de un cuadricóptero

Código de simulación en **MATLAB** del Trabajo Fin de Grado *«Guiado magnético autónomo
de un cuadricóptero siguiendo conductores de corriente»* — Grado en Física, Facultad de
Ciencias, Universidad de Granada (curso 2025/2026).

Un cuadricóptero de carga sigue una *«autopista magnética»* (conductores con corriente)
guiándose **únicamente por el campo magnético** —sin GPS ni trayectoria predefinida—,
midiendo el campo y su gradiente con un magnetómetro vectorial y una IMU. Las ganancias
del control PID se optimizan con **Cuckoo Search** sobre una integración **RK4**.

## Características

- Campo magnético por **Biot–Savart** (fórmula analítica cerrada de segmentos finitos).
- Estimación de la desviación por **gradiente** del campo (frente al método del ángulo).
- Distribución de conductores en **escalera con esquina** (raíles antiparalelos + travesaños).
- **Guiado autónomo**: rumbo por mínimo gradiente, centrado anulando la componente
  perpendicular, control de velocidad de avance con freno en curva y altura por la escala
  vertical del campo (`L_z`); INS acotada por el campo.
- Sensores reales modelados: magnetómetro **PNI RM3100** e IMU **ST LSM6DSO**.

## Estructura del repositorio

| Archivo | Descripción |
|---|---|
| `dron_escalera_autoguiado.m` | Guiado autónomo puro (programa principal). |
| `dron_escalera3D.m` | Seguimiento con corrección 3D sobre distribución en escalera. |
| `dron_L_3D.m` | Corrección por gradiente en 3 ejes sobre distribución en L. |
| `dron_L_1D_grad_vs_ang.m` | Comparación de corrección por gradiente vs por ángulo (1D). |
| `campos3.m` | Comparación de métodos de cálculo del campo (analítico / numérico / tabla). |
| `ejecutar_todos.m` | Lanzador secuencial de los programas. |

## Requisitos

- MATLAB R2019b o superior.
- *Parallel Computing Toolbox* (opcional): la optimización usa `parfor`; sin el toolbox,
  `parfor` se ejecuta como un bucle normal.

## Uso

```matlab
% Ejecutar un programa concreto:
dron_escalera_autoguiado

% O lanzarlos en secuencia:
ejecutar_todos
```
Cada programa crea su propia carpeta en `resultados/` con las figuras (PNG a 300 dpi).

## Autor y licencia

José Rosa Girona — Universidad de Granada.
Código publicado bajo licencia MIT (ver `LICENSE`).
