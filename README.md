# 2D FDTD Simulation of Room Acoustics with Mode Analysis

## Overview
This project implements a **2D Finite Difference Time Domain (FDTD)** simulation for room acoustics using MATLAB. The simulation models sound propagation in a rectangular room, capturing reflections, absorption at boundaries, and resonant room modes. The results are compared with **theoretical room modes** for validation.

The simulation outputs a **transfer function** between a source and a receiver in the room and overlays theoretical room modes for comparison.

---

## Features
- 2D FDTD simulation of acoustic wave propagation.
- Excitation with a **Ricker wavelet**.
- Absorbing boundary conditions using **impedance-based reflection**.
- Frequency domain analysis using **FFT**.
- Comparison with **theoretical room mode frequencies**.
- Visualization of room transfer function and resonant modes.

---

## Simulation Parameters
- **Room dimensions:** 4 m × 3 m
- **Spatial step (dh):** 0.01 m
- **Time step (dt):** Determined using Courant condition
- **Medium:** Air (density 1.21 kg/m³, speed of sound 341 m/s)
- **Absorption coefficient:** 0.1 (modifiable)
- **Excitation:** Ricker wavelet, central frequency 1000 Hz
- **Total time steps:** 10,000

---

## Usage
1. Open `FDTD_Room_Simulation.m` in MATLAB.
2. Adjust simulation parameters (room size, absorption coefficient, wavelet frequency) if needed.
3. Run the script.
4. The script will:
   - Simulate pressure and particle velocity fields over time.
   - Record pressure at a receiver point.
   - Compute FFT and transfer function.
   - Overlay theoretical room modes.
   - Save a figure `FDTD_Room_Transfer_Function.jpg`.

---

## Visualization
The output figure displays:
- **Blue line:** FDTD-simulated room transfer function in dB.
- **Red dashed lines:** Theoretical room modes for a rigid-walled room.
- **X-axis:** Frequency (Hz)
- **Y-axis:** Magnitude (dB)

---

## References
- Kinsler, L.E., Frey, A.R., Coppens, A.B., & Sanders, J.V. *Fundamentals of Acoustics*.
- Morse, P.M., & Ingard, K.U. *Theoretical Acoustics*.
- Taflove, A., & Hagness, S.C. *Computational Electrodynamics: The Finite-Difference Time-Domain Method* (adapted for acoustics).

---

## License
This project is open-source under the MIT License.
