# Bessel Filter Design & Simulation Software

This software provides a comprehensive environment for designing and simulating Bessel filters of various types (Low-pass, High-pass, Band-pass, Band-stop) in both analog and digital domains. The application features graphical visualization of filter characteristics including magnitude/phase responses, pole-zero plots, and time-domain responses.

## Overview
- Supports 4 filter types: Low-pass, High-pass, Band-pass, Band-stop
- Works in both analog and digital domains
- Automatic order determination or manual specification
- Visualizations:
  - Magnitude response (linear/dB scale)
  - Phase response and group delay
  - Pole-zero plots with unit circle for digital filters
  - Step and impulse responses
- Multiple coefficient representations:
  - Transfer function
  - State-space
  - Second-order sections (SOS)

## Mathematical Models

### Bessel Polynomials
The software uses Bessel polynomials to design the analog low-pass prototype. The recursive formula for Bessel polynomials is:
```
B₀(s) = 1
B₁(s) = s + 1
Bₙ(s) = (2n-1)Bₙ₋₁(s) + s²Bₙ₋₂(s)
```
### Filter Transformations
1. **Low-pass Prototype Normalization**:
   - Poles are scaled to achieve unity cutoff frequency: `s → s/ωc`
2. **Frequency Transformations**:
   - **Low-pass to High-pass**: `s → ωc/s`
   - **Low-pass to Band-pass**: `s → (s² + ω0²)/(Bs)`
   - **Low-pass to Band-stop**: `s → Bs/(s² + ω0²)`
   where ω0 is center frequency and B is bandwidth
3. **Bilinear Transform (Analog to Digital)**:
   `s → (2/T)(z-1)/(z+1)`
   where T is sampling period

### Pole-Zero Analysis
- Poles are roots of the Bessel polynomial
- Zeros are determined by filter type:
  - Low-pass: No finite zeros
  - High-pass: Multiple zeros at origin
  - Band-pass: Complex conjugate zeros at center frequency
  - Band-stop: Zeros at origin and complex frequencies

### Response Calculations
1. **Magnitude Response**:
   `|H(jω)| = 1/|D(jω)|`
   where `D(jω)` is denominator polynomial evaluated at `s=jω`
2. **Phase Response**:
   `φ(ω) = -arg(D(jω))` (in degrees)
3. **Group Delay**:
   `τ(ω) = -dφ/dω ≈ -Δφ/Δω` (numerical differentiation)
4. **Time Domain Responses**:
   - Step response: Output when input `u(t) = 1`
   - Impulse response: Output when input `δ(t) = 1` at `t=0`

### Screenshots
![](https://github.com/KMORaza/Bessel_Filter_Design_and_Simulation_Software/blob/main/Bessel%20Filter%20Design%20and%20Simulation%20Software/screenshots/screenshot.png)
