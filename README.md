# barycentric-nonlinear-fitting
Barycentric rational fitting for nonlinear vibration decay analysis

# Barycentric Rational Fitting for Nonlinear Vibration Decay

MATLAB implementation of barycentric rational interpolation with greedy 
support-point selection, applied to experimental vibration decay data 
for nonlinear system identification.

## Background

This work follows the approach proposed by Hugh Goyder (Cranfield University) 
for identifying differential equations from measured decay responses of 
nonlinear structures. The core idea: the Fourier spectrum of a nonlinear 
decay can be approximated by a sum of linear basis functions, fitted via 
barycentric rational interpolation on the complex plane.

## Files

- `src/real_data_barycentric_fit.m` — main MATLAB script
- `results/barycentric_fit_result.png` — example fitting result (100 support points)

## Usage

1. Place your `.mat` data file in the working directory
2. Adjust `Fs`, `fc`, `f_min`, `f_max`, and `Nsupport_target` in the script
3. Run `real_data_barycentric_fit.m` in MATLAB

## Method

1. FFT of the decay signal → complex spectrum H(s) on s = iω
2. Initial support points: edges + peak of |H|
3. Iterative greedy selection:
   - Solve weights via SVD null-space (Loewner-style matrix)
   - Evaluate barycentric model
   - Add support point at max absolute error
4. Stop when reaching noise floor (~1/1000 of max |H|)

## References

- Goyder, H. (2023). *Determining differential equations from measured decay 
  responses for nonlinear systems.* Tribomechadynamics 2023, Rice University.
- Goyder, H. (2023). *Exploring nonlinear systems on the complex plane by 
  looking into the forbidden zone.* Tribomechadynamics 2023, Rice University.
- Berrut, J.-P., & Trefethen, L. N. (2004). Barycentric Lagrange Interpolation. 
  *SIAM Review*, 46(3), 501–517.

## License

MIT (or specify your preferred license)
