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

@misc{goyder2023decay,
  author       = {Goyder, Hugh},
  title        = {Determining differential equations from measured 
                  decay responses for nonlinear systems},
  howpublished = {Presentation notes, Tribomechadynamics 2023},
  year         = {2023},
  institution  = {Cranfield University},
  note         = {Rice University, August 2023}
}

@misc{goyder2023forbidden,
  author       = {Goyder, Hugh},
  title        = {Exploring nonlinear systems on the complex plane 
                  by looking into the forbidden zone},
  howpublished = {Presentation notes, Tribomechadynamics 2023},
  year         = {2023},
  institution  = {Cranfield University},
  note         = {Rice University, August 2023}
}

@article{berrut2004barycentric,
  author  = {Berrut, Jean-Paul and Trefethen, Lloyd N.},
  title   = {Barycentric {L}agrange Interpolation},
  journal = {SIAM Review},
  volume  = {46},
  number  = {3},
  pages   = {501--517},
  year    = {2004},
  doi     = {10.1137/S0036144502417715}
}

@inproceedings{goyder2018extracting,
  author    = {Goyder, H. G. D. and Lancereau, D. P. T.},
  title     = {Extracting natural frequencies and damping from time 
               histories. Better than the frequency domain?},
  booktitle = {Proceedings of ISMA 2018},
  year      = {2018}
}

## License

MIT (or specify your preferred license)
