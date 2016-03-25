PIPE
====

Software package for analysis of photo-conversion microscopy movies.

This Matlab code was developed by Rotem Gura Sadovsky as a contribution to the scientific community in the laboratory of Dr. Jeremy England at Massachusetts Institute of Technology. See license file for usage and distribution guidelines.

Photo-converted Intensity Profile Expansion (PIPE) measures protein diffusion in live cells. PIPE works by pulsing photo-convertible fluorescent proteins, generating a peaked fluorescence signal at the pulsed region, and analyzing the spatial expansion of the signal as diffusion spreads it out. The width of the expanding signal is directly related to the protein ensemble mean-square displacement, from which the diffusion coefficient of the ensemble is calculated.

____________
Executables:

PIPE.m -- main program. Input data either as a folder path or as a matlab matrix read by PIPE_read_2d_movie.m. Run from within Matlab: see examples in file header.

PIPE_read_2d_movie.m -- create a data matrix from raw data files. By default, users are not required to run this script, but they may do so for more advanced usage.

____________
Core of analysis:

PIPE_calc_gaussian_widths.m -- fitting each intensity profile to a Gaussian and extracting its width.

PIPE_calc_D_of_expanding_gaussians.m -- fit widths vs. time to a normal diffusion model.

____________

Data characterization and pre-processing:


PIPE_identify_cell.m -- distinguish cell from cell exterior in image

PIPE_find_min_distance_to_boundary.m -- calculate distance from photo-conversion pulse to cell boundary

PIPE_correct_baseline.m -- subtract fluorescence baseline from signal

PIPE_find_phoc_peak.m -- find maximum of noisy time series

PIPE_find_pulse_center.m -- find center of 2D peak

PIPE_detect_pulse_end.m -- detect the last frame where the photo-conversion pulse is still on

PIPE_define_1d_profile_mask.m -- define a straight line passing through a given point within a 2D image 

PIPE_apply_profile_mask.m -- get intensity profile from a 2D image

PIPE_estimate_width0.m -- estimate the width of the initial photo-converted protein ensemble

PIPE_calc_time_points_by_noise_analysis.m -- determine optimal averaging scheme to reduce noise

PIPE_reduce_n_time_points.m -- reduce number of analyzed time points to a user-provided number

PIPE_correct_immobile_fraction.m -- subtract late-times intensity profile from rest of signal

PIPE_smooth_with_gaussian_kernel.m -- customized 1D or 2D smoothing procedure

_________
