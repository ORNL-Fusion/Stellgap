# Stellgap
Stellgap calculates the shear Alfvén gap structure for 3D configurations (stellarators, RFPs, 3D tokamaks)

These codes are used to calculate shear Alfven continua for 3D configurations, both with and without sound wave coupling effects. The associated paper is D. A. Spong, R. Sanchez, A. Weller, "Shear Alfvén continua in stellarators," Phys. Plasmas 10 (2003) 3217–3224.
The workflow is as follows:

(a) prepare a VMEC equilibrium for the case of interest

(b) Run xbooz_xform (from Stellopt code suite) to convert from VMEC coordinates to Boozer coordinates. This run needs the in_booz.* file, which tells it what range of m/n modes to use and the selected surfaces to use from the VMEC run. Typically, the first and last of the VMEC surfaces are removed, because they can sometimes have noisy data.

(c) Using the boozmn.* file produced in (b), run xmetric_ver* to extract needed data and calculate the metric elements used in the continuum calculation. This produces a file called tae_data_boozer which is used as input for xstgap.

(d) The continuum calculation is done by xstgap*. There are two versions: one with sound wave couplings (xstgap_snd_ver*) and one without (xstgap). In addition to ae_data_boozer, these use the input files fourier.dat and plasma.dat. Fourier.dat specifies the set of fourier modes used to represent the continuum eigenmodes and plasma.dat contains information/profiles for the the plasma. The first line of fourier.dat gives the field periods and the surface grid parameters ith and izt used for calculating theta/zeta dependent coefficients of the continuum equation. It is important that these values for ith and izt are exactly the same as were used in the metric_element_create.f code.

(e) After running xstgap, the post-processing code (either post_process.f or post_process_snd.f) should be run. This produces a text file called alfven_post which contains columns of radial coordinate, frequency, dominant m and dominant n. Typically, the frequency vs. radius is plotted as a scatter plot. The code stelgp_to_silo.f is provided to convert data from alfven_post to a Visit silo file. This allows plotting the continua using the Visit software with the option to color code the points with either the dominant m or n.


Compilation scripts are provided for the different codes in the files whose name begins with “bld”. These will need to be edited, depending on what compiler is being used and where the needed libraries are located. Stellgap is constructed so that it can be compiled either as a serial version or a parallel version using precompilation flags. When running in parallel the number of surfaces requested will be divided by the number of processors and each group of surfaces allocated to a different processor. At the end of the run, all processors will write their results out to separate files, with names containing the processor number. These must be concatenated together in the post-processing step (e.g., as done in the post_process_snd.f code).

