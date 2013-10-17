## impro

This repository contains an assortment of routines written in IDL for
reducing and analyzing multiwavelength imaging and spectroscopy of
galaxies.  

Some of the code highlights included with impro:

* iSEDfit - fully developed, documented, and publicly released software
  which uses Bayesian inference to extract the physical properties of
  galaxies from their observed broadband photometric spectral energy
  distributions.  Please see
  http://www.sos.siena.edu/~jmoustakas/isedfit for details.

* Routines to read and interact with the outputs from a broad range of
  stellar population synthesis models, including BC03, Pegase,
  Starburst99, and FSPS.

* Routines for inferring the physical conditions (metallicity,
  temperature, dust content, etc.) and star-formation rates of
  star-forming galaxies and HII regions from the strength of their
  nebular emission lines.

* Routines to construct and model the luminosity and stellar mass
  functions of galaxies using a variety of parameterizations (single,
  double Schechter, double power-law, etc.).

* Routines to fit and model the optical stellar continuum and
  emission-line spectra of galaxies, and to measure a variety of
  spectral indices (e.g., Lick indices, D4000, etc.).

As usual, the code is provided as-is, with no guarantee of accuracy or
success.  However, if you are interested in contributing actively to
impro send me a note!

### Installation

Download the bleeding-edge version of the repository from Github using
the following command:

    git clone http://github.com/moustakas/impro.git

Alternatively, you can download a tarball or zip file snapshot of the
code base by navigating directly to github.com:

    https://github.com/moustakas/impro

Next, add an ${IMPRO_DIR} environment variable to your appropriate
startup file. For example, if you use tcsh, add the following line to
your .tcshrc file: 

    setenv IMPRO_DIR /path/to/impro

Finally, add the 'pro' subdirectory of ${IMPRO_DIR} to your IDL path.
Using the example above, you would add:

    setenv IDL_PATH $IDL_PATH{:}+${IMPRO_DIR}/pro

to either your .tcshrc or .idlenv startup file.  

