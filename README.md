# Laser-interaction-in-fields-rate-equations-forces

Rate equations simulation of particles in fields

All details are found in the USER GUIDE. Here is simply the most relevant part. Modifications (new versions) are in Modif_code.txt

The program solves the rate equations to study laser excitation, forces (scattering + dipolar + magnetic + electric + coulombian
interactions). 


It has been developed under Code::Blocks and Windows. The inputs are 2 external files describing the levels (with information about their energy + linear or quadratic Stark, Zeeman effect) and the transitions lines (dipole transitions, photodetachement or photoionization cross sections).

A file named Liste Param.h contains parameters (nad lots of useful comments on them) to run the simulation such as:
sample size, temperature, magnetic fields and for the laser beams (waist size and position, polarisation, power, linewidth, wavelength, ...). 

A Kinetic Monte Carlo algorithm  evolves in motion (forces under extranal fields + optical dipole forces) and event (absorption or emission) 

The output is writen in a file containing relevant information such as population in given levels and statistics about velocities (temperature), potential energy ... Output is also performed through 3D snapshots. 

