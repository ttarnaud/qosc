SONIC and SECONIC solver

1) SONICtabulate will calculate the look-up tables based on the BLS-model
2) SONICrun will simulate the modified HH equations: 
	a. SONICrun: single compartment SONIC model
	b. SONICrun_nanoMC: nanoscale two compartment SONIC model
	c. SONICrun_nanoMC_Qosc: nanoscale two compartment SECONIC model

Input/output arguments 

1) TTime is the total simulation time
All inputs are strings (e.g. Tsim = '0.3' for 300 ms simulation duration)
2) Tsim: simulation duration [s]
3) MODE: simulation mode --  1: titration to excitation threshold (minimal ultrasonic intensity that induces an action potential)
						 --  2: Simulation of a fixed ultrasonic intensity
4) USpstart: Ultrasonic pulse start [s]
5) USpd: Ultrasonic pulse duration [s]
6) USfreq: Ultrasonic frequency [Hz]
7) USdc: Ultrasonic duty cycle [-]
8) USprf: Ultrasonic pulse repetition frequency [Hz]
9) USisppa: [ignored if MODE=1] Spatial peak pulsed average ultrasonic intensity [W/m^2]	
10) ESpstart,ESpd,ESdc,ESprf,ESisppa [s,s,s,Hz,A/m^2]: electrical waveform parameters
11) PLOT: plot mode (0: no plots, 1: plot results, 2: save plot data as .mat file)
12) model: Neuronal model: 1=Regular spiking cortical cell [1,7],2=fast spiking cortical cell [1,7],3=low threshold spiking cortical cell [1,7],
4=Thalamocortical cell [2,7], 5=Thalamus reticular neuron [2,7],
6=Regular spiking ferret visual cortex cell [1,7], 7=Fast spiking ferret visual cortex cell [1,7], 8=Low threshold spiking cat association cortex cell [1,7],
9=Subthalamic nucleus [3], 10=Thalamus Rubin Terman model [4], 11=Globus pallidus interna [4], 12=Globus pallidus externa [4],
13=Medium spiny neuron (striatum) [5], 14=Hodgkin–Huxley model [6].
13) USibegin, USiend: (ignored if MODE~=1) Lower and upper boundary for the excitation threshold. (Choose 0 if not known) 
14) SearchMode: (ignored if MODE~=1, boolean) - 0: Lower and/or upper boundary not known (USibegin == 0 || USiend == 0) = 1.
											  - 1: Both lower and upper boundary USibegin and USiend are known. 
15) aBLS: Sonophore radius [m]
16) fBLS: Sonophore coverage factor [-] (i.e. ratio of bilayer sonophore area and total area (including the protein islands))
17) proteinMode [Default 0]: coverage of proteins: (mode 0: full coverage; mode 1: partial coverage, leak current has full coverage; 
mode 2: all gates partial coverage including leak current) 
18) threshMode [Default 0]: (Boolean, 0: neuron excitation during Tsim for threshold determination; 1: only neuron excitation during stimulus duration).
19) gateMultip [Default 1]: multiplier, applied to the conductivity gains that are reduced if proteinMode ~= 0. 
20) CorrPec [Boolean, Default 0]: If 1, correct the electrostatic pressure (Pec) for fBLS < 1, to account for lateral charge redistribution

[1] Pospischil, M., Toledo-Rodriguez, M., Monier,... & Destexhe, A. (2008). 
Minimal Hodgkin–Huxley type models for different classes of cortical and thalamic neurons. Biological cybernetics, 99(4-5), 427-441.
[2] Destexhe, A., Contreras, D., & Steriade, M. (1998). Mechanisms underlying the synchronizing action of corticothalamic feedback through inhibition of thalamic relay cells. 
Journal of neurophysiology, 79(2), 999-1016.
[3] Otsuka, T., Abe, T., Tsukagawa, T., & Song, W. J. (2004). Conductance-based model of the voltage-dependent generation of a plateau potential in subthalamic neurons. 
Journal of neurophysiology, 92(1), 255-264.
[4] Rubin, J. E., & Terman, D. (2004). High frequency stimulation of the subthalamic nucleus eliminates pathological thalamic rhythmicity in a computational model.
Journal of computational neuroscience, 16(3), 211-235.
[5] McCarthy, M. M., Moore-Kochlacs, C., Gu, X., Boyden, E. S., Han, X., & Kopell, N. (2011). Striatal origin of the pathologic beta oscillations in Parkinson's disease. 
Proceedings of the National Academy of Sciences, 108(28), 11620-11625.
[6] Hodgkin, A. L., & Huxley, A. F. (1952). A quantitative description of membrane current and its application to conduction and excitation in nerve. 
The Journal of physiology, 117(4), 500-544.
[7] Plaksin, M., Kimmel, E., & Shoham, S. (2016). Cell-type-selective effects of intramembrane cavitation as a unifying theoretical framework for ultrasonic neuromodulation. eneuro, 3(3).

Example testrun: 
>> SONICrun('0.3','2','0.1','0.1','500e3','1','0','100','0','0','1','0','0','1','9','0','0','1','32e-9','1')

- Dependencies: 
interp1qr (Jose M. Mier), Odextend (Jacek Kierzenka), Nakeinterp1 (Bruno Luong)

- Email for bugs/questions: thomas.tarnaud@ugent.be
