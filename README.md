# Uncertainty Quantification software package in gas pipes

This is a software package written in C++ and Python for solving
Bayesian inverse problems that occur in the Gas industry, e.g., to
identify friction coefficients in gas pipes.

# Quick Installation on Mac OS

requirement: make, g++

download the package either through git,
```bash
	git clone URL
```
or download it directly from URL and decompress it. Go into the
package library and execute the following commands:
```bash
	make lib/CWrapper.dylib
	python setup.py install
```

# Testing
If everything goes well you should be able to test the package:
```bash
	python examples/test.py
```	
The output in the console looks like
```html
========================================== info ==
Omega _____________________________(0,1)
Length of the domain ______________1
Dt ________________________________0.005
Dx ________________________________0.005
Number of cells ___________________200
Current time ______________________0
Final time ________________________5
Epsilon for the boundary __________0.05
N_epsilon _________________________10
Size of friction coef. vec. _______7
==================================================
>> Computation Done
```
and also the following figure which shows continuous dependence of the
pressure drop at both ends of the pipe with respect to the friction
coefficient.
![Pressure drop](results/pressure_drop_test.png)
