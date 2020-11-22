@ Oskar Hoegberg 13/02/2020

## Usage
	
All files should be set to generate demo plots and 
contain instructions for creating custom initial
conditions.

pa_plots.py and pt_plots.py will generate the more
complicared plots and can take very long to run for 
advanced setups. pa stands for pressure-area and
investigates the relation between the two. pt means
pressure-temperature and is dedicated to this 
relationship.

runtime.py is used to simply run the 
simulation once with different configurations.

Generally, either temperatrue or container area is
varied to measure the effect this has on pressure.

## Remarks
The parameters ``separations`` in the run method
slows down the code a fair bit so do not have it 
set to true unless this data is specifically wanted.
