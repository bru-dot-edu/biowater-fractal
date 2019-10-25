Ref:
Almquist et. al. 2018 - Bioavailable Water in Coarse Soils: A Fractal
Approach. 

DOI:
https://doi.org/10.1016/j.geoderma.2018.02.036

This is a readme file to help get you going with the model so that you can
generate a thick film contribution pressure-saturation curve for your own
soil sample. 

To start, you're going to need your own particle size distribution data. This
model will be most suitable if you don't treat your sample with any hydrogen
peroxide prior to measuring the distribution. The file I use (which is in the
repository) is called 'particle_size_data.txt'.

If you have some moisture retention data, this code will also generate a
van Genuchten curve. But that is not necessary for the thick film model. My
version of that file is in the repository as 'quincy_moisture_data.txt'.

Once you have those files in your directory, run the .py files in this order:

1. particle_size_fit_and_bin.py
   a. this will output binned_PSD_mass_fraction.txt which you should put in your
   directory 
2. (optional) van_gen_fit.py
   a. You'll need to manually copy the van genuchten parameters into ffm.p. in
   the "constants" section depending on how popular the model is I may update
   this to make it more user friendly.
3. ffm.py
   a. this is the main model as described in the manuscript. if you've opted to
   ignore the van genuchten curve, I recommend just leaving the van genuchten
   parameters so you don't receive any errors. Of course, the van genuchten
   curve won't be applicable to your soil, and neither will the capillary
   contribution curve. But the thick film model will still be valid.  

Have fun!


