# Method-of-Gauss-Orbital-Determination
The following Python files were written to perform a variety of calculations to determine the orbital elements of an asteroid (or any sun-orbiting body) using .fit image files from a telescope

Note that these Python files use Python 2.7 as well as the numpy, matplotlib, vpython and astropy libraries

Below, I will list out the files of this repository and give a brief description

1. Gauss Orbit Determination - takes in a txt file (see 2002 KL6Final.txt as a template) that takes in the date of an observation, the time of the observation (UT), the RA and Declination of the object at that time, and the Sun-Earth vector in au. The code uses the method of gauss to calculate the orbital elements and uses differential correction to improve results. In addition, a vpython animation of the asteroid is created.

2. AstronomyTimeCalculator - calculates time conversions that are useful in Astronomy such as Sidereal time, Julian time, and UT

3. fitsInfoGetter - shows a .fit file as well as the .fit heading

4. LSPR - takes in a txt file with a list of information on background stars in a .fit file. Each line includes the x and y of the centroid of the star, and then the ra and declination of that star. LSPR then asks for the x and y of the centroid of the object of interest (asteroid). LSPR uses least-squares plate reduction to calculate the ra and declination of the asteroid which is the output.

5. centroidprogram - takes in a .fit file, approximate x and y positions of the centroid of an object, and then the dimensions of the centroid box. The program that calculates the exact centroid of the object.

6. magnitude circular Jas - This calculates the apparent magnitude of the asteroid using a least-squares method of background star magnitudes. Takes in a .fit file as well as information on stars of your choosing.

7. Ephemeris Generator - Knowing orbital elements, the time corresponding to the mean anamoly, and the earth Sun vector, you can find the ephemeris for an object.

8. randrDot to Orbitals - this piece of code calculates orbital elements given the velocity and position vectors of an object from the sun at a single moment in time.

9. JackKnife Uncertainty - Jacknife method is a simple method to calculate uncertainties in a data set
