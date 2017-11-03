from __future__ import division
from math import sqrt
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt


loop=True
while loop==True:
    fileName=raw_input("Enter file name: ")

    image=fits.getdata(fileName)
    print fits.getheader(fileName)
    def display(image):
        plt.imshow(image)
        plt.gray()
        plt.imshow(image,vmin=image.mean(),vmax=2*image.mean())
        plt.show()
    display(image)

    choice=int(raw_input("Enter 1 to loop or 2 to end: "))
    if choice==2:
        loop=False

