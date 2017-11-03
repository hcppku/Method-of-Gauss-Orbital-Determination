from __future__ import division
from math import*
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

#Returns true angle value based on sine and cosine values
def findQuadrant(sine, cosine):
    if cosine > 0 and sine > 0: #1
        return asin(sine)

    if cosine < 0 and sine > 0: #2
        return acos(cosine)

    if cosine < 0 and sine < 0: #3
        return pi - asin(sine)

    if cosine > 0 and sine < 0: #4
        return 2*pi + asin(sine)

#Makes angle between 0 and 2pi
def adjustAngle(angle):
    while angle>2*pi:
        angle=angle-2*pi
    while angle<0:
        angle=angle+2*pi
    return angle

#Finds magnitude of a vector
def magnitude(vector):
    squaresums=0
    for component in range(0,len(vector)):
        squaresums=squaresums+(vector[component])**2
    return sqrt(squaresums)

#Finds a
def aCalc(r,rdot):
    return 1/((2/magnitude(r))-rdot.dot(rdot))

#Finds e
def eCalc(r,rdot,a):
    return sqrt(1-(magnitude(np.cross(r,rdot))**2/a))

#Finds I*Depends
def ICalc(r,rdot):
    return adjustAngle(acos(np.cross(r,rdot)[2]/magnitude(np.cross(r,rdot))))

#Finds O
def OCalc(r,rdot,I):
    h=np.cross(r,rdot)
    sinO=h[0]/(magnitude(h)*sin(I))
    cosO=-h[1]/(magnitude(h)*sin(I))
        
    #if h[2]<0:
        #sinO=-h[0]/(magnitude(h)*sin(I))
        #cosO=h[1]/(magnitude(h)*sin(I))
    return adjustAngle(atan2(sinO,cosO))

#Finds w
def wCalc(r,rdot,a,e,I,O):
    sinwplusf=r[2]/(magnitude(r)*sin(I))
    coswplusf=(1/cos(O))*(r[0]/magnitude(r)+cos(I)*sinwplusf*sin(O))
    wplusf=atan2(sinwplusf,coswplusf)
    cosf=(1/e)*(a*(1-e**2)/magnitude(r)-1)
    sinf=(a*(1-e**2)/(magnitude(np.cross(r,rdot))))*(r.dot(rdot)/(e*magnitude(r)))
    f=findQuadrant(sinf,cosf)
    return adjustAngle(wplusf-f)

#Finds f
def fCalc(r,a,e,rdot):
    cosf=(1/e)*(a*(1-e**2)/magnitude(r)-1)
    sinf=(a*(1-e**2)/(magnitude(np.cross(r,rdot))))*(r.dot(rdot)/(e*magnitude(r)))
    f=findQuadrant(sinf,cosf)
    return f

#Finds M
def MCalc(r,rdot,a,e):
    f=fCalc(r,a,e,rdot)
    cosE=((1-(magnitude(r)/a))/e)
    sinE=sqrt(1-e**2)*sin(f)/(1+e*cos(f))
    E=findQuadrant(sinE,cosE)
    return E-e*sin(E)

#Returns orbital elements in this order: a,e,I,O,w,M
def orbElems(r,rdot):
    rot=np.array([[1,0,0],[0,cos(radians(23.437)),-sin(radians(23.437))],[0,sin(radians(23.437)),cos(radians(23.437))]])
    rot=np.linalg.inv(rot)
    r=rot.dot(r)
    rdot=(rot.dot(rdot))/(0.01720209894)
    a=aCalc(r,rdot)
    e=eCalc(r,rdot,a)
    I=ICalc(r,rdot)
    O=OCalc(r,rdot,I)
    w=wCalc(r,rdot,a,e,I,O)
    M=MCalc(r,rdot,a,e)
    I=degrees(I)
    O=degrees(O)
    w=degrees(w)
    M=degrees(M)
    return a,e,O,I,w,M

loop=True
while loop==True:

    print "Enter r and rdot for equatorial system: conversion to ecliptic is automatic"
    rx=float(raw_input("Enter x component of position vector: "))
    ry=float(raw_input("Enter y component of position vector: "))
    rz=float(raw_input("Enter z component of position vector: "))
    rdotx=float(raw_input("Enter x component of velocity vector: "))
    rdoty=float(raw_input("Enter y component of velocity vector: "))
    rdotz=float(raw_input("Enter z component of velocity vector: "))
    r=np.array([rx,ry,rz])
    rdot=np.array([rdotx,rdoty,rdotz])
    print "orbital elements in order of a,e,O,I,w,M"
    print orbElems(r,rdot)
    choice=int(raw_input("Enter 1 to do another calculation or 2 to end: "))
    if choice==2:
        loop=False

    
        
#Test:
##r=np.array([1,0,2])
##rdot=np.array([0,-0.01,0.01])
##print(orbElems(r,rdot))
##r2=np.array([2.54076127,0.84513781,0.3806295])
##r2dot=np.array([-0.003577,0.00943291,0.00429106])
##print(orbElems(r2,r2dot))
##r3=np.array([1.57,0.0677,0.253])
##r3dot=np.array([-0.00169,0.00985,0.00145])
##print(orbElems(r3,r3dot))
