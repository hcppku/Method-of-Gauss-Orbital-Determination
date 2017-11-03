from __future__ import division
import numpy as np
from math import*

def findQuadrant(sine, cosine):
    if cosine > 0 and sine > 0: #1
        return asin(sine)

    if cosine < 0 and sine > 0: #2
        return acos(cosine)

    if cosine < 0 and sine < 0: #3
        return pi - asin(sine)

    if cosine > 0 and sine < 0: #4
        return 2*pi + asin(sine)

def degToRad(degs):
    return degs*pi/180

def radToDeg(rads):
    return rads*180/pi

#Updates M based on time
def Mcalc(a,T,Mo,tOld):
    n=0.01720209894*sqrt(1.000000/a**3)
    return n*(T-tOld)+Mo

#Newton's Method for approximating E
def newton(M,e):
    Eguess = M
    Mguess = Eguess - e*sin(Eguess)
    while abs(Mguess - M) > 1e-004:
        Mguess = Eguess - e*sin(Eguess)
        Eguess = Eguess - (M-(Eguess-e*sin(Eguess)))/(e*cos(Eguess)-1)
    
    return Eguess

#The function f that is used within Newton's Method
def functionf(M,E,e):
    return M-(E-e*sin(E))

#Position Vector
def position(a,E,e):
    x=a*cos(E)-a*e
    y=a*sqrt(1-e**2)*sin(E)
    z=0
    pos=np.array([x,y,z])
    return pos

#Ecliptic position from cartesian position
def eclipticPos(pos,O,i,w):
    array1=np.array([[cos(O),-sin(O),0],[sin(O),cos(O),0],[0,0,1]])
    
    array2=np.array([[1,0,0],[0,cos(i),-sin(i)],[0,sin(i),cos(i)]])
    
    array3=np.array([[cos(w),-sin(w),0],[sin(w),cos(w),0],[0,0,1]])
    
    arr=array1.dot(array2).dot(array3)
    
    return arr.dot(pos)

#Equatorial position from ecliptic position
def equatorial(ecpos):
    array=np.array([[1,0,0],[0,cos(degToRad(23.4)),-sin(degToRad(23.4))],[0,sin(degToRad(23.4)),cos(degToRad(23.4))]])
    return array.dot(ecpos)

#Finds RA and Dec. Enter values in AU and degrees. RA is in hours and dec is in degrees
def generateRaDec(a,e,i,O,w,Mo,tOld,T,earthSunVec):
    i=degToRad(i)
    O=degToRad(O)
    w=degToRad(w)
    Mo=degToRad(Mo)
    r=equatorial(eclipticPos(position(a,newton(Mcalc(a,T,Mo,tOld),e),e),O,i,w))
    R=earthSunVec
    p=r+R
    pmag=sqrt(p[0]**2+p[1]**2+p[2]**2)
    punit=p/pmag
    dec=asin(punit[2])
    cosRA=punit[0]/cos(dec)
    sinRA=punit[1]/cos(dec)
    RA=findQuadrant(sinRA,cosRA)
    RA=radToDeg(RA)/15
    dec=radToDeg(dec)
    return RA,dec

#Test:
#print(generateRaDec(1.844166011011208,0.09619878680571502,23.66310678852568,132.1025646044957,129.0712937388337,315.8154532795419,2457600.5,2457753.5,[0.1622529802114096,-0.8898495228908435,-0.3857824315413861]))
loop=True
while loop==True:
    a=float(raw_input("Enter a: "))
    e=float(raw_input("Enter e: "))
    i=float(raw_input("Enter i: "))
    O=float(raw_input("Enter O: "))
    w=float(raw_input("Enter w: "))
    Mo=float(raw_input("Enter Mo: "))
    to=float(raw_input("Enter JD for original time: "))
    t=float(raw_input("Enter JD for desired time: "))
    x=float(raw_input("Enter x for earth Sun vector: "))
    y=float(raw_input("Enter y for earth Sun vector: "))
    z=float(raw_input("Enter z for earth Sun vector: "))
    print (generateRaDec(a,e,i,O,w,Mo,to,t,[x,y,z]))

    choice=int(raw_input("Enter 1 for another calculation or 2 to end: "))
    if choice==2:
        loop=False
        
    
    

