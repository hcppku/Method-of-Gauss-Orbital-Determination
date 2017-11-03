#Method of Gauss Orbit Determination
#Jason Kim: SSP 2016: 7/12/16

from __future__ import division
import visual
from visual import *
from math import*
import numpy as np

tol = 10**-12

#*****************************************************File Reader****************************************
def JDNumber(date,time):
    year=float(date[0:4])
    month=float(date[5:7])
    day=float(date[8:])
    hour=float(time[0:2])
    minute=float(time[3:5])
    second=float(time[6:])
    decimalTime=hour+minute/60+second/3600
    Jo=367*year-int(7*(year+int((month+9)/12))/4)+int(275*month/9)+day+1721013.5
    JD=Jo+decimalTime/24
    return JD

def decimalCoord(raDec):
    firstColon=raDec.index(":")
    hourDeg=float(raDec[0:firstColon])
    minute=float(raDec[firstColon+1:firstColon+3])
    second=float(raDec[firstColon+4:])
    return hourDeg+minute/60+second/3600

#*****************************************************METHOD OF GAUSS************************************
#From given right ascension and declination values, a UNIT vector (rhoHat or pHat) from Earth to an object orbiting the Sun is returned. RA is given in hours and Declination is given in degrees
def calcpHat(ra,dec):
    ra=radians(ra*15)
    dec=radians(dec)
    pHat=np.array([cos(ra)*cos(dec),sin(ra)*cos(dec),sin(dec)])

    return pHat

#Returns Dknot (Do) value given rhoHat(pHat) vectors
def calcDo(pHat1,pHat2,pHat3):
    return pHat1.dot(np.cross(pHat2,pHat3))

#Reurns array with all Dij values where i and j go from 1 to 3. i determines the row. j determines the column
def calcD(pHat1,pHat2,pHat3,R1,R2,R3):
    return np.array([[np.cross(R1,pHat2).dot(pHat3),np.cross(R2,pHat2).dot(pHat3),np.cross(R3,pHat2).dot(pHat3)],[np.cross(pHat1,R1).dot(pHat3),np.cross(pHat1,R2).dot(pHat3),np.cross(pHat1,R3).dot(pHat3)],[pHat1.dot(np.cross(pHat2,R1)),pHat1.dot(np.cross(pHat2,R2)),pHat1.dot(np.cross(pHat2,R3))]])

#Returns tao (Gaussian time) values in array in order of tao, tao1, tao3 given time values for t1,t2,t3 for three observations
def calcTao(t1,t2,t3):
    k=0.01720209894
    return np.array([k*(t3-t1),k*(t2-t1),k*(t3-t2)])

#Returns A1,A3,B1,B3 in an array in that order given tao, tao1, tao3
def calcAB13(tao,tao1,tao3):
    A1=tao3/tao
    A3=tao1/tao
    B1=A1*(tao**2-tao3**2)/6
    B3=A3*(tao**2-tao1**2)/6
    return np.array([A1,A3,B1,B3])

#Returns array with A,B values given A1,A3,B1,B3 in an array and D values and Do
def calcAB(AB13,D,Do):
    A=(AB13[0]*D[1,0]-D[1,1]+AB13[1]*D[1,2])/(-Do)
    B=(AB13[2]*D[1,0]+AB13[3]*D[1,2])/(-Do)
    return np.array([A,B])

#Returns E,F values in an array given pHat2,R2
def calcEF(pHat2,R2):
    E=-2*(pHat2.dot(R2))
    F=np.linalg.norm(R2)**2
    return np.array([E,F])

#Returns a,b,c values in an array given A,B and E,F
def calcr2Mag(AB,EF):
    a=-(AB[0]**2+AB[0]*EF[0]+EF[1])
    b=-(2*AB[0]*AB[1]+AB[1]*EF[0])
    c=-(AB[1]**2)
    return np.polynomial.polynomial.polyroots([c,0,0,b,0,0,a,0,1])

#Returns initial c1,c2,c3
def initialC(AB13,r2Mag):
    c1=AB13[0]+AB13[2]/(r2Mag**3)
    c2=-1
    c3=AB13[1]+AB13[3]/(r2Mag**3)
    return c1,c2,c3
    
#Returns r vector values for r1,r2,r3 and new t values based on light travel time correction
def calcr(Do,D,pHat1,pHat2,pHat3,R1,R2,R3,c1,c2,c3,t1,t2,t3):
    
    p1Mag=(c1*D[0,0]+c2*D[0,1]+c3*D[0,2])/(c1*Do)
    p2Mag=(c1*D[1,0]+c2*D[1,1]+c3*D[1,2])/(c2*Do)
    p3Mag=(c1*D[2,0]+c2*D[2,1]+c3*D[2,2])/(c3*Do)
    
    #Light travel time:
    ta=t1-p1Mag/173.144633 
    tb=t2-p2Mag/173.144633 
    tc=t3-p3Mag/173.144633 
    
    r1=p1Mag*pHat1-R1
    r2=p2Mag*pHat2-R2
    r3=p3Mag*pHat3-R3
    return r1,r2,r3,ta,tb,tc

#Returns r2Dot value
def calcr2Dot(r2Mag,r,tao1,tao3):
    u2=1/(r2Mag**3)
    f1=1-0.5*u2*(tao1**2)
    f3=1-0.5*u2*(tao3**2)
    g1=-tao1+u2*(tao1**3)/6
    g3=tao3-u2*(tao3**3)/6
    d1=-f3/(f1*g3-f3*g1)
    d3=f1/(f1*g3-f3*g1)
    r2Dot=d1*r[0]+d3*r[2]
    return r2Dot

#Remember to put in new taos for this function
#Loops to give more accurate values of r2Dot using Taylor Series
def loop(rOld,r2DotOld,tao1,tao3,Do,D,pHat1,pHat2,pHat3,R1,R2,R3,t1,t2,t3):
    r2Old=rOld[1]
    myvar=0
    
    while True:
        myvar+=1
        r2Mag=np.linalg.norm(r2Old)
        u=1/(r2Mag**3)
        z=(r2Old.dot(r2DotOld))/(r2Mag**2)
        q=(r2DotOld.dot(r2DotOld))/(r2Mag**2)-u
        f1=1-0.5*u*tao1**2-0.5*u*z*tao1**3+(1/24)*(3*u*q-15*u*z**2+u**2)*tao1**4
        f3=1-0.5*u*tao3**2+0.5*u*z*tao3**3+(1/24)*(3*u*q-15*u*z**2+u**2)*tao3**4
        g1=-(tao1-(1/6)*u*tao1**3+0.25*u*z*tao1**4)
        g3=tao3-(1/6)*u*tao3**3+0.25*u*z*tao3**4
        c1=g3/(f1*g3-g1*f3)
        c2=-1
        c3=-g1/(f1*g3-g1*f3)
        p1Mag=(c1*D[0,0]+c2*D[0,1]+c3*D[0,2])/(c1*Do)
        p2Mag=(c1*D[1,0]+c2*D[1,1]+c3*D[1,2])/(c2*Do)
        p3Mag=(c1*D[2,0]+c2*D[2,1]+c3*D[2,2])/(c3*Do)

        rNews=calcr(Do,D,pHat1,pHat2,pHat3,R1,R2,R3,c1,c2,c3,t1,t2,t3)
        
        newTaos=calcTao(rNews[3],rNews[4],rNews[5])
        tao=newTaos[0]
        tao1=newTaos[1]
        tao3=newTaos[2]
        
        r2=rNews[1]

        d1=-f3/(f1*g3-f3*g1)
        d3=f1/(f1*g3-f3*g1)
        r2Dot=d1*rNews[0]+d3*rNews[2]
        if abs(r2[0]-r2Old[0])<tol and abs(r2[1]-r2Old[1])<tol and abs(r2[2]-r2Old[2])<tol and abs(r2Dot[0]-r2DotOld[0])<tol and abs(r2Dot[1]-r2DotOld[1])<tol and abs(r2Dot[2]-r2DotOld[2])<tol:
            break
        r2Old=r2
        r2DotOld=r2Dot
    #print myvar
    return r2,r2Dot, p2Mag


def methodOfGauss(t1,ra1,dec1,R1,t2,ra2,dec2,R2,t3,ra3,dec3,R3):
    pHat1=calcpHat(ra1,dec1)
    pHat2=calcpHat(ra2,dec2)
    pHat3=calcpHat(ra3,dec3)
    
    Do=calcDo(pHat1,pHat2,pHat3)
    D=calcD(pHat1,pHat2,pHat3,R1,R2,R3)
    
    taos=calcTao(t1,t2,t3)
    
    tao=taos[0]
    tao1=taos[1]
    tao3=taos[2]
    AB13=calcAB13(tao,tao1,tao3)
    AB=calcAB(AB13,D,Do)
    EF=calcEF(pHat2,R2)
    
    roots=calcr2Mag(AB,EF)
    r2Mag=abs(roots[int(raw_input(str(roots)+"\nPlease enter which of these r2 magnitudes you would like to proceed with. Enter the number of the r2 magnitude you want. Enter 1 for the first, 2 for the second, and so on."))-1])
    #r2Mag=abs(roots[7])
    firstCs=initialC(AB13,r2Mag)
    c1=firstCs[0]
    c2=firstCs[1]
    c3=firstCs[2]
    rVectorsFirst=calcr(Do,D,pHat1,pHat2,pHat3,R1,R2,R3,c1,c2,c3,t1,t2,t3)
    
    ta=rVectorsFirst[3]
    tb=rVectorsFirst[4]
    tc=rVectorsFirst[5]
    taos=calcTao(ta,tb,tc)
    tao=taos[0]
    tao1=taos[1]
    tao3=taos[2]
    r2DotFirst=calcr2Dot(r2Mag,rVectorsFirst,tao1,tao3)
    
    
    iteration=loop(rVectorsFirst,r2DotFirst,tao1,tao3,Do,D,pHat1,pHat2,pHat3,R1,R2,R3,t1,t2,t3)
    r2=iteration[0]
    r2Dot=iteration[1]
    p2Mag=iteration[2]
    return r2,r2Dot,p2Mag


#****************************************BABY OD (Determining Orbital Elements)********************************
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
    
    a=aCalc(r,rdot)
    e=eCalc(r,rdot,a)
    rot=np.array([[1,0,0],[0,cos(radians(23.437)),-sin(radians(23.437))],[0,sin(radians(23.437)),cos(radians(23.437))]])
    rot=np.linalg.inv(rot)
    r=rot.dot(r)
    rdot=(rot.dot(rdot))/(0.01720209894)
    I=ICalc(r,rdot)
    O=OCalc(r,rdot,I)
    w=wCalc(r,rdot,a,e,I,O)
    M=MCalc(r,rdot,a,e)
    I=degrees(I)
    O=degrees(O)
    w=degrees(w)
    M=degrees(M)
    return a,e,O,I,w,M


#*************************************VPython Simulation***********************************
#Updates M based on time
def MVcalc(a,k,mu,T):
    n=k*sqrt(mu/a**3)
    return n*T

def Mcalc2(a,T,Mo,tOld):
    n=0.01720209894*sqrt(1.000000/a**3)
    return n*(T-tOld)+Mo

#Newton's Method for approximating E
def newton(M,e):
    Eguess = M
    Mguess = Eguess - e*sin(Eguess)
    while abs(Mguess - M) > 1e-004:
        Eguess = Eguess - functionf(M,Eguess,e)/(e*cos(Eguess)-1)
        Mguess = Eguess - e*sin(Eguess)
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
def equatorial(ecpos,ec):
    array=np.array([[1,0,0],[0,cos(ec),-sin(ec)],[0,sin(ec),cos(ec)]])
    return array.dot(ecpos)

def vPythSim(a,e,O,I,w,M):
    visual.scene.autoscale=False
    #Asteroid
    a=a
    e=e
    O=radians(O)
    I=radians(I)
    w=radians(w)
    M=radians(M)
    positionEq=equatorial(eclipticPos(position(a,newton(M,e),e),O,I,w),radians(23.437))
    rVector=visual.vector(0,0,0)
    rVector.x=positionEq[0]
    rVector.y=positionEq[1]
    rVector.z=positionEq[2]

    #Earth
    earthA=1.000703137122473
    earthe=1.599048849057274E-02
    earthO=radians(1.652211903265974E+02)
    earthI=radians(2.744957205807269E-03)
    earthW=radians(2.973932488195502E+02)
    earthM=radians(2.061751020668477E+02)
    
    earthPos=equatorial(eclipticPos(position(earthA,newton(earthM,earthe),earthe),earthO,earthI,earthW),radians(23.437))
    erVector=visual.vector(0,0,0)
    erVector.x=earthPos[0]
    erVector.y=earthPos[1]
    erVector.z=earthPos[2]

##    def MVcalc(a,k,mu,T):
##        n=k*sqrt(mu/a**3)
##        return n*T
##
##    def Mcalc2(a,T,Mo,tOld):
##        n=0.01720209894*sqrt(1.000000/a**3)
##        return n*(T-tOld)+Mo

    t=0
    #Where 0 corresponds to the JD of the central observation (where simulation starts)
    dt=0.5
    asteroid = visual.sphere(pos=rVector*150, radius=(5), color=color.white)
    asteroid.trail = visual.curve(color=color.white)
    Earth = visual.sphere(pos=erVector*150, radius=(15), material=materials.earth)
    Earth.trail = visual.curve(color=color.white)
    sun = visual.sphere(pos=(0,0,0), radius=(50), color=color.yellow)
    while(1==1):
        visual.rate(100)
        t=t+1
        M=MVcalc(a,0.01720209894,1,t)
        E=newton(M,e)
        positionEq=equatorial(eclipticPos(position(a,E,e),O,I,w),radians(23.437))
        rVector.x=positionEq[0]
        rVector.y=positionEq[1]
        rVector.z=positionEq[2]
        asteroid.pos = rVector*150
        asteroid.trail.append(pos=asteroid.pos)
        
        earthM=MVcalc(earthA,0.01720209894,1,t)
        earthE=newton(earthM,earthe)
        earthPos=equatorial(eclipticPos(position(earthA,earthE,earthe),earthO,earthI,earthW),radians(23.437))
        erVector.x=earthPos[0]
        erVector.y=earthPos[1]
        erVector.z=earthPos[2]
        Earth.pos = erVector*150
        Earth.trail.append(pos=Earth.pos)
    return True

#*************************************Ephemeris Generator***********************************
def generateRaDec(a,e,i,O,w,Mo,tOld,T,earthSunVec):
    i=radians(i)
    O=radians(O)
    w=radians(w)
    Mo=radians(Mo)
    
    r=equatorial(eclipticPos(position(a,newton(Mcalc2(a,T,Mo,tOld),e),e),O,i,w),radians(23.437))
    
    R=earthSunVec
    p=r+R
    pmag=sqrt(p[0]**2+p[1]**2+p[2]**2)
    punit=p/pmag
    dec=asin(punit[2])
    cosRA=punit[0]/cos(dec)
    sinRA=punit[1]/cos(dec)
    RA=findQuadrant(sinRA,cosRA)
    RA=degrees(RA)
    dec=degrees(dec)
    return RA,dec


#*************************************Orbit Improvement************************************
def OminusC1(ra,dec,radec):
    summ=0
    for j in range(0,len(ra)):
        summ=summ+(ra[j]*15-(radec[j])[0])**2+(dec[j]-(radec[j])[1])**2
    rms=sqrt(summ/(len(ra)*2))
    return rms

#respectTo is int where 0 is x,1 is y, 2 is z, 3 is xdot, 4 is ydot, 5 is zdot. raordec is 0 for ra, 1 for dec
def partial(j,raObserved,decObserved,radecInitial,r,rdot,respectTo,raordec,tOld,T,earthSunVec):
    rrdot=np.array([r[0],r[1],r[2],rdot[0],rdot[1],rdot[2]])#should concatenate these two
    delta=rrdot[respectTo]*10**-4
    rrdotMinus=np.array(rrdot)
    rrdotPlus=np.array(rrdot)
    rrdotMinus[respectTo]=rrdotMinus[respectTo]-delta
    rrdotPlus[respectTo]=rrdotPlus[respectTo]+delta
    orbMinus=orbElems(np.array(rrdotMinus[0:3]),np.array(rrdotMinus[3:6]))
    aminus=orbMinus[0]
    eminus=orbMinus[1]
    Ominus=orbMinus[2]
    Iminus=orbMinus[3]
    wminus=orbMinus[4]
    Mminus=orbMinus[5]
    
    orbPlus=orbElems(np.array(rrdotPlus[0:3]),np.array(rrdotPlus[3:6]))
    aplus=orbPlus[0]
    eplus=orbPlus[1]
    Oplus=orbPlus[2]
    Iplus=orbPlus[3]
    wplus=orbPlus[4]
    Mplus=orbPlus[5]
    alMinus=generateRaDec(aminus,eminus,Iminus,Ominus,wminus,Mminus,tOld,T,earthSunVec)[raordec]
    alPlus=generateRaDec(aplus,eplus,Iplus,Oplus,wplus,Mplus,tOld,T,earthSunVec)[raordec]
    return (alPlus-alMinus)/(2*delta),delta

##def deltAi(j,raObserved,decObserved,radecInitial,r,rdot,respectTo,raordec,tOld,T,solarVec):#respectTo value not used and raordec not used
##    dai_dx=partial(j,raObserved,decObserved,radecInitial,r,rdot,0,0,tOld,T,solarVec)
##    dai_dy=partial(j,raObserved,decObserved,radecInitial,r,rdot,1,0,tOld,T,solarVec)
##    dai_dz=partial(j,raObserved,decObserved,radecInitial,r,rdot,2,0,tOld,T,solarVec)
##    dai_dxdot=partial(j,raObserved,decObserved,radecInitial,r,rdot,3,0,tOld,T,solarVec)
##    dai_dydot=partial(j,raObserved,decObserved,radecInitial,r,rdot,4,0,tOld,T,solarVec)
##    dai_dzdot=partial(j,raObserved,decObserved,radecInitial,r,rdot,5,0,tOld,T,solarVec)
##    
##    ddi_dx=partial(j,raObserved,decObserved,radecInitial,r,rdot,0,1,tOld,T,solarVec)
##    ddi_dy=partial(j,raObserved,decObserved,radecInitial,r,rdot,1,1,tOld,T,solarVec)
##    ddi_dz=partial(j,raObserved,decObserved,radecInitial,r,rdot,2,1,tOld,T,solarVec)
##    ddi_dxdot=partial(j,raObserved,decObserved,radecInitial,r,rdot,3,1,tOld,T,solarVec)
##    ddi_dydot=partial(j,raObserved,decObserved,radecInitial,r,rdot,4,1,tOld,T,solarVec)
##    ddi_dzdot=partial(j,raObserved,decObserved,radecInitial,r,rdot,5,1,tOld,T,solarVec)
##    return dai_dx[0]*dai_dx[1]+dai_dy[0]*dai_dy[1]+dai_dz[0]*dai_dz[1]+dai_dxdot[0]*dai_dxdot[1]+dai_dydot[0]*dai_dydot[1]+dai_dzdot[0]*dai_dzdot[1]+ddi_dx[0]*ddi_dx[1]+ddi_dy[0]*ddi_dy[1]+ddi_dz[0]*ddi_dz[1]+ddi_dxdot[0]*ddi_dxdot[1]+ddi_dydot[0]*ddi_dydot[1]+ddi_dzdot[0]*ddi_dzdot[1],0
def deltAi(j,raObserved,decObserved,radecInitial,r,rdot,repectTo,raordec,tOld,T,solarVec):
    if raordec==0:
        return raObserved[j]-radecInitial[j,0],0
    else:
        return decObserved[j]-radecInitial[j,1],0
    
def summation(raObserved,decObserved,radecInitial,r,rdot,tMid,t,solar,function1,function2,respectTo1,respectTo2,nObs):
    sum=0
    for j in range (0,nObs):
        sum=sum+function1(j,raObserved,decObserved,radecInitial,r,rdot,respectTo1,0,tMid,t[j],solar[j])[0]*function2(j,raObserved,decObserved,radecInitial,r,rdot,respectTo2,0,tMid,t[j],solar[j])[0]+function1(j,raObserved,decObserved,radecInitial,r,rdot,respectTo1,1,tMid,t[j],solar[j])[0]*function2(j,raObserved,decObserved,radecInitial,r,rdot,respectTo2,1,tMid,t[j],solar[j])[0]
    return sum

def diffCorrection(raObserved,decObserved,radecInitial,r,rdot,tMid,t,solar,nObs):
    matrix1=np.array([[summation(raObserved,decObserved,radecInitial,r,rdot,tMid,t,solar,deltAi,partial,0,0,nObs)],[summation(raObserved,decObserved,radecInitial,r,rdot,tMid,t,solar,deltAi,partial,0,1,nObs)],[summation(raObserved,decObserved,radecInitial,r,rdot,tMid,t,solar,deltAi,partial,0,2,nObs)],[summation(raObserved,decObserved,radecInitial,r,rdot,tMid,t,solar,deltAi,partial,0,3,nObs)],[summation(raObserved,decObserved,radecInitial,r,rdot,tMid,t,solar,deltAi,partial,0,4,nObs)],[summation(raObserved,decObserved,radecInitial,r,rdot,tMid,t,solar,deltAi,partial,0,5,nObs)]])
    matrix2=np.empty([6,6])
    
    for row in range (0,6):
        matrix2[row,:]=np.array([summation(raObserved,decObserved,radecInitial,r,rdot,tMid,t,solar,partial,partial,row,0,nObs),summation(raObserved,decObserved,radecInitial,r,rdot,tMid,t,solar,partial,partial,row,1,nObs),summation(raObserved,decObserved,radecInitial,r,rdot,tMid,t,solar,partial,partial,row,2,nObs),summation(raObserved,decObserved,radecInitial,r,rdot,tMid,t,solar,partial,partial,row,3,nObs),summation(raObserved,decObserved,radecInitial,r,rdot,tMid,t,solar,partial,partial,row,4,nObs),summation(raObserved,decObserved,radecInitial,r,rdot,tMid,t,solar,partial,partial,row,5,nObs)])
    matrix2=np.linalg.inv(matrix2)
    changeR=matrix2.dot(matrix1)
    return changeR



#*************************************TEST CASES:******************************************
def run(inFileName,dataPerObs):
    inFile=open(inFileName,"r")
    fileList=inFile.read().split()
    #numObs is the number of data points in the file while dataperObs is the number of data points per observation.
    numObs=int(len(fileList)/7)
    obsList=np.empty([numObs,6])
    obsCount=0
    #This for loop populates obsList which is an numpy array with each row for each data point. Columns are JD,ra,dec,and solar vector values
    for j in range (0,numObs):
        obs=fileList[obsCount:obsCount+7]
        JD=JDNumber(obs[0],obs[1])
        ra=decimalCoord(obs[2])
        dec=decimalCoord(obs[3])
        solar=np.array([float(obs[4]),float(obs[5]),float(obs[6])])
        obsList[j,:]=np.array([JD,ra,dec,solar[0],solar[1],solar[2]])
        obsCount=obsCount+7

    rangeList=np.empty([numObs-2*dataPerObs,1])
    rVectors=np.empty([numObs-2*dataPerObs,6])
    oElemList=np.empty([numObs-2*dataPerObs,6])
    radecList=np.empty([numObs-2*dataPerObs,6])
    OminusCList=np.empty([numObs-2*dataPerObs])
    JDMiddles=np.empty([numObs-2*dataPerObs])
    #rangeList is a list of the range for each middle observation
    #rVectors is a list of all r and rdots for each middle observation
    #oElemList is a list of all orbital elements for each middle observation
    #radecList is a calculated list of all ra dec values for all observations
    #OminusCList is a list of all simplified OminusC values
    #JDMiddles is a collection of all middle JD dates
    for j in range (0,numObs-2*dataPerObs):
        JDFirst=obsList[0,0]
        raFirst=obsList[0,1]
        decFirst=obsList[0,2]
        solarFirst=obsList[0,3:6]
        
        JDLast=obsList[numObs-1,0]
        raLast=obsList[numObs-1,1]
        decLast=obsList[numObs-1,2]
        solarLast=obsList[numObs-1,3:6]
        
        JDMiddle=obsList[j+dataPerObs,0]
        raMiddle=obsList[j+dataPerObs,1]
        decMiddle=obsList[j+dataPerObs,2]
        solarMiddle=obsList[j+dataPerObs,3:6]
        JDMiddles[j]=JDMiddle
        
        randrDot=methodOfGauss(JDFirst,raFirst,decFirst,solarFirst,JDMiddle,raMiddle,decMiddle,solarMiddle,JDLast,raLast,decLast,solarLast)
        rangeList[j,0]=randrDot[2]
        rVectors[j,:]=np.array([randrDot[0][0],randrDot[0][1],randrDot[0][2],randrDot[1][0],randrDot[1][1],randrDot[1][2]])
        r=randrDot[0]
        rdot=randrDot[1]
        oElem=orbElems(r,rdot)
        oElemList[j,:]=oElem
        radec1=generateRaDec(oElem[0],oElem[1],oElem[3],oElem[2],oElem[4],oElem[5],JDMiddle,JDFirst,solarFirst)
        radec2=generateRaDec(oElem[0],oElem[1],oElem[3],oElem[2],oElem[4],oElem[5],JDMiddle,JDMiddle,solarMiddle)
        radec3=generateRaDec(oElem[0],oElem[1],oElem[3],oElem[2],oElem[4],oElem[5],JDMiddle,JDLast,solarLast)
        radecList[j,:]=[radec1[0]/15,radec1[1],radec2[0]/15,radec2[1],radec3[0]/15,radec3[1]]
        OminusCList[j]=(OminusC1([raFirst,raMiddle,raLast],[decFirst,decMiddle,decLast],[radec1,radec2,radec3]))

    choice=np.argmin(OminusCList)
    #randrDot, elements,ra/dec, and range for chosen middle observation
    randrDotFinal=rVectors[choice,:]
    elementsFinal=oElemList[choice,:]
    radecFinal=radecList[choice,:]
    rangeFinal=rangeList[choice,0]

    differentialYes=True
    if differentialYes:
        #Below are lists for JD,solarvector,ra,dec for all data points
        JDAll=obsList[:,0]
        solarAll=obsList[:,3:6]
        raObserved=obsList[:,1]
        decObserved=obsList[:,2]
        radecInitial=np.empty([numObs,2])#Initial ra dec calculations
        for j in range(0,numObs):
            radecInitial[j,:]=generateRaDec(elementsFinal[0],elementsFinal[1],elementsFinal[3],elementsFinal[2],elementsFinal[4],elementsFinal[5],JDMiddles[choice],JDAll[j],solarAll[j])
        #OminusC initial begins here:
        OminusCInitial=OminusC1(raObserved,decObserved,radecInitial)
        #print OminusCInitial

        #differential correction
        changeVectors=diffCorrection(raObserved*15,decObserved,radecInitial,randrDotFinal[0:3],randrDotFinal[3:6],JDMiddles[choice],JDAll,solarAll,numObs)
        changeVectors=np.hstack(changeVectors)
        randrDotFinal=np.array([randrDotFinal[0]+changeVectors[0],randrDotFinal[1]+changeVectors[1],randrDotFinal[2]+changeVectors[2],randrDotFinal[3]+changeVectors[3],randrDotFinal[4]+changeVectors[4],randrDotFinal[5]+changeVectors[5]])
        elementsFinal=np.array(orbElems(np.array(randrDotFinal[0:3]),np.array(randrDotFinal[3:6])))
        JDFirst=obsList[0,0]
        JDLast=obsList[numObs-1,0]
        JDMiddle=JDMiddles[choice]
        solarFirst=solarAll[0]
        solarLast=solarAll[numObs-1]
        solarMiddle=solarAll[choice+dataPerObs]
        radec1=generateRaDec(elementsFinal[0],elementsFinal[1],elementsFinal[3],elementsFinal[2],elementsFinal[4],elementsFinal[5],JDMiddles[choice],JDFirst,solarFirst)
        radec2=generateRaDec(elementsFinal[0],elementsFinal[1],elementsFinal[3],elementsFinal[2],elementsFinal[4],elementsFinal[5],JDMiddles[choice],JDMiddle,solarMiddle)
        radec3=generateRaDec(elementsFinal[0],elementsFinal[1],elementsFinal[3],elementsFinal[2],elementsFinal[4],elementsFinal[5],JDMiddles[choice],JDLast,solarLast)
        radecFinal=np.array([[radec1[0]/15,radec1[1]],[radec2[0]/15,radec2[1]],[radec3[0]/15,radec3[1]]])
        radecFinalCalc=np.empty([numObs,2])#Final ra dec calculations
        for j in range(0,numObs):
            radecFinalCalc[j,:]=generateRaDec(elementsFinal[0],elementsFinal[1],elementsFinal[3],elementsFinal[2],elementsFinal[4],elementsFinal[5],JDMiddles[choice],JDAll[j],solarAll[j])
        #OminusC Final begins here:
        OminusCLast=OminusC1(raObserved,decObserved,radecFinalCalc)
    print "\n\n***********************************\nCentral Observation: "+str(fileList[7*(choice+2):7*(choice+2)+2])
    print "\n***********************************\nFINAL Position and Velocity vectors for central observation:\n"
    print randrDotFinal[0:3],randrDotFinal[3:6]
    print"\n************************************\nRange to the asteroid for central observation in AU:\n"
    print rangeFinal
    print"\n************************************\nORBITAL ELEMENTS IN ORDER a,e,O,I,w,M:\n"
    print elementsFinal
    print "\n\n************************************\nRA (in hours) AND DEC (in degrees) FOR First, Second, Third Observation:\n"
    print radecFinal[0:2]
    print radecFinal[2:4]
    print radecFinal[4:6]
    if differentialYes:
        print "\n\n************************************\nOminusC before and after differential correction:\n"
        print OminusCInitial,OminusCLast
    vPythSim(elementsFinal[0],elementsFinal[1],elementsFinal[2],elementsFinal[3],elementsFinal[4],elementsFinal[5])
    return True
    
#run("1994PN Test Case.txt")
#run("Gauss_Input.txt",1)
#run("2002 KL6Final.txt",2)
#run("ZhengdongWangInput.txt",1)
#run("1999ML Test Case.txt",1)

fileName=raw_input("Enter file name: ")
dataPerObss=int(raw_input("Enter number of data points per observation: "))
run(fileName,dataPerObss)
    
    
    
    
