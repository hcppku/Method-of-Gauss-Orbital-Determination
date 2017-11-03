from __future__ import division
from math import*
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
   
#Formats file into 2D array with each row being x,y centroid, ra, and then dec for each star
def starArray(listStar):
    starInfo=[]
    for index in range(0,int(len(listStar)/4)):
        starLine=[]
        for smIndex in range(0,4):
            if smIndex==0 or smIndex==1:
                starLine.append(float(listStar[4*index+smIndex]))
            if smIndex==2:
                starLine.append(15*raDecFormat(listStar,index,smIndex))
            if smIndex==3:
                starLine.append(raDecFormat(listStar,index,smIndex))
        starInfo.append(starLine)
    starInfo=np.array(starInfo)
    return starInfo

#Turns ra or dec information from file into decimal value
def raDecFormat(listStar,index,smIndex):
    hrsdeg=""
    minutes=""
    raOrDec=listStar[4*index+smIndex]
    strIndex=0
    char=raOrDec[strIndex]
    while char!=":":
        hrsdeg=hrsdeg+char
        strIndex=strIndex+1
        char=raOrDec[strIndex]
    strIndex=strIndex+1
    char=raOrDec[strIndex]
    while char!=":":
        minutes=minutes+char
        strIndex=strIndex+1
        char=raOrDec[strIndex]
    seconds=raOrDec[strIndex+1:]
    return pi/180*(float(hrsdeg)+float(minutes)/60+float(seconds)/3600)
            
#returns plate constants for RA        
def chiRA(starInfo):
    sumRA=0
    for row in range(0,len(starInfo)):
        sumRA=sumRA+starInfo[row][2]
    sumRAX=0
    for row in range(0,len(starInfo)):
        sumRAX=sumRAX+starInfo[row][0]*starInfo[row][2]
    sumRAY=0
    for row in range(0,len(starInfo)):
        sumRAY=sumRAY+starInfo[row][1]*starInfo[row][2]
    sumX=0
    for row in range(0,len(starInfo)):
        sumX=sumX+starInfo[row][0]
    sumY=0
    for row in range(0,len(starInfo)):
        sumY=sumY+starInfo[row][1]
    sumXY=0
    for row in range(0,len(starInfo)):
        sumXY=sumXY+starInfo[row][0]*starInfo[row][1]
    sumXX=0
    for row in range(0,len(starInfo)):
        sumXX=sumXX+starInfo[row][0]**2
    sumYY=0
    for row in range(0,len(starInfo)):
        sumYY=sumYY+starInfo[row][1]**2
    arr1=np.array([[sumRA],[sumRAX],[sumRAY]])
    arr2=np.array([[len(starInfo),sumX,sumY],[sumX,sumXX,sumXY],[sumY,sumXY,sumYY]])
    invarr2=np.linalg.inv(arr2)
    return invarr2.dot(arr1)

#returns plate constants for Dec
def chiDec(starInfo):
    sumD=0
    for row in range(0,len(starInfo)):
        sumD=sumD+starInfo[row][3]
    sumDX=0
    for row in range(0,len(starInfo)):
        sumDX=sumDX+starInfo[row][0]*starInfo[row][3]
    sumDY=0
    for row in range(0,len(starInfo)):
        sumDY=sumDY+starInfo[row][1]*starInfo[row][3]
    sumX=0
    for row in range(0,len(starInfo)):
        sumX=sumX+starInfo[row][0]
    sumY=0
    for row in range(0,len(starInfo)):
        sumY=sumY+starInfo[row][1]
    sumXY=0
    for row in range(0,len(starInfo)):
        sumXY=sumXY+starInfo[row][0]*starInfo[row][1]
    sumXX=0
    for row in range(0,len(starInfo)):
        sumXX=sumXX+starInfo[row][0]**2
    sumYY=0
    for row in range(0,len(starInfo)):
        sumYY=sumYY+starInfo[row][1]**2
    arr1=np.array([[sumD],[sumDX],[sumDY]])
    arr2=np.array([[len(starInfo),sumX,sumY],[sumX,sumXX,sumXY],[sumY,sumXY,sumYY]])
    invarr2=np.linalg.inv(arr2)
    return invarr2.dot(arr1)

#returns array of plate constants in order: b1, b2, a11, a12, a21, a22
def plateConstants(decConstants,raConstants):
    return np.array([raConstants[0][0],decConstants[0][0],raConstants[1][0],raConstants[2][0],decConstants[1][0],decConstants[2][0]])

#returns ra in hours and dec in degrees for plate constants and x,y inputs
def findRADec(pltConstants,x,y):
    ra=pltConstants[0]+pltConstants[2]*x+pltConstants[3]*y
    dec=pltConstants[1]+pltConstants[4]*x+pltConstants[5]*y
    return (ra*180/pi)/15,(dec*180/pi)

#returns uncertainty for ra and dec in arcseconds
def uncRaDec(pltConstants,starInfo):
    sumRA=0
    sumDec=0
    for row in range(0,len(starInfo)):
        raDec=findRADec(pltConstants,starInfo[row][0],starInfo[row][1])
        sumRA=sumRA+(starInfo[row][2]-raDec[0]*15*pi/180)**2
        sumDec=sumDec+(starInfo[row][3]-raDec[1]*pi/180)**2
    uncRa=sqrt((1/(len(starInfo)-3))*sumRA)*180*3600/pi
    uncDec=sqrt((1/(len(starInfo)-3))*sumDec)*180*3600/pi
    return uncRa,uncDec

#formats ra and dec for output
def formatOutput(radec):
    ra=radec[0]
    dec=radec[1]
    minusPlus="+"
    if dec<0:
        minusPlus="-"
    hours=int(ra)
    minutes=int((ra-hours)*60)
    seconds=round((ra-hours-minutes/60)*3600,2)
    degrees=int(dec)
    arcminutes=int((dec-degrees)*60)
    arcseconds=round((dec-degrees-arcminutes/60)*3600,2)
    if hours<10:
        hours="0"+str(hours)
    else:
        hours=str(hours)
    if minutes<10:
        minutes="0"+str(minutes)
    else:
        minutes=str(minutes)
    if seconds<10:
        seconds="0"+str(seconds)
    else:
        seconds=str(seconds)
    if degrees<10:
        degrees="0"+str(degrees)
    else:
        degrees=str(degrees)
    if arcminutes<10:
        arcminutes="0"+str(arcminutes)
    else:
        arcminutes=str(arcminutes)
    if arcseconds<10:
        arcseconds="0"+str(arcseconds)
    else:
        arcseconds=str(arcseconds)
    
    return hours+":"+minutes+":"+seconds, minusPlus+degrees+":"+arcminutes+":"+arcseconds

#Runs and outputs file:
def run(infileName,posX,posY,outputFile):
    starsFile=open(infileName,"r")
    starsList=starsFile.read().split()
    starsFinal=starArray(starsList)
    plts=plateConstants(chiDec(starsFinal),chiRA(starsFinal))
    unc=uncRaDec(plts,starsFinal)
    fRADec=formatOutput(findRADec(plts,posX,posY))
    newfile=open(outputFile,"w")
    newfile.write("***************")
    newfile.write("\nplate constants")
    newfile.write("\n***************")
    newfile.write("\nb1: "+str(plts[0]))
    newfile.write("\nb2: "+str(plts[1]))
    newfile.write("\na11: "+str(plts[2]))
    newfile.write("\na12: "+str(plts[3]))
    newfile.write("\na21: "+str(plts[4]))
    newfile.write("\na22: "+str(plts[5]))
    newfile.write("\n***********")
    newfile.write("\nuncertainty")
    newfile.write("\n***********")
    newfile.write("\nRA: "+str(unc[0])+" arcsec")
    newfile.write("\nDec: "+str(unc[1])+" arcsec")
    newfile.write("\n***********************************")
    newfile.write("\nastrometry for (x,y)=("+str(posX)+".,"+str(posY)+".)")
    newfile.write("\n***********************************")
    newfile.write("\nRA= "+fRADec[0])
    newfile.write("\nDec= "+fRADec[1])
    newfile.close()
    return True

#Modified Header to incorporate RA and dec
def setHeader(infileName,imageName):
    
    starsFile=open(infileName,"r")
    starsList=starsFile.read().split()
    starsFinal=starArray(starsList)
    plts=plateConstants(chiDec(starsFinal),chiRA(starsFinal))
    image=fits.getdata(imageName)
    hdulist=fits.open(imageName,mode="update")
    hdu=hdulist[0].header
    hdu.set("CTYPE1","RA---TAN")
    hdu.set("CRPIX1",len(image[0])/2)
    hdu.set("CRVAL1",plts[0])
    hdu.set("CTYPE2","DEC--TAN")
    hdu.set("CRPIX2",len(image)/2)
    hdu.set("CRVAL2",plts[1])
    hdu.set("CD1_1",plts[2])
    hdu.set("CD1_2",plts[3])
    hdu.set("CD2_1",plts[4])
    hdu.set("CD2_2",plts[5])
    hdu.set("RADECSYS","FK5")
    hdu.set("EQUINOX",2000.0)
    hdu.set("MJD-OBS",0.0)
    hdulist.close()
    return True

#Test:
#run("LSPRparallax.txt",477.80701972853257,320.0889829754683,"lsprOutputparallax.txt")
#run("LSPRstars2.txt",498.00644650222625,526.99986342156296,"lsprOutput2(92).txt")
#run("LSPRstars1.txt",502.16667658789214,458.78689207690934,"lsprOutput1(59).txt")
#setHeader("LSPRstars1.txt","Obs1.fit")

loop=True
while loop==True:
    fileName=raw_input("Enter input file: ")
    outPut=raw_input("Enter write file: ")
    x=float(raw_input("Enter x centroid: "))
    y=float(raw_input("Enter y centroid: "))
    run(fileName,x,y,outPut)

    choice=int(raw_input("Enter 1 to do another calculation or 2 to end: "))
    if choice==2:
        loop=False
        



