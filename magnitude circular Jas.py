from __future__ import division
from math import*
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

#The Following Code is for Centroiding:
#returns sum of every element in array
def N(array):
    sum=0
    for rowIndex in range(0,len(array)):
        for columnIndex in range(0,len(array[0])):
            sum=sum+array[rowIndex,columnIndex]
    return sum

#returns array of sums of each column
def sumColumn(array):
    sumColumn=[]
    for columnIndex in range(0,len(array[0])):
        sum=0
        for rowIndex in range(0,len(array)):
            sum=sum+array[rowIndex,columnIndex]
        sumColumn.append(sum)
    return sumColumn
    
#returns array of sums of each row
def sumRow(array):
    sumRow=[]
    for rowIndex in range(0,len(array)):
        sum=0
        for columnIndex in range(0,len(array[0])):
            sum=sum+array[rowIndex,columnIndex]
        sumRow.append(sum)
    return sumRow

#returns x position for centroid
def xPosition(array):
    sc = sumColumn(array)
    sumBeforeMean=0
    divideSum=0
    
    for scIndex in range(0,len(sc)):
        sumBeforeMean=sumBeforeMean+scIndex*sc[scIndex]
        divideSum=divideSum+sc[scIndex]

    return sumBeforeMean/divideSum

#returns y position for centroid
def yPosition(array):
    sr=sumRow(array)
    sumBeforeMean=0
    divideSum=0
    for srIndex in range(0,len(sr)):
        sumBeforeMean=sumBeforeMean+srIndex*sr[srIndex]
        divideSum=divideSum+sr[srIndex]
    return sumBeforeMean/divideSum

#returns x and y coordinate for centroid
def centroid (array):
    return xPosition(array),yPosition(array)

#returns array of distance squared of each column to the x centroid coordinate
def distX(array):
    distanceX=[]
    xPos=xPosition(array)
    for column in range(0,len(array[0])):
        distanceX.append((xPos-column)**2)
    return distanceX

#returns array of distance squared of each column to the x centroid coordinate
def distY(array):
    distanceY=[]
    yPos=yPosition(array)
    for row in range(0,len(array)):
        distanceY.append((yPos-row)**2)
    return distanceY

#finds stdev for x coordinate
def stdevX(array):
    dX=distX(array)
    sCol=sumColumn(array)
    summation=0.0
    for col in range(0,len(array[0])):
        summation=summation+dX[col]*sCol[col]
    N1=N(array)/1000
    N2=(N(array)-1)/1000
    result=(1/1000)*sqrt(summation/(N1*N2))
    return result

#finds stdev for y coordinate
def stdevY(array):
    dY=distY(array)
    sRow=sumRow(array)
    summation=0
    for row in range(0,len(array)):
        summation=summation+dY[row]*sRow[row]
    N1=N(array)/1000
    N2=(N(array)-1)/1000
    result=(1/1000)*sqrt(summation/(N1*N2))
    return result

#return centroid, uncertainties of image
def centroidIm(picName, width, height, x, y):
    image=fits.getdata(picName)
    cutImage=image[y-int(height/2):y+int(height/2)+1,x-int(width/2):x+int(width/2)+1]
    uncX=stdevX(cutImage)
    uncY=stdevY(cutImage)
    cent=centroid(cutImage)
    newCent=[cent[0]+x-int(width/2),cent[1]+y-int(height/2)]
    return newCent,uncX,uncY










#The Following Code is for Photometry
#Returns array describing percent of pixel to be used in circular aperature
def percentArray(r):
    arr=np.zeros([2*r+1,2*r+1])
    for row in range(0,2*r+1):
        for col in range(0,2*r+1):
            dist=sqrt((r-col)**2+(r-row)**2)
            if dist>=r and dist<=r+1:
                arr[row,col]=1-(dist-r)
            elif dist<r:
                arr[row,col]=1
    return arr

#Returns array describing percent of pixel to be used in annulus
def annulusPercent(inR,outR):
    ann=np.zeros([2*outR+1,2*outR+1])
    for row in range(0,2*outR+1):
        for col in range(0,2*outR+1):
            dist=sqrt((outR-col)**2+(outR-row)**2)
            if dist>=outR and dist<=outR+1:
                ann[row,col]=1-(dist-outR)
            elif dist<outR and dist>inR:
                ann[row,col]=1
            elif dist<=inR and dist>=inR-1:
                ann[row,col]=1-(inR-dist)
    return ann

#Returns signal and skyAverage
def signal(posx,posy,r,inR,outR,imArray):
    aperaturePercent=percentArray(r)
    aperature=imArray[posy-r:posy+r+1,posx-r:posx+r+1]
    starSky=0
    for row in range(0,len(aperature)):
        for col in range(0,len(aperature[0])):
            starSky=starSky+aperaturePercent[row,col]*aperature[row,col]

    annPercent=annulusPercent(inR,outR)
    annulus=imArray[posy-outR:posy+outR+1,posx-outR:posx+outR+1]
    sky=0
    for row in range(0,len(annulus)):
        for col in range(0,len(annulus[0])):
            sky=sky+annulus[row,col]*annPercent[row,col]
    nAnn=np.sum(annPercent)
    skyAv=sky/nAnn
    nAp=np.sum(aperaturePercent)
    return starSky-(skyAv*nAp),skyAv

#Returns constant
def constant(starx,stary,r,inR,outR,imArray,starMag):
    sig=signal(starx,stary,r,inR,outR,imArray)[0]
    return starMag+2.5*log10(sig)





#Returns asteroid magnitude and sky brightness. Inputs are star x position, star y position, radius of aperature, inner annulus radius,
#outer annulus radius, image name in a string, magnitude of star, asteroid x position, asteroid y position
def astMag(starx,stary,r,inR,outR,picName,starMag,astx,asty):
    im=fits.getdata(picName)
    c=constant(starx,stary,r,inR,outR,im,starMag)
    sig=signal(astx,asty,r,inR,outR,im)
    skyMag=-2.5*log10(sig[1])+c
    return -2.5*log10(sig[0])+c, skyMag

#Returns sky brightness in magnitude per arcsecond squared. Arguments are focal length, pixel size (where both are in the same units), and the magnitude of the sky in a pixel
def skyBright(fLength,pixelSize,skyMag):
    plateScale=206265/fLength*pixelSize
    return skyMag/(plateScale**2)

#Remember centroidIm takes the name of the image, width and heigth, and then approximate coordinates of center

#Run takes picture name in string, output file name, star positions in x,y list, star magnitudes, ateroid position in x,y, focal length, and pixel size
def run(picName, outputFile,star1,star2,star3,star4,star5,mag1,mag2,mag3,mag4,mag5,astPos,fLength,pixelSize):
    newfile=open(outputFile,"w")
    centroid=centroidIm(picName,7,7,astPos[0],astPos[1])
    xy=centroid[0]
    newfile.write("Asteroid Centroid X,Y: "+str(xy))
    
    amag1=astMag(star1[0],star1[1],5,7,11,picName,mag1,int(round(xy[0])),int(round(xy[1])))
    sky1=skyBright(fLength,pixelSize,amag1[1])
    newfile.write("\n\nAsteroid Magnitude with Star 1: "+str(amag1[0]))
    newfile.write("\nSky Brightness in mag/arcseconds squared: "+str(sky1))
    
    amag2=astMag(star2[0],star2[1],5,7,11,picName,mag2,int(round(xy[0])),int(round(xy[1])))
    sky2=skyBright(fLength,pixelSize,amag2[1])
    newfile.write("\n\nAsteroid Magnitude with Star 2: "+str(amag2[0]))
    newfile.write("\nSky Brightness in mag/arcseconds squared: "+str(sky2))
    
    amag3=astMag(star3[0],star3[1],5,7,11,picName,mag3,int(round(xy[0])),int(round(xy[1])))
    sky3=skyBright(fLength,pixelSize,amag3[1])
    newfile.write("\n\nAsteroid Magnitude with Star 3: "+str(amag3[0]))
    newfile.write("\nSky Brightness in mag/arcseconds squared: "+str(sky3))

    amag4=astMag(star4[0],star4[1],5,7,11,picName,mag4,int(round(xy[0])),int(round(xy[1])))
    sky4=skyBright(fLength,pixelSize,amag4[1])
    newfile.write("\n\nAsteroid Magnitude with Star 4: "+str(amag4[0]))
    newfile.write("\nSky Brightness in mag/arcseconds squared: "+str(sky4))

    amag5=astMag(star5[0],star5[1],5,7,11,picName,mag5,int(round(xy[0])),int(round(xy[1])))
    sky5=skyBright(fLength,pixelSize,amag5[1])
    newfile.write("\n\nAsteroid Magnitude with Star 5: "+str(amag5[0]))
    newfile.write("\nSky Brightness in mag/arcseconds squared: "+str(sky5))

    avgAst=(amag1[0]+amag2[0]+amag3[0]+amag4[0]+amag5[0])/5
    avgSkyMag=(sky1+sky2+sky3+sky4+sky5)/5
    newfile.write("\n\nAverage asteroid magnitude: "+str(avgAst))
    newfile.write("\nAverage sky magnitude: "+str(avgSkyMag))
    newfile.close()
    return True

#run("Obs1.fit","photometry1.txt",[773,453],[592,299],[539,237],[487,499],[624,432],14.64,15.56,13.5,13.08,15.53,[500,454],2346,0.0068)
#run("Obs3SET2.00000138.fit","photometry3.txt",[704,249],[668,473],[521,621],[477,353],[598,209],14.58,12.33,15.32,14.17,13.44,[502,514],2346,0.0068)
#run("Obs1(82).fit","photometry1(82).txt",[358,477],[396,453],[699,361],[487,499],[624,432],11.52,13.32,13.56,13.08,15.88,[484,419],2346,0.0068)
#run("DARK.00000092.ENTERED COORDINATES.REDUCED.fit","photometry2(92).txt",[493,420],[460,309],[599,157],[571,593],[489,250],15.54,16.15,16.13,16.10,15.87,[498,527],2346,0.0068)
#run("Obs3set3(140).fit","photometry3(140).txt",[479,355],[597,405],[593,487],[571,302],[232,239],13.44,12.75,14.11,15.35,15.11,[502,514],2346,0.0068)

loop=True
while loop==True:
    fileName=raw_input("Enter name of file: ")
    outFile=raw_input("Enter name of write file: ")
    x1=float(raw_input("Enter x pixel for star 1: "))
    y1=float(raw_input("Enter y pixel for star 1: "))
    mag1=float(raw_input("Enter magnitude for star 1: "))
    x2=float(raw_input("Enter x pixel for star 2: "))
    y2=float(raw_input("Enter y pixel for star 2: "))
    mag2=float(raw_input("Enter magnitude for star 2: "))
    x3=float(raw_input("Enter x pixel for star 3: "))
    y3=float(raw_input("Enter y pixel for star 3: "))
    mag3=float(raw_input("Enter magnitude for star 3: "))
    x4=float(raw_input("Enter x pixel for star 4: "))
    y4=float(raw_input("Enter y pixel for star 4: "))
    mag4=float(raw_input("Enter magnitude for star 4: "))
    x5=float(raw_input("Enter x pixel for star 5: "))
    y5=float(raw_input("Enter y pixel for star 5: "))
    mag5=float(raw_input("Enter magnitude for star 5: "))
    astx=float(raw_input("Enter x pixel for object: "))
    asty=float(raw_input("Enter y pixel for object: "))
    fLength=float(raw_input("Enter focal length: "))
    pixelSize=float(raw_input("Enter pixel size: "))
    run(fileName,outFile,[x1,y1],[x2,y2],[x3,y3],[x4,y4],[x5,y5],mag1,mag2,mag3,mag4,mag5,[astx,asty],fLength,pixelSize)

    choice=int(raw_input("Enter 1 to make another calculation or 2 to end: "))
    if choice==2:
        loop=False
    
