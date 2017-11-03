from __future__ import division
import math
#Sidereal Calculator

def toDecimalHour(hour, minute, seconds):
    return hour+minute/60.0+seconds/3600.0

def julianDay(day, month, year):
    return 367*year-int(7*(year+int((month+9)/12))/4)+int(275*month/9)+day+1721013.5

def julianCenturies(julianDay):
    return (julianDay-2451545.0)/36525

def gmstMean(julianCenturies):
    return 100.46061837+36000.77053608*julianCenturies+3.87933*(10**-4)*(julianCenturies**2)-(julianCenturies**3)/(3.871*(10**7))

def gmstActual(gmstMean,ut):
    return gmstMean+360.985647366*(ut/24)

def locSidRaw(gmstActual,longitude, westOrEast):
    if westOrEast=="east":
        return gmstActual+longitude
    if westOrEast=="west":
        return gmstActual-longitude

def locSidHours(locSidRaw):
    while locSidRaw>360:
        locSidRaw=locSidRaw-360
    return locSidRaw/15

def convertToUt(localDecimal,z):
    return localDecimal+z


def run():
    loop=True
    while loop==True:
        print "1: decimal hour"
        print "2: UT time in hours"
        print "3: julian day"
        print "4: julian centuries"
        print "5: julian day number (incorporates time)"
        print "6: GM Sidereal Time Mean"
        print "7: GM Actual Sidereal Time"
        print "8: Raw Local Sidereal Time (in degrees)"
        print "9: Local Sidereal Time (in hours)"
        choice=int(raw_input("Enter a number to select something to calculate: "))
        if choice>9 or choice<1:
            print "That choice is not available. Try again."
        else:
            loop=False
    
    hour=float(raw_input("Enter in the hour: "))
    minutes=float(raw_input("Enter in the minutes: "))
    seconds=float(raw_input("Enter in the seconds: "))
    decHour=toDecimalHour(hour,minutes,seconds)
    print "The decimal hour value is: "+str(decHour)
    if choice<2:
        return True
    z=float(raw_input("Enter in the time zone shift from UT (eg 6 for etscorn daylights savings)"))
    UT=convertToUt(decHour,z)
    print "UT time is: "+str(UT)
    if choice<3:
        return True
    day=int(raw_input("Enter day number: "))
    month=int(raw_input("Enter month number: "))
    year=int(raw_input("Enter year number: "))
    jDay=julianDay(day,month,year)
    print "julian day is: "+str(jDay)
    if choice<4:
        return True
    jCentury=julianCenturies(jDay)
    print "julian century is: "+str(jCentury)
    if choice<5:
        return True
    JD=jDay+UT/24
    print "julian day number is: "+str(JD)
    if choice<6:
        return True
    gmSidMean=gmstMean(jCentury)
    print "GM Sidereal Time Mean is: "+str(gmSidMean)
    if choice<7:
        return True
    gmSidActual=gmstActual(gmSidMean,UT)
    print "GM Sidereal Time Actual is: "+str(gmSidActual)
    if choice<8:
        return True
    longitude=float(raw_input("Enter positive longitude: "))
    westEast=raw_input("Enter east or west: ")
    rawLocalSid=locSidRaw(gmSidActual,longitude,westEast)
    print "Raw local sidereal time: "+str(rawLocalSid)
    if choice<9:
        return True
    localSid=locSidHours(rawLocalSid)
    print "Local Sidereal time: "+str(localSid)


program=True
while program==True:
    run()
    endProgram=int(raw_input("Enter 1 to make another calculation or 2 to end program: "))
    if endProgram==2:
        program=False
        



