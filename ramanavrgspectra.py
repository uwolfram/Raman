#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on May 2025

@author: Uwe Wolfram

Code evaluates n individual Raman spectra taken from bone samples and averages the resulting fits 


Code was originally developed 2015 for: 
    [1] Mirzaali, M., Schwiedrzik, J., Thaiwichai, S., Best, J., Michler, J., Zysset, P. & Wolfram, U. 
        Mechanical properties of cortical bone and their relationships with age, gender, composition and 
        microindentation properties in the elderly. Bone 93, 196–211 (2016).

Spectra are evaluated for:
    v1PO4
    v2PO4
    Amide1
    Amide3
    PYD
    Lipids
    
Boundaries for these spectral regions are hard-encoded at the moment 

Check [1] for baclground info


Input
----------
spectra : 
    Directory holding a set of baseline corrected spectra consisting of only the counts.
    
nbrspec:
    Number of specimens tested.
    
nbrmeas:
    Number of repeated measurements per specimen
    
shifts:
    File containing counts and shifts. Do not ask why this was setup like this.
    

Returns
-------
averaged spectra for each specimen and 
ramanavrgfits-mean.dat && ramanavrgfits-std.dat: 
    2 ascii data files providing the mean and std (where possible) over nbrmeas for
    filename (donor)
    age
    gender
    v1PO4fwhm
    v1PO4int
    v1PO4peak
    v2PO4int
    amid1int
    amid1peak
    amid3int
    pyd
    CH3int
    lipid
    
    
        
Example run
-----------
python3 ramanavrgspectra.py -d Median-subtracted -n 41 -m 10 -s S18_F70L_1_AX_01.txt

"""


from sys import *
import time as time
from math import *
from numpy import *
#from random import random
from string import *
#from time import sleep
import os # os must be imported like that. otherwise open(...) cannot be used. 
import re
import glob
import argparse

#import matplotlib
#matplotlib.use('Agg') # plotting to this interface would avoid an image to be send to X
#import matplotlib.artist as art
import matplotlib.pyplot as plt

from scipy import ndimage
from scipy import signal
from scipy import stats
from scipy import integrate
from scipy.optimize import leastsq # Levenberg-Marquadt Algorithm
#from skimage import filter

# fix star import issue with min
import builtins as buin


class RamanData():
    
    __version__='V_May_2015'
    __author__='U. Wolfram'

    def __init__(self):

        # variables identifying the specimen
        self.files = []
        self.filename = []
        self.age = []
        self.gender = []
        self.shifts = []
        self.counts = []
        self.specimens = []

        self.v1PO4peak = []
        self.v1PO4fwhm = []
        self.v1PO4int  = []
        self.v2PO4int = []
        self.amid1int = []
        self.amid1peak = []
        self.amid3int = []
        self.pyd = []
        self.CH3int = []
        self.lipid = []        
        self.lbindex = []
        self.ubindex = []


        

    def getFilename(self):
        return self.filename

    def getAge(self):
        return self.age
        
    def getGender(self):
        return self.gender

    def getShifts(self):
        return self.shifts

    def getFiles(self):
        return self.files
        
    def getCounts(self):
        return self.counts        

    def getv1PO4(self):
        return self.v1PO4fwhm
    
    def getv2PO4(self):
        return self.v2PO4int

    def getamid1(self):
        return self.amid1int

    def getamid3(self):
        return self.amid3int
        
    def getpyd(self):
        return self.pyd

    def getCH3(self):
        return self.CH3int

    def getlipid(self):
        return self.lipid        

    def getlbindex(self):
        return self.lbindex

    def getubindex(self):
        return self.ubindex

    def breakPoint(self):

        stdout.write("\n ... B R E A K  P O I N T \n")
        stdout.flush()
        stop = "n"
        print("\n     Going on (y/n): ")
        stop = raw_input()            
        if (stop.lower() != "y"):
            exit(0)
        elif (stop.lower() == "y"):
            return        

    def getInputs(self):
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('-d', '--dir',
                            default='./', # i.e. local folder
                            type=str,
                            help='Directory containing spectra.'
                            )   
        parser.add_argument('-n', '--nbrspec',
                            default=1,
                            type=int,
                            help='Number of specimens tested.'
                            )
        parser.add_argument('-m', '--meas',
                            default=1,
                            type=int,
                            help='Number of repeated measurements per sample.'
                            )
        parser.add_argument('-s', '--shifts',
                            default='foo.txt',
                            type=str,
                            help='File to extract the shifts. Can be a sample from the specimens tested.'
                            )
        
        args = parser.parse_args()
        dir = args.dir
        spec = args.nbrspec
        meas = args.meas
        shifts = args.shifts
        
        return dir, spec, meas, shifts


    def readFiles(self, specdir, nbrspec, meas):
        
        stdout.write("\n \t... Read Force Files \n\n"); stdout.flush()        
        
        # listing subfolders in folder parentfolder
        files = []
        for item in os.listdir(specdir):
            files.append(item)

        files.sort()
        k = 0
        for i in range(nbrspec):
            #split = item.split('_')
        
            measurements = []
            for j in range(meas):
                
                if j == 0:
                    split = files[k].split('_')                
                    self.age.append(split[1][1:3])
                    self.gender.append(split[1][0])
                    self.filename.append(split[0] + "_" + split[1] + "_" + split[2] + "_" + split[3])
                
                measurements.append(files[k])
                k+=1
                        
            self.files.append(measurements)
    
    
    def readShifts(self, filename):
        
        stdout.write("\n \t... Read Counts \n\n"); stdout.flush()
        
        skip = 0
        for line in open(filename,'r'):
            line = re.sub(' +','',line)
            line = line.replace(",", ".")
            line = line.replace("\n", "")
            line = line.replace("\r", "")
            
            if skip >2:
                self.shifts.append(float(line.split(';')[0]))
                
            skip+=1
    
    def readCounts(self, specdir):
        
        stdout.write("\n \t... Read Corrected Counts \n\n"); stdout.flush()
                        
        for setofmeasuremnts in self.files:
            
            setofcounts = []
            for measurement in setofmeasuremnts:
                
                counts = []
                for line in open(specdir + '/' + measurement,'r'):
                    line = line.replace("\n", "")
                    line = line.replace("\r", "")
                    counts.append(float(line))

                # median filtering of signal, Gaussian overshoots

                #kernel = signal.gaussian(11, 1)
                #smoothed = signal.fftconvolve(counts, kernel, mode='same')
                #setofcounts.append(smoothed)
                
                #setofcounts.append(signal.medfilt(counts, kernel_size=3))
                setofcounts.append(counts)

            # averaging couts
            data = asarray(setofcounts)
            avrg = []
            for k in range(1024):
                avrg.append(median(data[0:10,k]))

            self.counts.append(avrg)
            
        
        # in case of self.filename[19] -> background correction lead to negative values -> shift by min[self.counts[19]]
        # self.plotFunction(self.shifts, self.counts[19], "Raman Shifts 1/cm", "Corrected Counts", "Raman Spectrum")
        self.counts[19][0:] = [x + abs(min(self.counts[19]))/2. for x in self.counts[19][0:]]
        
    
    
    def plotFunction(self, x, y, xlab, ylab, title, vlines = None, filename = None): 
           
        stdout.write("\n \t... Plot Function \n\n"); stdout.flush()
                    
        if (filename==None and vlines==None):
            plt.scatter(x, y, color='b')
            plt.grid(True)
            plt.ylabel(ylab)
            plt.xlabel(xlab)
            plt.title(title)
            #plt.axis(axislimits)            
            ax = plt.gca()
            ax.set_autoscale_on(False)
            plt.show()
        elif (filename):
            plt.ioff()
            fig = plt.figure(figsize=(4.5, 4.5), dpi=600)
            plt.scatter(x, y, color='b', s=5)#plt.plot(x, y, linewidth=2.0, color='b')
            plt.grid(True)
            plt.ylabel(ylab)
            plt.xlabel(xlab)
            plt.title(title)            
            #plt.axis(axislimits)
            ax = plt.gca()
            ax.set_autoscale_on(False)
            plt.savefig(filename, dpi = 600, format='pdf')
            plt.close(fig)
        
        if (vlines):
            plt.scatter(x, y, color='b')
            plt.grid(True)
            for i in range(len(vlines)):
                plt.axvline(x=vlines[i], ymin=0, ymax = 250, linewidth=2, color='k')            
            
            plt.ylabel(ylab)
            plt.xlabel(xlab)
            plt.title(title)
            #plt.axis(axislimits)
            ax = plt.gca()
            ax.set_autoscale_on(False)
            plt.show()
                

    def isolateBands(self):
        
        stdout.write("\n \t... Isolate Bands of Raman Shifts \n\n"); stdout.flush()

        """
        v1PO4   930 - 980, final version: 900 - 1000
        v2PO4   410 - 460, final version: 370 - 500 # note: upper bound for v1PO4 was changed from 460 to 459 to ease detection
        amid1   1620 - 1700, final version: 1620 - 1720
        amid3   1215 - 1300, final version: 1150 - 1350
        CH3     1365 - 1390
        pyd     1660
        lipid   1298
        """
        
        #lines = [930, 980, 410, 459, 1620, 1700, 1215, 1300, 1365, 1390, 1660, 1298]
        #lb = [930, 410, 1620, 1215, 1365]
        #ub = [980, 459, 1700, 1300, 1390]

        lines = [900, 1000, 370, 500, 1620, 1720, 1150, 1350, 1365, 1390, 1660, 1298]
        lb = [900, 370, 1620, 1150, 1365]
        ub = [1000, 500, 1720, 1350, 1390]

        
        #self.plotFunction(self.shifts, self.counts[19], "Raman Shifts 1/cm", "Corrected Counts", "Raman Spectrum", lines)

        # getting indizees to isolate peaks
        for k in lb:
            self.lbindex.append(self.shifts.index(buin.min(self.shifts, key=lambda x:abs(x - k))))
                
        for k in ub:
            self.ubindex.append(self.shifts.index(buin.min(self.shifts, key=lambda x:abs(x - k))))
                
        # identifying pyd and lipid bounds to interpolate at
        # find closest neighbour to given value in list
        # print min(self.shifts, key=lambda x:abs(x-1660))

        # pyd   1660        
        l = self.shifts.index(buin.min(self.shifts, key=lambda x:abs(x-1660)))-1
        u = self.shifts.index(buin.min(self.shifts, key=lambda x:abs(x-1660)))
        self.lbindex.append(l)
        self.ubindex.append(u)
        
        # lipids    1298
        l = self.shifts.index(buin.min(self.shifts, key=lambda x:abs(x-1298)))-1
        u = self.shifts.index(buin.min(self.shifts, key=lambda x:abs(x-1298)))
        self.lbindex.append(l)
        self.ubindex.append(u)
        
        
        
    def peakAnalysis(self):
        
        stdout.write("\n \t... Analyse peaks \n\n"); stdout.flush()        
        
        # analysing v1PO4   930 - 980
        for specimen in self.counts:
            
            #fit, integral = self.fwhm(self.shifts[self.lbindex[0]:self.ubindex[0]], self.counts[18][0][self.lbindex[0]:self.ubindex[0]])
            shifts = self.shifts[self.lbindex[0]:self.ubindex[0]]
            counts = specimen[self.lbindex[0]:self.ubindex[0]]
            fit, integral = self.fwhm(shifts, counts)
            
            self.v1PO4fwhm.append(2*fit[0])
            self.v1PO4int.append(integral[0])
            self.v1PO4peak.append(fit[2])
            

        # analysing v2PO4   410 - 460
        for specimen in self.counts:
            
            #fit, integral = self.fwhm(self.shifts[self.lbindex[0]:self.ubindex[0]], self.counts[18][0][self.lbindex[0]:self.ubindex[0]])
            shifts = self.shifts[self.lbindex[1]:self.ubindex[1]]
            counts = specimen[self.lbindex[1]:self.ubindex[1]]
            fit, integral = self.fwhm(shifts, counts)
                    
            self.v2PO4int.append(integral[0])
        

        # analysing amid1   1620 - 1700
        for specimen in self.counts:

            #fit, integral = self.fwhm(self.shifts[self.lbindex[0]:self.ubindex[0]], self.counts[18][0][self.lbindex[0]:self.ubindex[0]])
            shifts = self.shifts[self.lbindex[2]:self.ubindex[2]]
            counts = specimen[self.lbindex[2]:self.ubindex[2]]
            fit, integral = self.fwhm(shifts, counts)
                
            self.amid1int.append(integral[0])
            self.amid1peak.append(fit[2])
            

        # analysing amid3   1215 - 1300
        for specimen in self.counts:
            
            shifts = self.shifts[self.lbindex[3]:self.ubindex[3]]
            counts = specimen[self.lbindex[3]:self.ubindex[3]]
            fit, integral = self.fwhm(shifts, counts)

            self.amid3int.append(integral[0])
            
            
        # analysing CH3     1365 - 1390
        for specimen in self.counts:

            shifts = self.shifts[self.lbindex[4]:self.ubindex[4]]
            counts = specimen[self.lbindex[4]:self.ubindex[4]]
            fit, integral = self.fwhm(shifts, counts)
        
            self.CH3int.append(integral[0])


        # analysing pyd     1660
        for specimen in self.counts:
                
            # linear interpolation based on two-point-equation
            x0 = self.shifts[self.lbindex[5]]
            x1 = self.shifts[self.ubindex[5]]

            y0 = specimen[self.lbindex[5]]
            y1 = specimen[self.ubindex[5]]        
                
            self.pyd.append(y0 + (y1-y0)/(x1-x0)*(1660-x0))
        
        
        # analysing lipid   1298
        for specimen in self.counts:

            # linear interpolation based on two-point-equation
            x0 = self.shifts[self.lbindex[6]]
            x1 = self.shifts[self.ubindex[6]]

            y0 = specimen[self.lbindex[6]]
            y1 = specimen[self.ubindex[6]]        
                
            self.lipid.append(y0 + (y1-y0)/(x1-x0)*(1298-x0))
        


    # fitting full width at half maximum see http://mesa.ac.nz/mesa-resources/technical-tutorials/python-2/python-workshop-i-fitting-a-single-symmetric-peak/ visited May 11, 2015
    # the old version was done with a lorentzian    
    def fwhmLorentzian(self, x, y):
        
        # initial values [hwhm, peak center, intensity]
        # see http://en.wikipedia.org/wiki/Full_width_at_half_maximum for the hwhm estimation
        p = [2.355*std(x)/2.0, (max(x)+min(x))/2.0, max(y)]
        
        # optimization
        pbest = leastsq(self.residuals, p, args=(y,x),full_output=1)
        best_parameters = pbest[0]
        
        # fit to data
        fit = self.lorentzian(x, best_parameters)
        integral = integrate.quad(lambda x: self.lorentzian(x, best_parameters), min(x), max(x))
        
        # plot overlay
        plt.plot(x,y,'wo')
        plt.plot(x,fit,'r-',lw=2)
        plt.xlabel(r'$\omega$ in cm$^{-1}$', fontsize=18)
        plt.ylabel('Intensity (a.u.)', fontsize=18)
        plt.show()
        
        return best_parameters, integral
        
    def fwhm(self, x, y):
    
        # p = [intensity, peak center, std]
        p = [max(y), (max(x)+min(x))/2.0, 1.0]
        #get Gaussian parameters
        pbest = leastsq(self.residualsGaussian, p, args=(y,x),full_output=1)        
        best_parameters = pbest[0]
        
        # get fit to data and integral
        #fit = self.lorentzian(x, best_parameters)
        #integral = integrate.quad(lambda x: self.lorentzian(x, best_parameters), min(x), max(x))
        fit = self.gaussian(x, best_parameters)        
        integral = integrate.quad(lambda x: self.gaussian(x, best_parameters), min(x), max(x))
        
#        # plot overlay
#        plt.plot(x,y,'wo')
#        plt.plot(x,fit,'r-',lw=2)
#        plt.xlabel(r'$\omega$ in cm$^{-1}$', fontsize=18)
#        plt.ylabel('Intensity (a.u.)', fontsize=18)
#        plt.show()
        
        # recast best parameters to [fwhm, peak center, intensity]        
        para = [2*sqrt(2*log(2))*best_parameters[2], best_parameters[1], best_parameters[0]]
        
        return para, integral
        
        
    def fwhm2(self, x, y):
       
        # initial values [hwhm, peak center, intensity]
        # see http://en.wikipedia.org/wiki/Full_width_at_half_maximum for the hwhm estimation
        p = [2.355*std(x)/2.0, (max(x)+min(x))/2.0, max(y)]
        
        # optimization
        #get Lorentz parameters
        pbest = leastsq(self.residuals, p, args=(y,x),full_output=1)        
        best_parameters = pbest[0]

        # get array from left peak - fwhm to end
        #print best_parameters[1] - best_parameters[0], (x.index(min(x, key=lambda omega:abs(omega - (best_parameters[1] - best_parameters[0]))))), x[(x.index(min(x, key=lambda omega:abs(omega - (best_parameters[1] - best_parameters[0])))))]
        
        # defining x and y new
        xnew = x[(x.index(buin.min(x, key=lambda omega:abs(omega - (best_parameters[1] - best_parameters[0]))))):len(x)]
        ynew = y[(x.index(buin.min(x, key=lambda omega:abs(omega - (best_parameters[1] - best_parameters[0]))))):len(y)]

        # p = [intensity, peak center, std]
        p = [max(y), best_parameters[1], 1.0]
        #get Gaussian parameters
        pbest = leastsq(self.residualsGaussian, p, args=(ynew,xnew),full_output=1)        
        best_parameters = pbest[0]
        
        # get fit to data and integral
        #fit = self.lorentzian(x, best_parameters)
        #integral = integrate.quad(lambda x: self.lorentzian(x, best_parameters), min(x), max(x))
        fit = self.gaussian(x, best_parameters)        
        integral = integrate.quad(lambda x: self.gaussian(x, best_parameters), min(x), max(x))
        
        # plot overlay
#        plt.plot(x,y,'wo')
#        plt.plot(x,fit,'r-',lw=2)
#        plt.xlabel(r'$\omega$ in cm$^{-1}$', fontsize=18)
#        plt.ylabel('Intensity (a.u.)', fontsize=18)
#        plt.show()
        
        # recast best parameters to [fwhm, peak center, intensity]        
        para = [2*sqrt(2*log(2))*best_parameters[2], best_parameters[1], best_parameters[0]]
        
        return para, integral
        
                
    def lorentzian(self, x, p):
        
        if len(p) == 3:
            numerator =  ( p[0]**2 )
            denominator = ( x - (p[1]) )**2 + p[0]**2
            y = p[2]*(numerator/denominator)
        
        elif len(p) == 6:            
            numerator1 =  ( p[0]**2 )
            denominator1 = ( x - (p[1]) )**2 + p[0]**2
            
            numerator2 =  ( p[3]**2 )
            denominator2 = ( x - (p[4]) )**2 + p[3]**2
            
            y = p[2]*(numerator1/denominator1) + p[5]*(numerator2/denominator2)
                        
        return y
        
    def gaussian (self, x, p):
        
        numerator = -1.0*(x - p[1])**2
        denominator = 2.0*p[2]**2
        y = p[0]*exp(numerator/denominator)
        return y
        
    def residuals(self, p, y, x):        
        err = (y - self.lorentzian(x, p))**2
        return err
        
    def residualsGaussian(self, p, y, x):
        err = (y - self.gaussian(x, p))**2
        return err


    # write tabseparated file
    def writeData(self):
            
        stdout.write("\n \t... Write results \n\n"); stdout.flush()
        
        f = open('ramanavrgspectra-mean.dat', 'w')
        f.write('filename\tage\tgender\tv1PO4fwhm\tv1PO4int\tv1PO4peak\tv2PO4int\tamid1int\tamid1peak\tamid3int\tpyd\tCH3int\tlipid\n')

        for k in range(len(self.filename)):

            f.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
            self.filename[k],
            self.age[k],
            self.gender[k],
            str(self.v1PO4fwhm[k]),
            str(self.v1PO4int[k]),
            str(self.v1PO4peak[k]),
            str(self.v2PO4int[k]), 
            str(self.amid1int[k]), 
            str(self.amid1peak[k]), 
            str(self.amid3int[k]), 
            str(self.pyd[k]), 
            str(self.CH3int[k]), 
            str(self.lipid[k])))

        f.close()
        
        
        # Create a new directory if it does not exist 
        isExist = os.path.exists('averagedSpectra')
        if not isExist:
            os.makedirs('averagedSpectra')
            stdout.write("\n \t    New output directory for averaged spectra is created!\n\n"); stdout.flush()
        
        
        for k in range(len(self.filename)):
            
            avrgfilename = 'averagedSpectra/' + self.filename[k] + '_avrg_spectra.dat'            
            
            f = open(avrgfilename, 'w')
            f.write('%s\t%s\n' % ( 'shifts', 'counts' ))
            
            for i in range(len(self.shifts)):
                f.write('%s\t%s\n' % ( str(self.shifts[i]), self.counts[k][i]))
            
            f.close()

if __name__ == '__main__':

    ################################################################
    # initialising classes to extract data
    raman = RamanData()
    
    ################################################################
    stdout.write( '\n S T A R T   %s \n\n' % "EXTRACT RAMAN DATA" )
    stdout.flush()
    startTime = time.perf_counter()
    ctime1    = time.time()

    ################################################################
    # read in put parameters
    dir, spec, meas, shifts = raman.getInputs()

    ################################################################
    # read raman data
    raman.readFiles(dir, spec, meas)
    raman.readShifts(shifts)
    raman.readCounts(dir)
    
    ################################################################
    # compute raman data
    raman.isolateBands()
    raman.peakAnalysis()

    ################################################################
    # write raman data
    raman.writeData()
    
    endTime = time.perf_counter()
    ctime2 = time.time()
    stdout.write( '\n E N D E D  SUCCESSFULLY in  CPU/TOT : %8.1f/%8.1f sec \n\n' % (endTime-startTime, ctime2-ctime1) )
    stdout.flush()
    ################################################################
