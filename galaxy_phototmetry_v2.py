#!/usr/bin/env python3

# Galaxy surface photometry 
# Provides isophotal analysis of the target centered in the frame or specified by certain coordinates. 
# The script attempts a 4 step ellipse fitting and returns different plots for verification and an output file with the information about the script execution steps and physical parameters of the target.

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
#plt.style.use("viki");

import numpy as np
import sys
import matplotlib as mpl   #setup the cmap colours
from copy import deepcopy
import argpar
import argparse
import warnings
import os.path 


from astropy.modeling import models, fitting    # Fitting specific function to the intensities 
from scipy.optimize import curve_fit

from astropy.io import fits
from astropy.wcs import WCS

from photutils.isophote import EllipseGeometry, Ellipse, build_ellipse_model 
from photutils import Background2D, MedianBackground
#from itertools import groupby    #calculating the cutoff index
#from operator import itemgetter, sub   #  "itemgetter" used only for cutoff index

import tkinter.simpledialog   # getting an input parameter while the script is running 
import tkinter as tk

import pickle   # prints results into file. 




######################################## Auxiliary FUNCTIONS ####################################################################################################
############################################################################################################################################################

# Create a 1/error list for calculating the weights
def error_calculate(error_list):   

    minval = np.nanmin(error_list[np.nonzero(error_list)])
    return map(lambda err: 1. / minval if err == 0 or np.isnan(err) else 1. / (err), error_list)
    
# Getting the weighted avg value, based on each points error range
def get_weighted_avg(input_list, error_list): #  weighted average  gives more importance to the central part of the target 

    err = list(error_calculate(error_list))
    indices = np.where(np.logical_not(np.isnan(input_list)))[0]
    return np.average(np.array(input_list)[indices], weights=np.array(err)[indices])        

def get_string(text): #if some header parameters are missing an input window appears for asking the certain value    
        application_window = tk.Tk()
        application_window.withdraw()
        a = tk.simpledialog.askstring("Input", "Ennetr new "+str(text)+"value", parent=application_window, initialvalue="Input")
        application_window.after(60, application_window.destroy)
        return a

def get_png(png_name): #checks if the input png is correct otherwise asks for a new one
        while not os.path.isfile(png_name.replace(".png", "")+'.png') :    
            png_name=get_string("PNG name:")
        return png_name.replace(".png", "")
                    
#Magnitude flux conversion   
def mag_star(zp, flux):
    return zp - 2.5 * np.log10(flux)

#Fitting functions
def func_Sersic_4(x, a, b, c):
    """ https://astronomy.swin.edu.au/cosmos/S/Surface+Brightness+Profiles """ 
    return a * np.exp(b * -(x)**(1/4)) + c
    
def func_Vaucouleurs(x,a,b,c):
    return a * np.exp(-7.67*(x/b)**(1/4)-1)
    

def func_Sersic(x,a,b,c):
    k= models.Sersic1D(a,b,c)
    return k(x)
    
def func_exp(x,a,b,c):
     return a * np.exp(-x/b)


######################################## CLASS BEGINING ####################################################################################################
############################################################################################################################################################
class galaxy_photometry:

    ### Setting up the parameters
    def  __init__(self,fits_name,zp, target_name):
        self.fits_name=fits_name
        self.zp=zp    #zero point, external parameter
        self.target_name= target_name
        
        self.hdul = fits.open(self.fits_name)                    # full open fits file
        self.data = self.hdul[0].data                            # image part of the fits file
        self.header=WCS(self.hdul[0].header)                     #coordinates from the fits file
        self.x0 = float(self.data.shape[0]) / 2.                # galaxy center
        self.y0 = float(self.data.shape[1]) / 2.                # galaxy center
        self.maxsma = self.x0 if self.x0 > self.y0 else self.y0  #self.data.shape[0]  #  the maximum semi major axis value is half of the image size 
        self.file_name=self.fits_name.split('.fits')[0]
        self.png_name=self.file_name                            # if png name is not specified it will search for it with the same name as the fits file

        self.isolist=None  # object returned after the elliptical isophote fit
        self.sma = 3       # parameter  used as  EllipseGeometry sma  ("https://photutils.readthedocs.io/en/stable/api/photutils.isophote.IsophoteList.html#photutils.isophote.IsophoteList")
        self.ellip = 0.2   # parameter used as EllipseGeometry eps
        self.pa = 0        # parameter used as EllipseGeometry pa
        self.astep= 0.2    # parameter used as EllipseGeometry astep     
        self.fix_center=False # parameter used as EllipseGeometry fix_center
        self.fix_pa=False     # parameter used as EllipseGeometry fix_pa
        self.fix_eps=False    # parameter used as EllipseGeometry fix_eps 
        self.failed=False     # verifies if a certain step failed 

        self.check_intens_fluctuate=True    # search for the 3rd inflection point after reaching the background level,  in case it is fails it stops at the background level. Most of the cases was better results obtained with this parameter kept false. 
        
        self.bkg=None                       #background parameters
        self.background_scale = 10.         #  division factor of the image to obtain the box size  (https://photutils.readthedocs.io/en/stable/api/photutils.background.Background2D.html#photutils.background.Background2D)
        self.i = 0  # calculating the repetition for the images- this parameter is used for saving the plots name 
       
    px_scale= None        #parameters from the header
    filter_name=None
    seeing_header=None
    combining_method=None   

    cmap=None  # imshow jet profile parameters
    vmin=50
    
    iraf_result_file=None  # the results can be compared with  results from other sources by inputting a 3 column file sma, mag_isophotal, mag_total

### Ignoring warning during the "try" process

    warnings.filterwarnings("error", category=UserWarning)
    warnings.simplefilter("ignore")

    def get_data_from_header(self): # THELI headers can give us extra information in case of using other data reduction method may require verifying these values
     
    #Pixel scale from header or external input   
        try:  
            if self.hdul[0].header["INST_TEL"] == 'VIRCAM@VISTA':
                        #    print 'VIRCAM@VISTA'
                                self.px_scale = 0.34
            elif self.hdul[0].header["INST_TEL"] == 'LIRIS@WHT':
                        #    print  'LIRIS@WHT'
                                self.px_scale = 0.25
            elif self.hdul[0].header["INST_TEL"] == 'WIRCam-2@CFHT':
                        #    print  'WIRCam-2@CFHT'
                                self.px_scale = 0.307                  
        except:    
            if not (self.px_scale):       self.px_scale=float(get_string("px_scale"))
        print("Pixel scale:"+str(self.px_scale))    
    #filetr name from header or external input   
        try:
            self.filter_name = self.hdul[0].header["FILTER"]
        except:
            if not self.filter_name:      self.filter_name=get_string("filter_name")
    #seeing value from header or external input       
        try:
            self.seeing_header = self.hdul[0].header["FWHM"]   #["FWHM"] 
        except:
            try: 
                self.seeing_header = self.hdul[0].header["FWHM_W"]   
            except:
                if not self.seeing_header:    self.seeing_header=float(get_string("seeing_header"))
        print("Seein header:", self.seeing_header)
    # the  coaddition method used during data reduction from header or external input   
        try:
            self.combining_method = self.hdul[0].header["COMBINET"]
        except:        
            if not self.combining_method: self.combining_method=get_string("combining_method")
        print("Combining method:"+self.combining_method)
            
    def setup_cmap(self):   # for getting a better representation customized jet colourmap was used  
        
        upper = mpl.cm.jet(np.arange(246))
        lower = np.ones((int(246/4),4))
        for i in range(3):
          lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
        cmap = np.vstack(( lower, upper ))
        self.cmap = mpl.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

    def calculate_bkg_level(self):  # calculating the background level
        
        print("Backgroud /scale:", self.background_scale)
        bkg_estimator = MedianBackground()
        self.bkg = Background2D(self.data, (int(float(self.data.shape[0]) / self.background_scale), int(float(self.data.shape[0]) / self.background_scale))) 
        print("Background", np.min(self.bkg.background), np.max(self.bkg.background))
        
    def substract_bkg(self): # set the image background to 0, by subtracting the median background and plots the result. 
   
        fig, (ax1, ax2, ax3) = plt.subplots(figsize=(24, 7),nrows=1, ncols=3)     
       
        #original image
        im1=ax1.imshow(self.data, origin='lower', cmap= self.cmap,vmin=-0.5, vmax=0.5)
        plt.colorbar(im1, ax=ax1)
        ax1.set_title('Reduced image',pad=40)
        ax1.grid(color='white', ls='--')
       
        self.data=self.data-self.bkg.background                    #  DATA-BACKGROUND MODEL
        
        background_hdu=fits.PrimaryHDU(self.data)
        background_hdu.writeto("Background_substacted_"+self.file_name+".fits",overwrite=True)  # write the background sutracted image into a fits file
        
        #background model 
        im2=ax2.imshow(self.bkg.background, origin='lower', cmap= self.cmap)
        plt.colorbar(im2, ax=ax2)
        ax2.set_title('Background',pad=40)
        ax2.grid(color='white', ls='--')
        
        #background subtracted data
        im3=ax3.imshow(self.data, origin='lower', cmap= self.cmap,vmin=-0.5, vmax=0.5)
        plt.colorbar(im3, ax=ax3)
        ax3.set_title('Background substracted',pad=40)
        ax3.grid(color='white', ls='--')
        plt.savefig(self.file_name+str(self.i)+"_background_model.png")
        plt.show(block=False)  
        plt.pause(3)
        plt.close()

    def create_ellipses(self):   #  most important function it generates the isolist file which contains the fitted ellipses and their details.  https://photutils.readthedocs.io/en/stable/isophote.html
        
        geometry = EllipseGeometry(x0=self.x0, y0=self.y0, sma=self.sma, eps=self.ellip, pa=self.pa, fix_center=self.fix_center, fix_pa=self.fix_pa, fix_eps=self.fix_eps, linear_growth=True, astep=self.astep)
        print("grometry ",geometry.sma, geometry.astep,self.x0,self.y0,self.ellip, self.pa)
        ellipse = Ellipse(self.data, geometry)
        params=ellipse.fit_isophote(5) #helping for decision of the new parameters in case of low s/n, it fits only one isophote 
        print("parameters for sma 5",params.x0,params.y0,params.eps,params.pa )
        self.isolist = ellipse.fit_image(maxsma=self.maxsma,minsma=self.sma,step=self.astep) #,maxrit=1) #,integrmode='median' 
        if  (self.isolist.sma.size <1) : 
                print("Failed:", self.isolist.sma, "Check the fits image for new input parameters!")
                exit(0)
        self.i = self.i + 1
     #  data = pickle.load(open("save.p", "rb"))        # Open saved data files and read into a variable
     #  self.isolist = data


#Gets the last valid index and verify it, also creates verification plots. Based on what the other parameters will be calculated.     
    def plot_verify_last_valid_index(self):
            
            self.last_valid_index=self.get_intens_fluctuate()  #  generates the last valid index value 
            
            min_x0_err = max(self.isolist.x0) - min(self.isolist.x0)   # plotting with error bar sets the error bar limit low enough so the structure of the official plot still visible
            min_x0_err = 5.  if min_x0_err < 1 else min_x0_err
            min_y0_err = max(self.isolist.y0) - min(self.isolist.y0)
            min_y0_err = 5.  if min_y0_err < 1 else min_y0_err
           
            fig, (ax1,ax2) = plt.subplots(1,2)  
             
            #plots the weighted average in case for the central coordinates and their fluctuations
            ax1.errorbar(self.isolist.sma, self.isolist.x0, yerr=np.array(list(map(lambda x: x if x < min_x0_err else min_x0_err, self.isolist.x0_err))), fmt='o',markersize=4, ecolor='skyblue',zorder=1)
            ax2.errorbar(self.isolist.sma, self.isolist.y0, yerr=np.array(list(map(lambda x: x if x < min_y0_err else min_y0_err, self.isolist.y0_err))), fmt='o',markersize=4, ecolor='skyblue',zorder=1)
            if self.last_valid_index >0:
              #  self.x0=np.median(self.isolist.x0[:self.last_valid_index])   # by using the median the outer part of the galaxy gets more effect in the calculations
                self.x0 = get_weighted_avg(self.isolist.x0[:self.last_valid_index],self.isolist.x0_err[:self.last_valid_index])  #x center coordinate; central  region more effect on the center of the galaxy calculation than the outer region
                ax1.hlines(y=self.x0,xmin=self.isolist.sma[0], xmax=self.isolist.sma[self.last_valid_index], linestyle='--', color='m', zorder=3)
                ax1.annotate(str(self.x0), (0, self.x0))

                #self.y0=np.median(self.isolist.y0[:self.last_valid_index])
                self.y0 = get_weighted_avg(self.isolist.y0[:self.last_valid_index], self.isolist.y0_err[:self.last_valid_index]) #y center coordinate
                ax2.hlines(y=self.y0,xmin=self.isolist.sma[0], xmax=self.isolist.sma[self.last_valid_index],  linestyle='--', color='m', zorder=4, label="intensity")
                ax2.annotate(str(self.y0), (0, self.y0))
                
            ax1.set_xlabel('Semimajor Axis Length (pix)')
            ax1.set_ylabel('x0')            
            ax2.set_xlabel('Semimajor Axis Length (pix)')
            ax2.set_ylabel('y0')

            plt.savefig(self.file_name+str(self.i)+"_1.png")
            
            #plot isophotal intensity            
            plt.figure()  

            plt.plot(self.isolist.sma*self.px_scale,self.isolist.intens,"o")
            plt.plot(self.isolist.sma[self.last_valid_index]*self.px_scale, self.isolist.intens[self.last_valid_index],"o", color='m' )
            plt.axvline(x=self.isolist.sma[self.last_valid_index]*self.px_scale, color='m', linestyle='--')
    
            plt.xlabel('Semimajor Axis Length (arcsec)')
            plt.ylabel('Isophotal intensity')
            
            plt.savefig(self.file_name+str(self.i)+"_2.png")     
            
            #plot isophotal magnitude          
            plt.figure()
            if (self.iraf_result_file) :   # in case verification file exists 
                X =[]
                Y= []
                Z=[]
                for line in open(self.iraf_result_file, 'r').readlines():
                     values = [float(s) for s in line.split()]
                     if (values[1] != self.zp):
                        X.append(values[0])
                        Y.append(values[1])
                        Z.append(values[2])
                plt.plot(np.array(X)*self.px_scale, np.array(Y),zorder=4)
                
            mag_list=(np.array(mag_star(self.zp, self.isolist.intens/ self.px_scale / self.px_scale)))  #converting intensity to magnitude, including the instrumental magnitude
            mag_err_list=2*1.08/self.isolist.intens*self.isolist.int_err            #converting the intensity errors into magnitude errors using the error propagation formula 
            plt.errorbar(self.isolist.sma*self.px_scale, mag_list, yerr= mag_err_list, fmt='o',   markersize=4, ecolor='skyblue',zorder=3) 
            plt.plot(self.isolist.sma[self.last_valid_index]*self.px_scale,mag_list[self.last_valid_index],"o", color='m',zorder=2 )
            plt.axvline(x=self.isolist.sma[self.last_valid_index]*self.px_scale, color='m', linestyle='--',zorder=1)
            plt.xlabel('Semimajor Axis Length (arcsec)')
            plt.ylabel('Isophotal magnitude')
            extra=(max(mag_list) -min(mag_list) )*0.05   
            plt.ylim([max(mag_list)+extra, min(mag_list)-extra])   # set correct inverse plotting range
            plt.savefig(self.file_name+str(self.i)+"_3.png")     
                
            # plotting  the total magnitude/ellipses in function of SMA
            plt.figure()
            plt.gca().invert_yaxis()
            magnitude_list = mag_star(self.zp,self.isolist.tflux_e / self.px_scale / self.px_scale)   #converting intensity to magnitude, including the instrumental magnitude
            final_mag=magnitude_list[self.last_valid_index]
            t_mag_error=np.mean(1.08/self.isolist.intens[:self.last_valid_index]*self.isolist.int_err[:self.last_valid_index] )  #,np.mean(1.08/self.isolist.intens[self.last_valid_index-5:self.last_valid_index+5]*self.isolist.int_err[self.last_valid_index-5:self.last_valid_index+5] )
            print("final_magnitude",final_mag,"intens_error in intens and mean magnitude error, magnitude error:", self.isolist.int_err[self.last_valid_index], t_mag_error)   # printing  the total magnitude and the  error value   
            plt.plot(self.isolist.sma*self.px_scale, magnitude_list)
            plt.plot(self.isolist.sma[self.last_valid_index:]*self.px_scale , [final_mag] * len(self.isolist.sma[self.last_valid_index:]), '-.')
            plt.xlabel('Semimajor Axis Length (arcsec)')
            plt.ylabel('Total magnitude')
            plt.annotate(str(final_mag), (self.isolist.sma[self.last_valid_index]*self.px_scale,final_mag ))
            plt.savefig(self.file_name+str(self.i)+"_4.png")
            
            plt.close()
      
    def  get_intens_fluctuate(self): #intensity value under the  noise level and we are keeping the third inflection point; generates the last valid index 
            self.calculate_bkg_level()
            bg_noise=np.max(self.bkg.background) #self.bkg.background_rms_median  # getting the background noise level
            print("Noise level", bg_noise, self.bkg.background_median )
            low_intens=list(filter(lambda x: x<bg_noise, self.isolist.intens)) #find elements where the intensity is low, comparative with the noise level
            if (self.check_intens_fluctuate) :   #searching for the third inflection point 
                loc_changes=np.where(np.diff(np.sign(np.diff(low_intens) )))[0] #locations where the intensity ch)nges. # sign: -1 if x < 0, 0 if x==0, 1 if x > 0. # diff: x1-x
                
                new_loc_changes=np.diff(loc_changes)
           
                for k in range(len(new_loc_changes)-2):
                    if all(i < 5 for i in new_loc_changes[k:k+10]) :  # the  max distance between the inflection points are 5 data points
                        final=loc_changes[new_loc_changes[k]+1]+2                          #adding 2 to compensate the np.diff
                        print("Intensity successfully found")    
                        return self.isolist.intens.tolist().index(low_intens[final])      # searching for the index of the last valid point in the big list
                        break
                        
            try:             
                return self.isolist.intens.tolist().index(low_intens[2])   # if the search for the third inflection point fails or skipped the 3rd value under the noise level is used. 
            except:
                print("!!!!!Bigger CROP needed!!!!")                    # in case the noise level cannot be reached one possibly needs a bigger crop or different input parameters
                return len(self.isolist.intens)-1

###########################################################Creating the plots######################################################################################################
    def plot_4_quadrant(self): # Plots the 4 most important parameter and their mean values what will be used later for fixing the parameters
    
       plt.figure(figsize=(8, 8), dpi=1600)
       plt.subplots_adjust(hspace=0.35, wspace=0.35)
       
        #Ellipticity 
       ax=plt.subplot(2, 2, 1)     
 #      self.ellip=np.median(self.isolist.eps[:self.last_valid_index] )     # the outer region’s ellipticity more effect in the final ellipticity than the center region
       self.ellip = get_weighted_avg(self.isolist.eps[:self.last_valid_index], self.isolist.ellip_err[:self.last_valid_index]) #  define  the average ellipticity, later on used as fixed value 
         
       plt.hlines(y=self.ellip,xmin=self.isolist.sma[0], xmax=self.isolist.sma[self.last_valid_index],  linestyle='--', color='m', zorder=3) # average ellipticity value  inside the galaxy 
       plt.annotate(str(self.ellip), (0, self.ellip))
       plt.errorbar(self.isolist.sma, self.isolist.eps, yerr=self.isolist.ellip_err, fmt='o', markersize=4, ecolor='skyblue')
       plt.xlabel('Semimajor Axis Length (pix)')
       plt.ylabel('Ellipticity')
       plt.locator_params(axis='y', nbins=6)
       plt.locator_params(axis='x', nbins=10)
       extra=(max(self.isolist.eps) -min(self.isolist.eps) )*0.05   
       plt.ylim([ min(self.isolist.eps)-extra,max(self.isolist.eps)+extra])
       plt.gca().set_xlim(left=0)
       xleft, xright = ax.get_xlim()
       ybottom, ytop = ax.get_ylim()
       ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)))
       
       #PA
       ax1=plt.subplot(2, 2, 2) 
       min_pa0_err = (max(self.isolist.pa) - min(self.isolist.pa))
       if min_pa0_err < 1:
            min_pa0_err = 3.
       plt.errorbar(self.isolist.sma, self.isolist.pa, yerr= self.isolist.pa_err, fmt='o', markersize=4, ecolor='skyblue')
       self.pa=np.arctan2(get_weighted_avg(np.sin(self.isolist.pa[:self.last_valid_index]),self.isolist.pa_err[:self.last_valid_index]),  get_weighted_avg(np.cos(self.isolist.pa[:self.last_valid_index]),self.isolist.pa_err[:self.last_valid_index]))  #  average PA value  inside the galaxy size : https://en.wikipedia.org/wiki/Mean_of_circular_quantities
       plt.hlines(y=self.pa,xmin=self.isolist.sma[0], xmax=self.isolist.sma[self.last_valid_index],  linestyle='--', color='m', zorder=3)
       plt.annotate(str(self.pa), (0, self.pa))
       plt.xlabel('Semimajor Axis Length (pix)')
       plt.ylabel('PA (deg)')
       plt.locator_params(axis='y', nbins=6)
       plt.locator_params(axis='x', nbins=10)
       extra=(max(self.isolist.pa) -min(self.isolist.pa) )*0.05
       plt.ylim([max(self.isolist.pa)+extra, min(self.isolist.pa)-extra])
       plt.gca().set_xlim(left=0)
       xleft, xright = ax1.get_xlim()
       ybottom, ytop = ax1.get_ylim()
       ax1.set_aspect(abs((xright-xleft)/(ybottom-ytop)))

       #x center coord
       ax2=plt.subplot(2, 2, 3) 
       min_x0_err = max(self.isolist.x0) - min(self.isolist.x0)
       if min_x0_err < 1:
           min_x0_err = 5.
       plt.errorbar(self.isolist.sma, self.isolist.x0, yerr=self.isolist.x0_err, fmt='o',  markersize=4, ecolor='skyblue') 
       plt.hlines(y=self.x0,xmin=self.isolist.sma[0], xmax=self.isolist.sma[self.last_valid_index],  linestyle='--', color='m', zorder=3) # the average x was calculated in step  plot_verify_last_valid_index
       plt.annotate(str(self.x0), (0, self.x0))
       plt.xlabel('Semimajor Axis Length (pix)')
       plt.ylabel('x0')
       plt.locator_params(axis='y', nbins=6)
       plt.locator_params(axis='x', nbins=10)
       extra=(max(self.isolist.x0) -min(self.isolist.x0) )*0.05
       plt.ylim([max(self.isolist.x0)+extra, min(self.isolist.x0)-extra])
       plt.gca().set_xlim(left=0)
       xleft, xright = ax2.get_xlim()
       ybottom, ytop = ax2.get_ylim()
       ax2.set_aspect(abs((xright-xleft)/(ybottom-ytop)))

       #y center coord
       ax3=plt.subplot(2, 2, 4) 
       min_y0_err = max(self.isolist.y0) - min(self.isolist.y0)
       if min_y0_err < 1:
            min_y0_err = 5.
       plt.errorbar(self.isolist.sma, self.isolist.y0, yerr=self.isolist.y0_err, fmt='o',   markersize=4, ecolor='skyblue')  
       plt.hlines(y=self.y0,xmin=self.isolist.sma[0], xmax=self.isolist.sma[self.last_valid_index], linestyle='--', color='m', zorder=3 ) # the average y was calculated in step  plot_verify_last_valid_index
       plt.annotate(str(self.y0), (0, self.y0))
       plt.xlabel('Semimajor Axis Length (pix)')
       plt.ylabel('y0')
       plt.locator_params(axis='y', nbins=6)
       plt.locator_params(axis='x', nbins=10)
       extra=(max(self.isolist.y0) -min(self.isolist.y0) )*0.05
       plt.ylim([max(self.isolist.y0)+extra, min(self.isolist.y0)-extra])
       plt.gca().set_xlim(left=0)
       xleft, xright = ax3.get_xlim()
       ybottom, ytop = ax3.get_ylim()
       ax3.set_aspect(abs((xright-xleft)/(ybottom-ytop)))
       
       
       plt.savefig(self.file_name + "result"+str(self.i) + ".png")
       
    def plot_mag_iso_functions(self):  #Test fits with basic functions: exp, sersic, vaucouleurs  
        plt.figure()
        plt.gca().invert_yaxis()
        plt.xlabel('Semimajor Axis Length (arcsec)')
        plt.ylabel('Isophotal intensity')

        sma=np.array(self.isolist.sma) * self.px_scale # converting pixel to arcsec 
        plt.plot(sma,self.isolist.intens ,'.',label="isophot_mag")
        length=self.last_valid_index+int(3./self.px_scale/self.astep) # it attempts the fit till the galaxy size + 3"
        try:  # attempts exponential function fitting 
            popt, pcov = curve_fit(func_exp,sma[1:length], np.nan_to_num(self.isolist.intens[1:length]),p0 = [11,1,10])
            print("exp fit parameters",popt)   
            plt.plot(sma[1:length],func_exp(sma[1:length],*popt),label="exp") 
        except:
            print("No exp fit")
        
        try: # attempts exponential Sersic  function fitting 
            popt, pcov = curve_fit(func_Sersic,sma[1:length], np.nan_to_num(self.isolist.intens[1:length])) #bounds=([np.NINF,np.NINF,1.],[np.inf,np.inf,20.])
            print("sersic fit parameters", popt)
            plt.plot(sma[1:length],func_Sersic(sma[1:length],*popt),label="Sersic")
        except:
            print("No  Sersic fit")
        
        try: # attempts exponential Sersic 1/4 function fitting 
            popt, pcov = curve_fit(func_Vaucouleurs,sma[1:length], np.nan_to_num(self.isolist.intens[1:length]))
            print("vaucouleurs fit parameters",popt)
            plt.plot(sma[1:length],func_Vaucouleurs(sma[1:length],*popt),label="Vaucou")
        except:
            print("No Vaucouleurs fit")
        
        plt.legend(framealpha=1, frameon=True)     
        plt.savefig(self.file_name+str(self.i)+"fitted"+".png", dpi=400)
       
            
    def plot_model_image(self): #plots the reziduals  in an area 2x size of galaxy     
        last_ellips=2*self.last_valid_index*self.astep-1 if 2*self.last_valid_index*self.astep <float(self.data.shape[0]) / 2. else float(self.data.shape[0]) / 2.-1    # verify if the galaxy size is outside of the borders of the image 
        model_image = build_ellipse_model(self.data.shape, self.isolist) # reconstruct the model based on the calculated intensities inside consecutive ellipses.  
        residual = self.data - model_image # original image minus model 
        
        residual_hdu=fits.PrimaryHDU(residual)
        residual_hdu.writeto("Residual_"+self.file_name+".fits",overwrite=True)  # writes out the residuals into fits 

        fig, (ax1, ax2, ax3) = plt.subplots(figsize=(24, 7),nrows=1, ncols=3)        
        fig.subplots_adjust(left=0.1, right=0.9, bottom=0.02, top=0.9,wspace=0.6)
        
        #original image
        ax1.imshow(self.data[int(self.y0-last_ellips):int(self.y0+last_ellips),int(self.x0-last_ellips):int(self.x0+last_ellips)], origin='lower', cmap= self.cmap,vmin=np.percentile(self.data, self.vmin),  aspect='equal')
        ax1.set_title('Reduced image',pad=40)
        ax1.ticklabel_format(style='plain')
        smas = np.linspace(int(self.sma), int(last_ellips), 6)
        color=np.array([1.,1.,1.])
        for sma in smas:
            iso = self.isolist.get_closest(sma)
            x, y, = iso.sampled_coordinates()
            ax1.plot(x-int(self.x0-last_ellips), y-int(self.y0-last_ellips), color=color)  
            color=color-[0.06,0.06,0]
        ax1.set_xlabel('x0')
        ax1.set_ylabel('y0')
        plt.grid(color='white', ls='dotted')

        #model
        ax2.imshow(model_image[int(self.y0-last_ellips):int(self.y0+last_ellips),int(self.x0-last_ellips):int(self.x0+last_ellips)], origin='lower', cmap= self.cmap,vmin=np.percentile(self.data, self.vmin),  aspect='equal')
        ax2.set_title('Ellipse Model',pad=40)
        ax2.ticklabel_format(style='plain')    
        ax2.set_xlabel('x0')
        ax2.set_ylabel('y0')
        plt.grid(color='white', ls='dotted')

        #residuals
        ax3.imshow(residual[int(self.y0-last_ellips):int(self.y0+last_ellips),int(self.x0-last_ellips):int(self.x0+last_ellips)], origin='lower', cmap= self.cmap,vmin=np.percentile(self.data, self.vmin),  aspect='equal')
        ax3.set_xlabel('x0')
        ax3.set_ylabel('y0')
        plt.grid(color='white', ls='dotted')
        ax3.set_title('Residual',pad=40)
        
        plt.savefig(self.file_name+str(self.i)+"_model_elipses.png")
        
        

    def plot_ellipses(self):    # MAIN PLOTIING SECTION

        last_ellips_index=self.last_valid_index+3./self.px_scale/self.astep if self.last_valid_index+3./self.px_scale/self.astep < len(self.isolist.sma) -1 else len(self.isolist.sma) -1 #setting the border of the image 3 arcsec  bigger than the last valid index if that value is inside the image bourders  (  float(self.data.shape[0]/self.astep) / 2.+1; the +1 is  last index of the list +1 to obtain all the elemnets  with the [:x] task) 
        last_ellips= int(self.isolist.sma[int(last_ellips_index)])   # sma of the last ellips
        new_intensity = self.isolist.intens[:int(last_ellips_index)]  # list of intensities  until the galaxy size+3" 
        iso_mag=(np.array(mag_star(self.zp, new_intensity / self.px_scale / self.px_scale)))     #converting intensity to magnitude, including the instrumental magnitude
        
        #changing the ra dec  header values to the  values at the center of the frame 
        new_header= deepcopy(self.header)  
        new_coord= self.header.pixel_to_world(int(self.data.shape[0]/2), int(self.data.shape[1]/2))  
        new_header.wcs.crpix = [last_ellips,last_ellips]  
        new_header.wcs.crval=[new_coord.ra.deg,new_coord.dec.deg]
     
        plt.figure(figsize=(30, 7 ))
        
        plt.rc('font', size=14.5) 
        plt.subplots_adjust(hspace=0.2, wspace=0.35)
 
        #plot original image        
        plt.subplot(141)
        img=get_png(self.png_name)
        im2 = mpimg.imread(img+'.png')
        plt.imshow(im2,aspect='equal')
        plt.xlabel(self.target_name,horizontalalignment='left', x=0.0, )
        plt.xticks([])
        plt.yticks([])
        
        #isophotal magnitudes
        ax=plt.subplot(142)
        plt.gca().invert_yaxis()
        mag_err_list=2*1.08/self.isolist.intens[:int(last_ellips_index)]*self.isolist.int_err[:int(last_ellips_index)]  #2x magnitude errors based on  error propagation formula 
        plt.errorbar(self.isolist.sma[:int(last_ellips_index)]*self.px_scale, iso_mag, yerr=mag_err_list, fmt='o',   markersize=4, ecolor='skyblue') 
        plt.plot(self.isolist.sma[self.last_valid_index]*self.px_scale,       iso_mag[self.last_valid_index],"o", color='m' )
        plt.axvline(x=self.isolist.sma[self.last_valid_index]*self.px_scale, color='m', linestyle='--', label = "$sma_{TK}$")  # marking for the  calculated size of the galaxy 
        plt.xlabel('Semimajor Axis (arcsec)')
        plt.ylabel('Isophotal magnitude (Ks mag/sq.arcsec)')
        extra=(max(iso_mag) -min(iso_mag) )*0.05 
        plt.ylim([max(iso_mag)+extra, min(iso_mag)-extra])
        plt.gca().set_xlim(left=0)
        xleft, xright = ax.get_xlim()
        ybottom, ytop = ax.get_ylim()
        ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)))
        ax.legend()

        #    input image with the overplotted ellipses 
        ax1= plt.subplot(143, projection=new_header)
        ax1.coords['ra'].set_major_formatter('hh:mm:ss.ss')
        ax1.coords['ra'].set_ticks(number=4)                      
        plt.imshow(self.data[int(self.y0-last_ellips):int(self.y0+last_ellips),int(self.x0-last_ellips):int(self.x0+last_ellips)], origin='lower', cmap= self.cmap,vmin=np.percentile(self.data, self.vmin),  aspect='equal')   
        plt.xlabel('RA')
        plt.ylabel('Dec')      
        plt.grid(color='gray', ls='dotted')
        # plotting 10 ellipse between the first and last calculated isophots inside the galaxy size
        smas = np.linspace(int(self.sma), last_ellips, 10 )  
        color=np.array([1.,1.,1.])
        for sma in smas[:-1]:
            iso = self.isolist.get_closest(sma)    # finding the closest  ellips to the input sma 
            x, y, = iso.sampled_coordinates()
            plt.plot(x-int(self.x0-last_ellips), y-int(self.y0-last_ellips),  color=color)  # plot the ellipse centered to the galaxy x0, y0 center coordinates 
            color=color-[0.06,0.06,0]
        last_valid_ellips=self.isolist.get_closest((self.isolist.sma[self.last_valid_index])) # mark the calculated galaxy size 
        x, y, = last_valid_ellips.sampled_coordinates()
        plt.plot(x-int(self.x0-last_ellips), y-int(self.y0-last_ellips),linestyle='--', color='m')
              
        plt.arrow(2, 2, 10./self.px_scale,0, head_width=0, head_length=0,fc='white', ec='black', width=0.003)  # plot a 10"  size bar for scaling 
        plt.text(2, 5, '10"',  color='black')    

        # Input image - modelled ellipses 
        ax2=plt.subplot(144, projection=new_header)
        ax2.coords['ra'].set_major_formatter('hh:mm:ss.ss')
        ax2.coords['ra'].set_ticks(number=4) 
        model_image = build_ellipse_model(self.data.shape, self.isolist)
        residual = self.data - model_image
        plt.imshow(residual[int(self.y0-last_ellips):int(self.y0+last_ellips),int(self.x0-last_ellips):int(self.x0+last_ellips)], origin='lower', cmap= self.cmap,vmin=np.percentile(self.data, self.vmin),  aspect='equal')
        plt.xlabel('RA')
        plt.ylabel('Dec') 
        print(self.target_name.replace(" ", "")) 
        plt.savefig("Finals_"+self.target_name.replace(" ", "")+str(self.i)+"_elipses.png")
        
########################################################### MAIN ######################################################################################################
if __name__ == '__main__':

# Getting the parameters

    parser = argparse.ArgumentParser(description='Fits ellipses to an  image and returns the central coordinates, PA, Ellipticity. A table with the photometry results.')
    parser.add_argument('-fits_name', required=True, help='Name of the fits image.')
    parser.add_argument('-target_name', required=True, help='Name of the Galaxy')
    parser.add_argument('-png_name',help='Name of the png image. Def. same as the fits image name')
    parser.add_argument("-zp", required=True, type=float, help='The zero point of the image')
    parser.add_argument("-iraf_file", help='In case one want to compare the python and iraf or other results. File format SMA E_MAG MAG.')
    parser.add_argument("-x0", type=float, help='Possible center coordinates. Def. center of the image ')
    parser.add_argument("-y0", type=float, help='Possible center coordinates. Def. center of the image. ')
    parser.add_argument("-fix_center",type=bool, help='Force the center to be fixed.')
    parser.add_argument("-sma", type=int, help='The semimajor axis of the ellipse in pixels. Def. 2 ')
    parser.add_argument("-ellip", type=float, help='Possible ellipticity. Def. 0 ')
    parser.add_argument("-fix_ellip",type=bool, help='Force the ellipticity to be fixed.')
    parser.add_argument("-pa", type=float, help='Possible position angle. Def. 0 ')
    parser.add_argument("-fix_pa",type=bool, help='Force the PA to be fixed.')
    parser.add_argument("-maxsma", type=float, help='Maximum SMA value.  Def. equal to x0  ')
    parser.add_argument("-astep", type=float, help='The step value for growing/shrinking the semimajor axis. Def. 0.2')
    parser.add_argument("-show", type=bool, help='Show the plots. (boolean)  ')
    parser.add_argument("-bkg_scale", type=float, help='The size of the boxes used to calculate the background value. This value is connection with the image size. Def. 10  ')
    parser.add_argument("-chk_infl", action='store_false', help='Skip checking for the third inflection point under the noise level')
    parser.add_argument("-no_fix_ellip", action='store_false', help='Skip trying to fix Ellipticity')
    parser.add_argument("-no_fix_pa", action='store_false', help='Skip trying to fix Position Angle')
    
    args = parser.parse_args()
    fits_name = args.fits_name
    iraf_resut_file = args.iraf_file
    x0 = (args.x0)
    y0 = (args.y0)
    sma = (args.sma)
    ellip = (args.ellip)
    pa = (args.pa)
    maxsma = (args.maxsma)
    astep = (args.astep)
    zp = (args.zp)
    show_image = args.show
    background_scale = args.bkg_scale
    fix_center=args.fix_center
    fix_eps=args.fix_ellip
    fix_pa=args.fix_pa
    png_name=args.png_name
    target_name=args.target_name
    chk_infl=args.chk_infl
    no_fix_ellip=args.no_fix_ellip
    no_fix_pa=args.no_fix_pa
 
    g_p=galaxy_photometry(fits_name, zp,target_name)
    g_p.setup_cmap()
    
    sys.stdout=open(g_p.file_name+".out","w") # save the print commands output into a file 
    
    # Setting up the parameters
    if (x0):
        g_p.x0=x0
    if (y0):
         g_p.y0=y0
    if  (sma):
        g_p.sma=sma
    if (ellip):
         g_p.ellip=ellip
    if  (pa):
         g_p.pa=pa
    if  (maxsma):
         g_p.maxsma=maxsma
    if (astep):
         g_p.astep=astep
    if  (background_scale):
         g_p.background_scale=background_scale
         
    if (fix_center):
            g_p.fix_center=fix_center
    if (fix_eps):
            g_p.fix_eps=fix_eps
    if (fix_pa):
            g_p.fix_pa=fix_pa   
        
    if (iraf_resut_file):
         g_p.iraf_result_file=iraf_resut_file

    if (png_name):
        g_p.png_name=png_name
    if (chk_infl is not None ):
        g_p.check_intens_fluctuate=chk_infl
          
    g_p.get_data_from_header() # read parameters from headers 
    
    
    # isophotal analysis with variable input parameters 
    g_p.calculate_bkg_level()
    g_p.substract_bkg()
    g_p.create_ellipses()
    
    # saving initial parameters 
    x0_backup=g_p.x0
    y0_backup=g_p.y0
    ellip_backup=g_p.ellip
    pa_backup=g_p.pa

    #initial guess for the size of the galaxy 
    g_p.plot_verify_last_valid_index()
        
    # plot initial results 
    g_p.plot_4_quadrant()
    g_p.plot_model_image()
    
    print("Min/ Max x0, ellip, pa", min(g_p.isolist.x0),max(g_p.isolist.x0), min(g_p.isolist.eps),max(g_p.isolist.eps),  min(g_p.isolist.pa),max(g_p.isolist.pa), )  
    
    # attending to fit the physical parameters of the galaxy, if fails resets the saved parameters otherwise it updates it 
    if not (g_p.fix_center):
        try:
            x0_backup,y0_backup,ellip_backup,pa_backup =g_p.x0,g_p.y0,g_p.ellip,g_p.pa # saving initial parameters
            pickle.dump(g_p.isolist, open(g_p.file_name+"1.p", "wb"))  # saving isolist values 
            g_p.fix_center=True
            g_p.create_ellipses() # attends recalculation of the  isophotal intensities 
            g_p.plot_4_quadrant()
        except:
            g_p.fix_center=False
            g_p.x0,g_p.y0,g_p.ellip,g_p.pa=x0_backup,y0_backup,ellip_backup,pa_backup # resets initial parameters 
            g_p.isolist = pickle.load(open(g_p.file_name+"1.p", "rb")) # resets isolist values 
            print("Fix center failed")
            g_p.failed=True 
    if no_fix_ellip:       
        if not (g_p.fix_eps):
            try:
                x0_backup,y0_backup,ellip_backup,pa_backup =g_p.x0,g_p.y0,g_p.ellip,g_p.pa # saving initial parameters
                pickle.dump(g_p.isolist, open(g_p.file_name+"2.p", "wb")) # saving isolist values 
                g_p.fix_eps=True
                g_p.create_ellipses() # attends recalculation of the  isophotal intensities 
                g_p.plot_4_quadrant()    
            except:
                g_p.fix_eps=False
                g_p.x0,g_p.y0,g_p.ellip,g_p.pa=x0_backup,y0_backup,ellip_backup,pa_backup # resets initial parameters 
                g_p.isolist = pickle.load(open(g_p.file_name+"2.p", "rb"))  # resets isolist values 
                print("Fix ellipticity failed")
                g_p.failed=True
      
    if no_fix_pa:     
        if not (g_p.fix_pa):
            try:    
                x0_backup,y0_backup,ellip_backup,pa_backup =g_p.x0,g_p.y0,g_p.ellip,g_p.pa # saving initial parameters
                pickle.dump(g_p.isolist, open(g_p.file_name+"3.p", "wb")) # saving isolist values 
                g_p.fix_pa=True
                g_p.create_ellipses() # attends recalculation of the  isophotal intensities 
                g_p.plot_4_quadrant()   
            except:
                g_p.fix_pa=False
                g_p.x0,g_p.y0,g_p.ellip,g_p.pa=x0_backup,y0_backup,ellip_backup,pa_backup # resets initial parameters 
                g_p.isolist = pickle.load(open(g_p.file_name+"3.p", "rb")) # resets isolist values 
                print("Fix PA failed")
                g_p.failed=True  
   
    g_p.plot_verify_last_valid_index()  # recalculates galaxy size 
    # cerate final plots 
    g_p.plot_mag_iso_functions()        
    g_p.plot_ellipses()
    g_p.plot_model_image()
    
    #print out final list for further analysis
    pickle.dump(g_p.isolist, open(g_p.file_name+"_Final.p", "wb"))  
    
    # print out physical parameters and important information 
    print("Min/ Max x0, ellip, pa", min(g_p.isolist.x0),max(g_p.isolist.x0), min(g_p.isolist.eps),max(g_p.isolist.eps),  min(g_p.isolist.pa),max(g_p.isolist.pa), )      
    print("Final geometry:",g_p.x0,g_p.y0,g_p.ellip,g_p.pa,g_p.i) 
    print("x0", g_p.x0)
    print("y0",g_p.y0 ) 
    coord=g_p.header.pixel_to_world(g_p.x0, g_p.y0)
    print("coords degree ",coord.to_string('hmsdms') )
    print("pa", g_p.pa)
#    print("pa_error",g_p.isolist.pa_err[min(range(len(g_p.isolist.pa)), key=lambda i: abs(g_p.isolist.pa[i]-g_p.pa))], g_p.isolist.ellip_err[g_p.last_valid_index])
#    print("ellip_error",g_p.isolist.ellip_err[min(range(len(g_p.isolist.eps)), key=lambda i: abs(g_p.isolist.eps[i]-g_p.ellip))], g_p.isolist.pa_err[g_p.last_valid_index])
    print("ellip", g_p.ellip)
    print("px scale", g_p.px_scale)
    print("seeing_header", g_p.seeing_header)
    print("sma", g_p.isolist.sma[g_p.last_valid_index])
    print("sma in arcsec", (g_p.isolist.sma[g_p.last_valid_index]*g_p.px_scale))
    print("file name", g_p.file_name )
    print("Error ellip  median", np.mean(g_p.isolist.ellip_err[:g_p.last_valid_index])) # mean ellip error value inside the galaxy
    print("Error PA median", np.mean(g_p.isolist.pa_err[:g_p.last_valid_index])) # mean PA error value inside the galaxy
    print("Zero Point", g_p.zp)
    print("isophotal mag",np.nanmean(mag_star(zp, g_p.isolist.intens[g_p.last_valid_index-5:g_p.last_valid_index+5]/g_p.px_scale/g_p.px_scale)))
    #Calculating SMA error     
    intens_error=np.mean(g_p.isolist.int_err[:g_p.last_valid_index])
    intens_low=g_p.isolist.intens[g_p.last_valid_index]-intens_error
    intens_high=g_p.isolist.intens[g_p.last_valid_index]+intens_error
    sma_low_index=[n for n,i in enumerate(g_p.isolist.intens[:g_p.last_valid_index][::-1]) if i>=intens_high][0]
    sma_low=g_p.isolist.sma[:g_p.last_valid_index]
    sma_low_v=sma_low[-1*sma_low_index]
    try: 
        sma_high_index=[ n for n,i in enumerate(g_p.isolist.intens[g_p.last_valid_index:]) if i<=intens_low][0]
        sma_high=g_p.isolist.sma[g_p.last_valid_index:]
        sma_high_v=sma_high[sma_high_index]
    except: 
         sma_high_v=-1
    print("SMA error unused", sma_low_v, sma_high_v, (sma_high_v-sma_low_v)/2)
    print("sma error", g_p.isolist.sma[g_p.last_valid_index]-sma_low_v)
    # calculating 2MASS comparison parameters    
    mag=(np.array(mag_star(g_p.zp, g_p.isolist.intens / g_p.px_scale / g_p.px_scale))) 
    index=[ n for n,i in enumerate(mag) if i>=20 ][0]
    print("2MASS Comparison MAG", mag_star(g_p.zp,g_p.isolist.tflux_e[index] / g_p.px_scale / g_p.px_scale))
    print("2MASS Comparison SMA", g_p.isolist.sma[index])
    print("2MASS Comparison SMA in arcsec", g_p.isolist.sma[index]*g_p.px_scale)
    
    sys.stdout.close() 

    

