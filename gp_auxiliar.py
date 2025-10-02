#!/usr/bin/env python3
# Reads the output of the GPv2 script and organizes for furthure usage. Creates also the input file used for fitting teh surfacee brightness profiles. 
import os
import glob 
import re
import argpar
import argparse
import numpy as np
import csv
import pickle 
import sys

#magnitude flux conversion   
def mag_star(zp, flux):
    return zp - 2.5 * np.log10(flux)
    
class galaxy_photometry_aux:
    def  __init__(self,latest_p="",latest_out=""):   #setup requested parameters 
        self.latest_p=latest_p
        self.latest_out=latest_out
        self.isolist=None
        self.last_valid_index=None
        self.px_scale=None
        self.zp=None
        self.astep=None
        self.name=""
        self.seeing=None
        self.bg_diff=0.
        
        
    def find_read_file(self):     # search for last modified  .p and .out files generated in the folder 
        if self.latest_p=="":
            list_of_files = glob.glob(os.getcwd()+'/*.p') 
            self.latest_p = max(list_of_files, key=os.path.getctime)
            print(self.latest_p)
        
        if self.latest_out=="":
            list_of_files = glob.glob(os.getcwd()+'/*.out') 
            self.latest_out = max(list_of_files, key=os.path.getctime)
            print(self.latest_out)
    def read_obj(self):
           self.isolist  = pickle.load(open(self.latest_p, "rb"))        # Open saved data files and read into a variable

    def read_params(self,fName=""):
	    fName=self.latest_out
	    try:
		     open(fName,'r')
	    except:
		    return 0		
	    name_param=0
	    for line in open(fName,'r'):
	        if line.find("coords degree")>-1: 
	            coords=line.split("coords degree")[1].replace('\n','').replace('m',':').replace('d',':').replace('h',':')
	            ra=coords.split()[0]
	            dec=coords.split()[1]
	        if line.find("seeing_header")>-1: 
	            self.seeing=line.split("seeing_header")[1]
	        if line.find("Pixel scale")>-1: 
	            self.px_scale=np.float(line.split("Pixel scale:")[1])   
	        if line.find("Zero Point")>-1: 
	            self.zp=np.float(line.split("Zero Point ")[1])
	        if line.find("sma")>-1: 
	            try:
	                a_T=float(line.split("sma")[1])
	                self.last_valid_index=np.where(abs(self.isolist.sma-a_T) <0.0001)[0][0]
	               # print(a_T, self.isolist.sma[self.last_valid_index])
	               # print("Last valid index",self.last_valid_index)
	            except: 
	                print(line)
	        if line.find("sma error ")>-1: 
	            sigma_a=line.split("sma error")[1]
	        if line.find("2MASS Comarison SMA")>-1: 
	            a_20=line.split("2MASS Comarison SMA")[1]
	        if line.find("final_magnitude")>-1: 
	            m_T=re.findall(r"[-+]?\d*\.\d+|\d+", line )[0]
	            sigma_m=re.findall(r"[-+]?\d*\.\d+|\d+", line )[2]
	        if line.find("geometry")>-1: 
	            astep_1=re.findall(r"[-+]?\d*\.\d+|\d+", line )[0]
	            self.astep=np.float(re.findall(r"[-+]?\d*\.\d+|\d+", line )[1])
	            
	        if line.find("2MASS Comarison MAG ")>-1: 
	            m_20=line.split("2MASS Comarison MAG ")[1]
	        if line.find("isophotal mag")>-1: 
	            m_sky=line.split("isophotal mag")[1]  
	        if line.find("fit parameters")>-1:
	                name_param=1
	        if name_param ==1 and  not line.find("fit parameters")>-1: 
	              self.name=line
	              name_param=0
	        if line.find("Background")>-1:
	              txt, bg_min, bg_max= line.split()
	              bg_diff=float(bg_max)-float(bg_min)
	              if (self.bg_diff<bg_diff):
	                    self.bg_diff=bg_diff
	                        
    def m_22(self): # total magnitude if we define the size of the galaxy at the isophotal intensity 22mag. 
        mag=(np.array(mag_star(self.zp, self.isolist.intens / self.px_scale / self.px_scale))) 
        index=[ n for n,i in enumerate(mag) if i>=22 ][0]
        print("OV Comarison MAG", mag_star(self.zp,self.isolist.tflux_e[index] / self.px_scale / self.px_scale))
        print("OV Comarison SMA", self.isolist.sma[index]*self.px_scale)	 
        
#    def m_20(self):     #Comparing results 
#        mag=(np.array(mag_star(self.zp, self.isolist.intens / self.px_scale / self.px_scale)))    
#        index=[ n for n,i in enumerate(mag) if i>=20 ][0]
#        print("2MASS Comarison MAG", mag_star(self.zp,self.isolist.tflux_e[index] / self.px_scale / self.px_scale))
#        print("2MASS Comarison SMA", self.isolist.sma[index])
#        print("2MASS Comarison SMA", self.isolist.sma[index]*self.px_scale)	
#        print("sma", self.isolist.eps[self.last_valid_index])

#---------------------------------------------------------------generating .csv -------------------------------------------------------------------------------------------------	
                 	               	                
    def gen_csv(self):  
        #Parameters used during fitting to obtain a better fitting results:           
        total_flux=self.isolist.tflux_e[self.last_valid_index]
        half_flux=self.isolist.tflux_e[self.last_valid_index]/2
        t_mag=mag_star(self.zp,total_flux/self.px_scale/self.px_scale)
        h_mag=mag_star(self.zp,half_flux/self.px_scale/self.px_scale)
        try:
             eff_r_index=(abs(self.isolist.tflux_e-half_flux)).argmin()
        except:  
            eff_r_index=int(self.last_valid_index/2.)
            print("eff_r_index failed")
        h_sma=self.isolist.sma[eff_r_index]
        print("Total magnitude:",t_mag )
        print("Half magnitude:", h_mag)
        print(np.where(abs(self.isolist.tflux_e-half_flux) <50) ) 
        print("SMA:", h_sma )
    
    
    #Sampling data  to obtain stepsize coparable with seeing
        size = int(float(self.seeing) / self.px_scale / self.astep) + 1  # convert arcsecond to index 
        sma_px = np.array([sum(group) / size for group in zip(*(self.isolist.sma[x::size] for x in range(size)))]) #creating avg sma 
        intens_avg_px = [sum(group) / size for group in zip(*(self.isolist.intens[x::size] for x in range(size)))]  # creating avg intens 
        int_err_avg_px = [sum(group) / size for group in zip(*(self.isolist.int_err[x::size] for x in range(size)))]  # creating avg intens_error 
        int_err_modified=(np.array(int_err_avg_px)+self.bg_diff)/intens_avg_px

    # Setting the  cutoff of the output 
        last_ellips_index=self.last_valid_index #+int(self.last_valid_index*10/100+1) if self.last_valid_index+int(self.last_valid_index*10/100+1) <len(self.isolist.sma)-1 else len(self.isolist.sma)-1
        
    #CREATING INTENS CSV
        Details = ['Semimajor Axis (arcsec)', 'Isophotal intensity',  'Modified intensity error', self.name+"/"+str(self.isolist.tflux_e[eff_r_index])+"/"+str( self.isolist.sma[self.last_valid_index])+"/"+str(h_sma)+"/"+str(mag_star(self.zp, 1. / self.px_scale / self.px_scale))]  
        rows=np.column_stack((sma_px[:int(last_ellips_index)]*self.px_scale,intens_avg_px[:int(last_ellips_index)],int_err_modified[:int(last_ellips_index)], int_err_avg_px[:int(last_ellips_index)]))
        with open('intens_fit.csv', 'w') as f: 
            write = csv.writer(f) 
            write.writerow(Details) 
            write.writerows(rows) 
            

    #######Converting intensity to magnitude
        mag=(np.array(mag_star(self.zp, self.isolist.intens / self.px_scale / self.px_scale))) 
        mag_err_list=1.08/self.isolist.intens[:int(last_ellips_index)]*self.isolist.int_err[:int(last_ellips_index)]
                
        mag_avg_px = [sum(group) / size for group in zip(*(mag[x::size] for x in range(size)))]  # creating avg intens 
        mag_err_avg_px = [sum(group) / size for group in zip(*(mag_err_list[x::size] for x in range(size)))]  # creating avg intens_error    
        
    #CREATING MAGNITUDE CSV
        Details = ['Semimajor Axis (arcsec)', 'Isophotal magnitude (Ks mag/sq.arcsec)', self.name+"/"+str(t_mag)+"/"+str(self.isolist.tflux_e[eff_r_index])+"/"+str( self.isolist.sma[self.last_valid_index])+"/"+str(self.isolist.intens[0])]  
        rows=np.column_stack((sma_px[:int(last_ellips_index)]*self.px_scale,mag_avg_px[:int(last_ellips_index)],mag_err_avg_px[:int(last_ellips_index)]))
        with open('mag_fit.csv', 'w') as f: 
            write = csv.writer(f) 
            write.writerow(Details) 
            write.writerows(rows) 
        

#        intens_error=np.mean(self.isolist.int_err[:self.last_valid_index])
#        intens_low=self.isolist.intens[self.last_valid_index]-intens_error
#        intens_high=self.isolist.intens[self.last_valid_index]+intens_error
#        sma_low_index=[n for n,i in enumerate(self.isolist.intens[:self.last_valid_index][::-1]) if i>=intens_high][0]
#        sma_low=self.isolist.sma[:self.last_valid_index]
#        sma_low_v=sma_low[-1*sma_low_index]
#        try: 
#            sma_high_index=[ n for n,i in enumerate(self.isolist.intens[self.last_valid_index:]) if i<=intens_low][0]
#            sma_high=self.isolist.sma[self.last_valid_index:]
#            sma_high_v=sma_high[sma_high_index]
#        except: 
#             sma_high_v=-1
#        print("SMA error", sma_low_v, sma_high_v, (sma_high_v-sma_low_v)/2)
#        print("sma error final", self.isolist.sma[self.last_valid_index]-sma_low_v)
	
    
if __name__ == '__main__': 

    parser = argparse.ArgumentParser(description='give table values for .out')
        
    parser.add_argument("-name", help="Name of the file ")
    args = parser.parse_args()
    
    
    sys.stdout=open("additionals.aux","w")
    
    name= args.name
    g_p=galaxy_photometry_aux() 
    g_p.find_read_file() 
    g_p.read_obj()
    g_p.read_params()
    g_p.m_22()
#    g_p.m_20()
    g_p.gen_csv()
    sys.stdout.close()

    
    
