#!/usr/bin/env python3
#Fit surface brightness profiles. Returns the fitting parameters, the chi-square and chi-score.
#Run the command in multiple folders find . -name "intens_fit.csv" -execdir fit_isophote_intens_v2.py {} \;


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import copy
from scipy.stats import chi2
import sys
from astropy.modeling.models import Sersic1D, KingProjectedAnalytic1D

from scipy.optimize import curve_fit


sys.stdout=open("fitting_params_intens_multiple.aux","w")
def mag_star(zp, flux):
    return zp - 2.5 * np.log10(flux)

def f_exp(x,I,a):   
    return I * np.exp(-x/a)

def f_vau(x,I,a): 
    return I * np.exp(-7.67*((np.sign(x/a) *np.abs((x/a))**(1/4)-1)))

def f_Sersic(x,I,a,n): 
    return Sersic1D(amplitude=I, r_eff=a,n=n)(x)
    
def f_King(x,I,a,a2): 
    return KingProjectedAnalytic1D(amplitude = I, r_core = a, r_tide = a2)(x)    

def f_sech(x,I,a):   
    return I/np.cosh(x/a)

def f_csch(x,I,a): 
    return I/np.sinh(x/a)
    
def f_csch_n(x,I,a,n): 
    return I/(np.sinh(x/a)**n)    

def f_gaus(x,I,a):   
   return  I* np.exp(-1/2*(x/a)**2)

def f_gaus2(x,I,a,b=0):   
   return  I* np.exp(-1/2*((x-b)/a)**2)   
      
def combine(func1, func2):     
    def f_(x,I1,a1,I2,a2):
        return func1(x,I1,a1)+func2(x,I2,a2)
    return f_

def plotting():
    dictionary  = {f_exp:"exponential",
                   f_vau:"de Vaucouleurs", 
                   f_sech:"sech",
                   f_csch:"csch",
                   f_gaus:"Gaussian",
                   f_Sersic:"Sersic",
                   f_King:"King",
                   f_csch_n:"$csch^n$ "
                   }
    
    
    fig = plt.figure(figsize=(6,8),dpi=800)
    
    ax = fig.subplots(nrows=1, ncols=1, sharex=False, sharey=False,gridspec_kw={'height_ratios':[1]})
    
    plt.gca().invert_yaxis()
    magnitudes=mag_star(zp,isophots)
    mag_err_list=2*1.08/isophots*err_isophots
    plt.errorbar(sma, magnitudes, yerr=mag_err_list, color='skyblue', alpha=0.7, capsize=2,zorder=1)
    try:
            plt.plot(sma,mag_star(zp, f(sma,*result_copy)), '-',label=dictionary[f1]+"+"+dictionary[f2]+'\n'+r'I1= {:.2f} '.format(result_copy[0]) + 
                                                                                                            r'r1= {:.2f} '.format(result_copy[1])+'\n'+
                                                                                                            r'I2= {:.2f} '.format(result_copy[2]) + 
                                                                                                            r'r2= {:.2f}'.format(result_copy[3]))
            I1,a1,I2,a2=result_copy
            plt.plot(sma,mag_star(zp, f1(sma,I1,a1)), '--',label=dictionary[f1], color='slategray')
            plt.plot(sma,mag_star(zp, f2(sma,I2,a2)), '-.',label=dictionary[f2], color='lightslategray')
    except:
        if f == f_Sersic or f ==f_csch_n:
                plt.plot(sma,mag_star(zp, f(sma,*result_copy)), '-',label=dictionary[f]+'\n'+r'I1= {:.2f} '.format(result_copy[0]) + 
                                                                                              r'r1= {:.2f} '.format(result_copy[1])+
                                                                                              r'n= {:.2f}'.format(result_copy[2]))
        elif f == f_King:
                 plt.plot(sma,mag_star(zp, f(sma,*result_copy)), '-',label=dictionary[f]+'\n'+r'I1= {:.2f} '.format(result_copy[0]) + 
                                                                                               r'r1= {:.2f} '.format(result_copy[1])+
                                                                                               r'r2= {:.2f}'.format(result_copy[2]))
        else:
            plt.plot(sma,mag_star(zp, f(sma,*result_copy)), '-',label=dictionary[f]+'\n'+r'I1= {:.2f} '.format(result_copy[0]) + 
                                                                                            r'r1= {:.2f} '.format(result_copy[1]))
            
    extra=(max(magnitudes) -min(magnitudes) )*0.05
    plt.ylim([max(magnitudes)+extra, min(magnitudes)-extra])
    plt.gca().set_xlim(left=0)
    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()
    ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)))
    ax.set_ylabel('Isophotal magnitude', labelpad=10, fontsize=12)
    ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) 
    res=magnitudes-mag_star(zp,f(sma,*result_copy))    
    ax.legend(loc='upper left',bbox_to_anchor =(0.1, 1))
    plt.title(name_target,loc='left') 
#-----------------------------------------------------------------------------------------------    
    
    ax3 = fig.add_axes([0.125, 0.115, 0.775, 0.1])
    plt.gca().invert_yaxis()
    plt.plot(sma, res,zorder=1,color='steelblue')    
        
    ax3.set_xlabel('sma (")', labelpad=10,fontsize=12)
    ax3.set_ylabel('Fit residual (mag)',fontsize=12)

#-----------------------------------------------------------------------------------------------  

    ax2 = fig.add_axes([0.72, 0.63, 0.15, 0.15])

    # Plot the resulting gaussians
    bins=int(len(res)/10)
    hist3, bins3 = np.histogram(res[~np.isnan(res)],bins=bins)
    centers3 = bins3[:-1] + np.diff(bins3)/2
    yerr3 = np.sqrt(hist3)
    xerr3 = (np.abs(centers3[1])-np.abs(centers3[0]))/2
    min2 = 2*min(np.min(centers3),np.min(centers3))
    max2 = 2*max(np.max(centers3),np.max(centers3))
    x_gauss = np.linspace(min2, max2, 100)
    # Remove zero bins for Chi2 fitting
    x = centers3[ hist3>0 ]
    y = hist3[ hist3>0 ]
    sy = yerr3[ hist3>0 ]

    try:
         popt, pcov= curve_fit(f_gaus2,x,y,p0 =[10,1,0],sigma = sy)
         N3, sig3,m =popt[:]
    except:
        popt, pcov= curve_fit(f_gaus,x,y,p0 =[1,1],sigma = sy)
        N3, sig3 =popt[:]
        m=0
    

    # Reaction time corrected
 #   labelFit = 'Gaussian fit\n' + r'$N$ = {:.2f}'.format(N3) + '\n' + r'$\sigma$ = {:.2f}'.format(sig3)+ '\n' + r'$\mu$ = {:.2f}'.format(m) + '\n' + r'$p_\chi$$_2$ = {:.2f}'.format(chisq)
    ax2.hist(res, bins=bins,color='steelblue', alpha=0.4)
    ax2.step(centers3, hist3, linestyle='-', color='lightblue', where='mid')
 #   ax2.errorbar(centers3, hist3, xerr=xerr3, yerr=yerr3, color='lightblue', drawstyle='steps-mid', linestyle='None', capsize=2)
    ax2.set_xlim([-2.5,2.5])
    y_gauss3 = f_gaus2(x_gauss, *popt)
    ax2.plot(x_gauss, y_gauss3, '--') #alpha=.8, label=labelFit)

    ax2.set_xlabel('Residuals')
    ax2.set_ylabel('Frequency')
    #ax2.legend(loc='upper left', fontsize='x-small', framealpha=0., markerscale=2)
  #  ax2.text(0.4, 0.55, labelFit, transform=plt.gca().transAxes)
#-----------------------------------------------------------------------------------------------
    ax4 = fig.add_axes([0.23, 0.28, 0.15, 0.15])

    # Plot the resulting gaussians
    plt.errorbar(sma, isophots, yerr=err_isophots, color='skyblue', alpha=0.7, capsize=2,zorder=1)
    plt.plot(sma, f(sma,*result_copy), '-')
    extra=(max(isophots) -min(isophots) )*0.05
    plt.ylim([ min(isophots)-extra,max(isophots)+extra])
    plt.gca().set_xlim(left=0)
    xleft, xright = ax.get_xlim()
    ybottom, ytop = ax.get_ylim()
    ax4.set_ylabel('Isophotal intensity') 
    ax4.set_xlabel('sma (")')
    plt.savefig("_fit_func_result_"+f.__name__+"_"+name_target+".png")
    plt.close()
    
    file_out.write(str(chisq_score)+" : "+ str(sig3)+"\n")
#    Ndof = len(t) - 2
#    chi2Score = minuit_lin.fval
#    p_P = chi2.sf(minuit_lin.fval, Ndof)

#    sig_a = minuit_lin.errors[0]
#    sig_b = minuit_lin.errors[1]

#    s1 = r'$\mathbf{Resulting\ fit:}$'
#    s2 = r'Offset = {:.4f}$\pm$'.format(b) + '{:.4f}s'.format(sig_b)
#    s3 = r'Period = {:.4f}$\pm$'.format(2*a) + '{:.4f}s'.format(sig_a)
#    s4 = r'$p(\chi^2={:.2f}$,'.format(chi2Score)
#    s5 = r'$N_{dof}=$' + '{:.2f}) = '.format(Ndof) + '{:.2f}'.format(p_P)

#    ax.text(2,90, s1 + '\n' + s4+s5 + '\n' + s2 + '\n' + s3, family='monospace')

file_out = open('fitting_func_godness.txt', 'w')

   
dataframe = pd.read_csv("intens_fit.csv", header=0, index_col=None)
dataframe1=dataframe.dropna()
sma, isophots,err_mod_isophots,  err_isophots = dataframe1.values.T


names = dataframe.columns
name_target,I_half,a_t,a_half,zp=names[3].split("/")


I_half=float(I_half)
a_t=float(a_t)
a_half=float(a_half)
zp=float(zp)
name_target=name_target.replace('\n','')

list_I=[sum(isophots),I_half,isophots[0],isophots[0]/2,10,1] #np.linspace(2*I_half, 11, 15) #
list_a=[a_t,a_half,1,0.1] #np.linspace(a_t, 0.1, 15) #

list_functions=[f_exp, f_vau, f_sech,f_csch,f_gaus,f_Sersic,f_King,f_csch_n]

#---------------------------------------------------Simgle function fitting---------------------------------------------------------------------------------------------
for f in list_functions: 
    chisq=0.
    var1=100000
    for I in list_I:
        for a in list_a:
            
            start = [I, a]
            bounds=([0.1,0.1],[4*I_half,2*a_t])
            if f == f_Sersic or f ==f_csch_n:
                start = (I, a,1.)
                bounds=[(1.,0.1,0.5),(4*I_half,2*a_t,10.)]
            if f == f_King:
                start = (I, a,a)
                bounds=[(1.,0.1,0.1),(4*I_half,2*a_t,2*a_t)]
                
            try:    
                popt, pcov = curve_fit(f,sma,isophots,sigma = err_isophots,p0 = start,bounds=bounds) #method='trf'
                res = f(sma, *popt)
                chisq_test = np.sum(((isophots-res)/err_isophots*1.5)**2)
                dof=len(sma)-len(popt)
                res_chisq=chi2.sf(chisq_test,dof )
                var=np.var((isophots-res))
            except:
                    var=100001
            if var <= var1:
                        var1=var
                        chisq=res_chisq
                        chisq_score=chisq_test
                        result_copy = copy.deepcopy(popt)
                        result_err= np.sqrt(np.diag(pcov))
                        print(pcov,result_err)

                             

    
#    fig = plt.figure()
#    ax = fig.add_subplot(111)
#    plt.gca().invert_yaxis()
#    plt.errorbar(sma, isophots, yerr=err_mod_isophots, color='skyblue', alpha=0.7, capsize=2,zorder=1)
#    plt.plot(sma, f(sma,*result_copy), '-',)
#    extra=(max(isophots) -min(isophots) )*0.05
#    plt.ylim([max(isophots)+extra, min(isophots)-extra])
#    plt.gca().set_xlim(left=0)
#    xleft, xright = ax.get_xlim()
#    ybottom, ytop = ax.get_ylim()
#    ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)))
#    plt.legend()
    
    try: 
        print("Fitting result:",f.__name__,result_copy,result_err)
        plotting()
        chisq1=sum(((isophots-f(sma, *result_copy))/(err_isophots/1.5))**2)   #TODO CHECK how this error values would need to be calculated. #Plot also  intensities not only magnitudes, check the variance 
        dof=len(sma)-len(result_copy)
        res_chisq=chi2.sf(chisq1,dof )
        print("Function chisq:"+f.__name__+" : ",chisq_score,dof,chisq)  
    except:
        file_out.write("Nan : Nan \n")
        print("\nFAILED\n")
    
#---------------------------------------------------Combined function fitting---------------------------------------------------------------------------------------------  
list_functions2=[f_exp, f_vau, f_sech,f_csch,f_gaus]  
i=0
for f1 in list_functions2:
    i=i+1
    list_functions3= list_functions2[i:]
    if list_functions3:
        for f2 in list_functions3:     
        
            chisqr=0.
            var1=100000.
            for I in list_I:
                for a in list_a:
                    start = [I, a,I,a]
                    bounds=([0.1,0.1,0.1,0.1],[4*I_half,2*a_t,4*I_half,2*a_t])
                    f=combine(f1,f2)  
                    try:                                           
                        popt, pcov = curve_fit(f,sma,isophots,sigma = err_isophots,p0 = start,bounds=bounds)
                        res = f(sma, *popt)
                        chisq_test = np.sum(((isophots-res)/err_isophots/1.5)**2)
                        dof=len(sma)-len(popt)
                        res_chisq=chi2.sf(chisq_test,dof )
                        var=np.var((isophots-res))
                    except:
                        var=100001
                    if var <= var1:
                                var1=var
                                chisq=res_chisq
                                chisq_score=chisq_test
                                result_copy = copy.deepcopy(popt)
                                result_err= np.sqrt(np.diag(pcov))

                    
            
            try:                
                f.__name__=str(f1.__name__+f2.__name__)
                print("Fitting result:",f.__name__,result_copy,result_err)
                print("Function chisq:"+f.__name__+" : ",chisq_score,dof,chisq)
                plotting()
                
            except:
                file_out.write("Nan : Nan \n")
                print("\nFAILED\n")
                  
#            fig = plt.figure()
#            ax = fig.add_subplot(111)
#            plt.gca().invert_yaxis()
#            plt.errorbar(sma, isophots, yerr=err_mod_isophots, color='skyblue', alpha=0.7, capsize=2,zorder=1)
#            plt.plot(sma, f(result_copy.params, sma, 0.,1.), '-',)
#            extra=(max(isophots) -min(isophots) )*0.05
#            plt.ylim([max(isophots)+extra, min(isophots)-extra])
#            plt.gca().set_xlim(left=0)
#            xleft, xright = ax.get_xlim()
#            ybottom, ytop = ax.get_ylim()
#            ax.set_aspect(abs((xright-xleft)/(ybottom-ytop)))
#            plt.legend()
#            plt.savefig("fit_func_result_"+f.__name__+".png")
#                    #plt.show()
#            plt.close()
#            chisq1=sum(f(result_copy.params, sma,isophots,err_mod_isophots)**2)
#            print("Chisq_verification",chisq1,result_copy.chisqr)
#            dof=len(sma)-2. #len(result_copy.params)

        

sys.stdout.close()
file_out.close()
