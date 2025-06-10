/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.clas.detector.clas12calibration.dc.calt2d;

//import java.io.ByteArrayOutputStream;
//import java.io.PrintStream;
import org.freehep.math.minuit.FunctionMinimum;
import org.freehep.math.minuit.MnMigrad;
import org.freehep.math.minuit.MnScan;
import org.freehep.math.minuit.MnUserParameters;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.IDataSet;
import org.jlab.groot.fitter.FitterFunction;
import org.jlab.groot.math.F1D;
import org.jlab.groot.math.Func1D;
import org.jlab.groot.math.FunctionFactory;
import org.jlab.groot.math.UserParameter;


/**
 *
 * @author gavalian
 */
public class DataFitter {
    
    public static Boolean FITPRINTOUT = true;
    
    public DataFitter(){
        
    }
    
    public static boolean fit(Func1D func, IDataSet  data, String options, double[] errors){
        
        if(options.contains("Q")==true){
            DataFitter.FITPRINTOUT = false;
        } else {
            DataFitter.FITPRINTOUT = true;
        }
        
        FunctionMinimum scanmin = null;
        FunctionMinimum min = null;
        FitterFunction funcFitter =null;
        try{
	        funcFitter = new FitterFunction(func,
	                data,options);
	        
	        int npars = funcFitter.getFunction().getNPars();
	        
	        MnUserParameters upar = new MnUserParameters();
	        for(int loop = 0; loop < npars; loop++){
	            UserParameter par = funcFitter.getFunction().parameter(loop);
	            upar.add(par.name(),par.value(),errors[loop]);
	            if(par.getStep()<0.0000000001){
	                upar.fix(par.name());
	            }
	            if(par.min()>-1e9&&par.max()<1e9){
	                upar.setLimits(par.name(), par.min(), par.max());
	            }
	        }
	        
	        
	        MnScan  scanner = new MnScan(funcFitter,upar);
	        scanmin = scanner.minimize(); 
                
	        /*
	        System.err.println("******************");
	        System.err.println("*   SCAN RESULTS  *");
	        System.err.println("******************");
	        System.out.println("minimum : " + scanmin);
	        System.out.println("pars    : " + upar);
	        System.out.println(upar);
	        System.err.println("*******************************************");
	        */
	        MnMigrad migrad = new MnMigrad(funcFitter, upar);
	        
	        min = migrad.minimize();
	        func.setFitValid(min.isValid());
                
	        MnUserParameters userpar = min.userParameters();
	        
	        for(int loop = 0; loop < npars; loop++){
	            UserParameter par = funcFitter.getFunction().parameter(loop);
	            par.setValue(userpar.value(par.name()));
	            par.setError(userpar.error(par.name()));
	        }
	        
	        if(options.contains("V")==true){
	            System.out.println(upar);
	            System.err.println("******************");
	            System.err.println("*   FIT RESULTS  *");
	            System.err.println("******************");
	            
	            System.err.println(min);
	        }
	        return min.isValid();
	        //System.out.println(funcFitter.getBenchmarkString());
        }catch(Exception e){
	       // e.printStackTrace();
                return false;
        }
       
    }
    
}
