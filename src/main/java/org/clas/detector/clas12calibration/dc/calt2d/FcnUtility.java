/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.dc.calt2d;


/**
 *
 * @author ziegler
 */
public class FcnUtility {
    
    /**
     * 
     * @param alpha
     * @return the reduced angle from 0 to 30 deg.
     */
    public static double getReducedAngle(double alpha) {
        double ralpha = 0;

        ralpha = Math.abs(Math.toRadians(alpha));

        while (ralpha > Math.PI / 3.) {
            ralpha -= Math.PI / 3.;
        }
        if (ralpha > Math.PI / 6.) {
            ralpha = Math.PI / 3. - ralpha;
        }

        return Math.toDegrees(ralpha);
    }  

    /**
     * 
     * @param x
     * @param beta
     * @param distbeta
     * @param v_0
     * @return the beta-dependent timewalk
     */
    public static double getDeltaTimeBeta(double x, double beta, double distbeta, double v_0) {
       
        double value = (0.5*Math.pow(beta*beta*distbeta,3)*x/(Math.pow(beta*beta*distbeta,3)+x*x*x))/v_0;
        
        return value;
    }
    
    
}
