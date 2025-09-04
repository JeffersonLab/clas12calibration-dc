/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.clas.detector.clas12calibration.dc.mctuning.analysis.wireineff;
import org.jlab.rec.dc.trajectory.SegmentTrajectory;

/**
 *
 * @author veronique
 */


public class CalSegmentTrajectory extends SegmentTrajectory {


    /**
     * @return the trkBfield
     */
    public double[] getTrkBfield() {
        return trkBfield;
    }

    /**
     * @param trkBfield the trkBfield to set
     */
    public void setTrkBfield(double[] trkBfield) {
        this.trkBfield = trkBfield;
    }

    private double[] trkBfield = new double[6];			// list of Bfield values of trajectory hits for all layers in a superlayer
	
}
