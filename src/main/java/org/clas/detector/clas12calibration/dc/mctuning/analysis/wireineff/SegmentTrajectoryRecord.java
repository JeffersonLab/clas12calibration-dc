/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.clas.detector.clas12calibration.dc.mctuning.analysis.wireineff;

/**
 *
 * @author veronique
 */
public class SegmentTrajectoryRecord {
    private int segmentID;
    private int sector;
    private int superlayer;
    private int layer;
    private int matchedHitID;
    private double trkDoca;
    private double trkBfield;

    // constructor + getters
    public SegmentTrajectoryRecord(int segmentID, int sector, int superlayer, int layer,
                                   int matchedHitID, double trkDoca, double trkBfield) {
        this.segmentID = segmentID;
        this.sector = sector;
        this.superlayer = superlayer;
        this.layer = layer;
        this.matchedHitID = matchedHitID;
        this.trkDoca = trkDoca;
        this.trkBfield = trkBfield;
    }

    /**
     * @return the segmentID
     */
    public int getSegmentID() {
        return segmentID;
    }

    /**
     * @param segmentID the segmentID to set
     */
    public void setSegmentID(int segmentID) {
        this.segmentID = segmentID;
    }

    /**
     * @return the sector
     */
    public int getSector() {
        return sector;
    }

    /**
     * @param sector the sector to set
     */
    public void setSector(int sector) {
        this.sector = sector;
    }

    /**
     * @return the superlayer
     */
    public int getSuperlayer() {
        return superlayer;
    }

    /**
     * @param superlayer the superlayer to set
     */
    public void setSuperlayer(int superlayer) {
        this.superlayer = superlayer;
    }

    /**
     * @return the layer
     */
    public int getLayer() {
        return layer;
    }

    /**
     * @param layer the layer to set
     */
    public void setLayer(int layer) {
        this.layer = layer;
    }

    /**
     * @return the matchedHitID
     */
    public int getMatchedHitID() {
        return matchedHitID;
    }

    /**
     * @param matchedHitID the matchedHitID to set
     */
    public void setMatchedHitID(int matchedHitID) {
        this.matchedHitID = matchedHitID;
    }

    /**
     * @return the trkDoca
     */
    public double getTrkDoca() {
        return trkDoca;
    }

    /**
     * @param trkDoca the trkDoca to set
     */
    public void setTrkDoca(double trkDoca) {
        this.trkDoca = trkDoca;
    }

    /**
     * @return the trkBfield
     */
    public double getTrkBfield() {
        return trkBfield;
    }

    /**
     * @param trkBfield the trkBfield to set
     */
    public void setTrkBfield(double trkBfield) {
        this.trkBfield = trkBfield;
    }
}
