/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.dc.calt2d;

import java.util.ArrayList;


/**
 *
 * @author ziegler
 */
public class SegmentProperty {

    public int id;
    private ArrayList<Integer> _wires;
    private int _size; 
    private int _sumLR;
    private float _aveWire;
    private int _maxDeltaW  ;
    private int _numWireWithinDW;
    public SegmentProperty() {
        
    }
    
    public SegmentProperty(int sid, ArrayList<Integer> wires, int maxDeltaW) {
        _wires       = wires;
        id           = sid;
        _size = wires.size();
        _maxDeltaW = maxDeltaW;
        this.setProperties(_maxDeltaW);
    }
    
    /**
     * @return the _size
     */
    public int getSize() {
        return _size;
    }

    /**
     * @param _size the _size to set
     */
    public void setSize(int _size) {
        this._size = _size;
    }

    /**
     * @return the _aveWire
     */
    public float getAveWire() {
        return _aveWire;
    }

    /**
     * @param _aveWire the _aveWire to set
     */
    public void setAveWire(float _aveWire) {
        this._aveWire = _aveWire;
    }
    
    /**
     * 
     * @return sum LR
     */
    public int getSumLR() {
        return _sumLR;
    }

    /**
     * 
     * @param _sumLR 
     */
    public void setSumLR(int _sumLR) {
        this._sumLR = _sumLR;
    }

    /**
     * @return the _maxDeltaW
     */
    public int getMaxDeltaW() {
        return _maxDeltaW;
    }

    /**
     * @param _maxDeltaW the _maxDeltaW to set
     */
    public void setMaxDeltaW(int _maxDeltaW) {
        this._maxDeltaW = _maxDeltaW;
    }

    /**
     * @return the _numWireWithinDW
     */
    public int getNumWireWithinDW() {
        return _numWireWithinDW;
    }

    /**
     * @param _numWireWithinDW the _numWireWithinDW to set
     */
    public void setNumWireWithinDW(int _numWireWithinDW) {
        this._numWireWithinDW = _numWireWithinDW;
    }

    private void setProperties(int maxDeltaW) {
        if(this._wires==null || this._wires.size()==0)
            return;
        int avgWire = 0;
        for(int i = 0; i < this._wires.size(); i++) {
            avgWire+=Math.abs(this._wires.get(i));
        }
        avgWire = avgWire/this._wires.size();
        
        int countWithinAveWirNum =0;
        int sumLR=0;
        for(int i = 0; i < this._wires.size(); i++) {
            sumLR+=(int)Math.signum(this._wires.get(i));
            if(Math.abs(Math.abs(this._wires.get(i))-avgWire)<=maxDeltaW)
                countWithinAveWirNum++;
        }
        this.setAveWire(avgWire);
        this.setMaxDeltaW(maxDeltaW);
        this.setNumWireWithinDW(countWithinAveWirNum);
        this.setSumLR(sumLR);
    }
   
}
