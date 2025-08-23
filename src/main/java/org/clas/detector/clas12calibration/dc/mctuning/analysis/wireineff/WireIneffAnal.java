/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.clas.detector.clas12calibration.dc.mctuning.analysis.wireineff;

import org.clas.detector.clas12calibration.dc.mctuning.viewer.AnalysisMonitor;
import org.jlab.detector.calib.utils.CalibrationConstants;
import org.jlab.detector.calib.utils.ConstantsManager;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataEvent;
import org.jlab.utils.groups.IndexedList;

/**
 *
 * @author veronique
 */
public class WireIneffAnal extends AnalysisMonitor {
    private final HistogramManager histogramMgr = new HistogramManager();
    private final FitManager fitMgr ;
    private final EventProcessor eventProc = new EventProcessor();
    private OutputManager outputMgr;
    private int[][] totLayA;   // Total layer hits
    private int[][] effLayA;   // Effective layer hits
    public static String variation = "default";
    
    public WireIneffAnal(String name, ConstantsManager ccdb) {
        super(name, ccdb);
        this.setAnalysisTabNames("Inefficiency vs TrackDoca");
        this.init(false, "p0:p1:p2:p3:p4:p5:p6"); 
        fitMgr = new FitManager(this);
    }

    final int nBins = 100;
    @Override
    public void createHistos() {
        histogramMgr.initialize(6, nBins); // 6 superlayers, nBins bins
        totLayA = new int[6][nBins];  // Initialize total layer hits array
        effLayA = new int[6][nBins];  // Initialize effective layer hits array
        DataGroup tr = new DataGroup(6,1);
        for(int i =0; i<6; i++) {
            tr.addDataSet(histogramMgr.getHistogram(i), 0);
            this.getDataGroup().add(tr, 0, i+1, 0);
        }
        
        this.getDataGroup().add(tr, 0,0,0);
        this.getCalib().fireTableDataChanged();
    }
    
    @Override
    public void processEvent(DataEvent ev) {
        eventProc.processEvent(ev, totLayA, effLayA,nBins,this.fitMgr);
    }

    @Override
    public void analysis() {
        histogramMgr.fillFromArrays(totLayA, effLayA);

        // Fit and profile each superlayer
        for (int sl = 0; sl < 6; sl++) {
            fitMgr.runFit(sl, histogramMgr.getProfiles()); // Fit the profile data
            //outputMgr.writeFitParams(sl, fitMgr.getFitParams(sl)); // Write fit parameters to output
            plot(sl);
        }
        fitMgr.updateTable();
    }

    public void plot(int sl) {
        // Plot the profiles and the fits for each superlayer
        GraphErrors profile = histogramMgr.getProfile(sl);
        FitLine fitLine = fitMgr.getFitLine(sl);

        getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").cd(0);
        getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").draw(profile);  // Draw the profile
        getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").draw(fitLine, "same");  // Draw the fit line
        
    }

    public void finish() {
        // Finalize and close output
        if (outputMgr != null) {
            outputMgr.close();
        }
    }

    // Set the output manager
    public void setOutputManager(OutputManager outputMgr) {
        this.outputMgr = outputMgr;
    }

    // Renaming of output file with run number and timestamp
    public void renameOutput(int runNumber, int iteration) {
        if (outputMgr != null) {
            outputMgr.renameOutput(runNumber, iteration);
        }
    }
    
    @Override
    public void constantsEvent(CalibrationConstants cc, int col, int row) {
        String str_sector    = (String) cc.getValueAt(row, 0);
        String str_layer     = (String) cc.getValueAt(row, 1);
        String str_component = (String) cc.getValueAt(row, 2);
        System.out.println(str_sector + " " + str_layer + " " + str_component);
        IndexedList<DataGroup> group = this.getDataGroup();

       int sector    = Integer.parseInt(str_sector);
       int layer     = Integer.parseInt(str_layer);
       int component = Integer.parseInt(str_component);

       if(group.hasItem(sector,layer,component)==true){
           this.plot(layer-1);
       } else {
           System.out.println(" ERROR: can not find the data group");
       }
       
   
    }
}
