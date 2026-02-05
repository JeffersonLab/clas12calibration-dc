/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.clas.detector.clas12calibration.dc.mctuning.analysis.wireineff;

import java.awt.BorderLayout;
import java.awt.Color;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.table.DefaultTableModel;
import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import org.clas.detector.clas12calibration.dc.mctuning.viewer.AnalysisMonitor;
import org.freehep.math.minuit.MnUserParameters;
import org.jlab.detector.calib.utils.CalibrationConstants;
import org.jlab.detector.calib.utils.ConstantsManager;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.ui.TCanvas;
import org.jlab.io.base.DataEvent;
import org.jlab.utils.groups.IndexedList;

/**
 *
 * @author veronique
 **/
public class WireIneffAnal extends AnalysisMonitor {
    
    //public static double[] Bsq = new double[] {0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.25,5.75,6.25}; 
    public static double[] Bsq = new double[] {0.5,1.5,2.5,3.5,4.5};
    //public static double[] Bsq = new double[] {0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25};
    private final HistogramManager histogramMgr = new HistogramManager();
    private final FitManager fitMgr ;
    private final EventProcessor eventProc = new EventProcessor();
    private OutputManager outputMgr;
    private int[][][] totLayA;   // Total layer hits
    private int[][][] effLayA;   // Effective layer hits
    private int[][][] totLayB;   // Total layer hits
    private int[][][] effLayB;   // Effective layer hits
    public static String variation = "default";
    
    public WireIneffAnal(String name, ConstantsManager ccdb) {
        super(name, ccdb);
        this.setAnalysisTabNames("Inefficiency vs TrackDoca", 
                "Sector1 Efficiency / Layer","Sector2 Efficiency / Layer","Sector3 Efficiency / Layer",
                "Sector4 Efficiency / Layer","Sector5 Efficiency / Layer","Sector6 Efficiency / Layer");
        this.init(false, "p0:p1:p2:p3:p4:p5:p6:p7:p8:p9:p10:p11"); 
        fitMgr = new FitManager(this);
    }

    public static final int nBins = 50;//100;
    @Override
    public void createHistos() {
        histogramMgr.initialize(6, nBins); // 6 superlayers, nBins bins
        totLayA = new int[6][nBins][Bsq.length];  // Initialize total layer hits array
        effLayA = new int[6][nBins][Bsq.length];  // Initialize effective layer hits array
        totLayB = new int[6][6][6];  // Initialize  hits array
        effLayB = new int[6][6][6];  // Initialize effective  hits array
        DataGroup tr = new DataGroup(6,2);
        for(int i =0; i<6; i++) {
            tr.addDataSet(histogramMgr.getNormalizedHistogram(i), 0);
            if(i==2 || i==3) {
                for(int bbin = 0; bbin<Bsq.length; bbin++) {
                    tr.addDataSet(histogramMgr.getBHistogram(i-2,bbin), 0);
                }
            }
            this.getDataGroup().add(tr, 0, i+1, 0); 
        }
        
        this.getDataGroup().add(tr, 0,0,0);
        this.getCalib().fireTableDataChanged();
    }
    
    @Override
    public void processEvent(DataEvent ev) {
        eventProc.processEvent(ev, totLayA, effLayA,nBins,this.fitMgr, 
                totLayB, effLayB);
    }
    final Map<Coordinate, FitLine2D> fitLines2d = new HashMap<>();
    final Map<Coordinate, FitLine2D> fitLines2dcolor = new HashMap<>();
    @Override
    public void analysis() {
        histogramMgr.fillFromArrays(totLayA, effLayA, totLayB, effLayB);
        
        // Fit and profile each superlayer
        for (int sl = 0; sl < 2; sl++) {
            fitMgr.runFit(sl, histogramMgr.getNormalizedProfiles()); // Fit the profile data
            fitMgr.runFit(sl+4, histogramMgr.getNormalizedProfiles()); // Fit the profile data
            //outputMgr.writeFitParams(sl, fitMgr.getFitParams(sl)); // Write fit parameters to output
            //plot(sl);
//            if(sl==2 || sl==3) {
//                for(int bbin = 0; bbin<Bsq.length; bbin++) {
//                    H1F hb = histogramMgr.getBHistogram(sl-2, bbin);
//                    if(hb.getIntegral()>0) { System.out.println("OK to fit....");
//                        GraphErrors bprofile = histogramMgr.getBProfile(sl-2, bbin);
//                        MnUserParameters mnp = fitMgr.r2Fit(sl,bprofile);
//                        
//                        for(int p =0; p<6; p++) { 
//                             histogramMgr.getParsvvsB().get(new Coordinate(sl-2,p)).addPoint(Bsq[bbin], mnp.value(p+1), 0, mnp.error(p+1));
//                        }
//                    }
//                }
//            }
        }
        for(int sli = 0; sli<2; sli++) {
            fitMgr.run2DFit(sli, histogramMgr.getN2dhistograms(), 
                histogramMgr.getN2dhistogramserrs());
            MnUserParameters params = fitMgr.getFitParams2d(sli);
            for(int j = 0; j<Bsq.length; j++) {
                fitLines2d.put(new Coordinate(sli, j), 
                        new FitLine2D("f2d"+sli+j, sli, j, params));
                fitLines2d.get(new Coordinate(sli, j)).setLineWidth(4);
                if(j==0) fitLines2d.get(new Coordinate(sli, j)).setLineColor(0);
                fitLines2dcolor.put(new Coordinate(sli, j), 
                        new FitLine2D("f2d"+sli+j, sli, j, params));
                fitLines2dcolor.get(new Coordinate(sli, j)).setLineStyle(3);
                fitLines2dcolor.get(new Coordinate(sli, j)).setLineWidth(4);
                fitLines2dcolor.get(new Coordinate(sli, j)).setLineColor(j+1);
                
            }
            
        } 
        
        fitMgr.updateTable();
        TCanvas can1 = new TCanvas("SL 3 B-field dependence of parameters", 800, 800);
        TCanvas can2 = new TCanvas("SL 4 B-field dependence of parameters", 800, 800);
        
        can1.cd(0);
        can1.draw(histogramMgr.getN2dhistograms().get(new Coordinate(0)));
        
   
        can2.cd(0);
        can2.draw(histogramMgr.getN2dhistograms().get(new Coordinate(1)));
        
      /*  
        for(int p =0; p<6; p++) {
            TCanvas can01 = new TCanvas("SL 3 B-field dependence of parameters", 800, 800);
            can01.cd(0);
            can01.draw(histogramMgr.getParsvvsB().get(new Coordinate(0,p)));
        }
         for(int p =0; p<6; p++) {
            TCanvas can02 = new TCanvas("SL 4 B-field dependence of parameters", 800, 800);
            can02.cd(0);
            can02.draw(histogramMgr.getParsvvsB().get(new Coordinate(1,p)));
        }  
        */
        fitMgr.updateTable();
        for(int s =0; s< 6; s++){
            String stg="Sector"+(s+1)+" Efficiency / Layer";
            getAnalysisCanvas().getCanvas(stg).divide(2, 3);
            for(int sl =0; sl< 6; sl++){
                getAnalysisCanvas().getCanvas(stg).cd(sl);
                double min = histogramMgr.getN1dhistograms().get(new Coordinate(s,sl)).getMin()-5;
                double max = histogramMgr.getN1dhistograms().get(new Coordinate(s,sl)).getMax()+5;
                if(max>103)
                    max =103;
                getAnalysisCanvas().getCanvas(stg).getPad(sl).getAxisY().setRange(min, max);
                getAnalysisCanvas().getCanvas(stg).draw(histogramMgr.getN1dhistograms().get(new Coordinate(s,sl)), "E");
                        
            }
        }
    }
    
    private boolean reRunFit(int sl, int[] fixedPar) {
        if(sl<2 || sl>3) {
            fitMgr.runFit(sl, histogramMgr.getNormalizedProfiles(), fixedPar, false); // Fit the profile data
        } else {
            int sli = sl-2;
            fitMgr.run2DFit(sli, histogramMgr.getN2dhistograms(), 
                histogramMgr.getN2dhistogramserrs(), fixedPar, false);
            MnUserParameters params = fitMgr.getFitParams2d(sli);
            for(int j = 0; j<Bsq.length; j++) {
                fitLines2d.put(new Coordinate(sli, j), 
                        new FitLine2D("f2d"+sli+j, sli, j, params));
                fitLines2d.get(new Coordinate(sli, j)).setLineWidth(4);
                if(j==0) fitLines2d.get(new Coordinate(sli, j)).setLineColor(0);
                fitLines2dcolor.put(new Coordinate(sli, j), 
                        new FitLine2D("f2d"+sli+j, sli, j, params));
                fitLines2dcolor.get(new Coordinate(sli, j)).setLineStyle(3);
                fitLines2dcolor.get(new Coordinate(sli, j)).setLineWidth(4);
                fitLines2dcolor.get(new Coordinate(sli, j)).setLineColor(j+1);
                
            }
        } 
        
        fitMgr.updateTable();
        
        return fitMgr.isFitValid;
    }
    
    public void openFitParsPanel(int li) {
    // Example fit parameters (replace with your actual values)
        String[] columnNames = {"Parameter", "Value", "Error", "Min", "Max", "Fix"};
        String[] paramNames = fitMgr.parNames2d;
         // p0 ... p11
        if(li<2 || li>3) paramNames = fitMgr.parNames;
        
        int nPars = paramNames.length; 
        Object[][] data = new Object[nPars][6];

        for (int i = 0; i < nPars; i++) {
            String parName = paramNames[i];
            double value   = 0;
            double error   = 0;
            if(li<2 || li>3) {
                value = fitMgr.getFitParams(li).value(i); 
                error = fitMgr.getFitParams(li).error(i); 
            } else {
                value = fitMgr.getFitParams2d(li-2).value(i); 
                error = fitMgr.getFitParams2d(li-2).error(i); 
            }
            double min     = value - 3*error;   // default (adjust per param if needed)
            double max     = value + 3*error;
            boolean fixed  = false;

            data[i][0] = parName;
            data[i][1] = value;
            data[i][2] = error;
            data[i][3] = min;
            data[i][4] = max;
            data[i][5] = fixed;
        }
        // Create table model allowing editing of ranges and fix column
        DefaultTableModel model = new DefaultTableModel(data, columnNames) {
            @Override
            public boolean isCellEditable(int row, int col) {
                // Value & Error not editable, only Min, Max, Fix
                return col == 1 || col == 2 || col == 3 || col == 4 || col == 5;
            }

            @Override
            public Class<?> getColumnClass(int col) {
                if (col == 5) return Boolean.class;
                if (col == 1 || col == 2 || col == 3 || col == 4) return Double.class;
                return String.class;
            }
        };

        JTable table = new JTable(model);
        table.setFillsViewportHeight(true);

        JScrollPane scrollPane = new JScrollPane(table);
         // ✅ Fit convergence info
        JTextArea convergenceInfo = new JTextArea("Fit status: not run yet");
        convergenceInfo.setEditable(false);
        convergenceInfo.setLineWrap(true);
        convergenceInfo.setWrapStyleWord(true);
        convergenceInfo.setBackground(new Color(240, 240, 240));
        convergenceInfo.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));

        // ✅ Re-run button
        // ================
        JButton rerunButton = new JButton("Re-run Fit");
        rerunButton.addActionListener(e -> {
            // Extract updated parameters
            List<Integer> fixedParsList = new ArrayList<>();
            for (int i = 0; i < model.getRowCount(); i++) {
                String name   = (String) model.getValueAt(i, 0);
                double value  = (double) model.getValueAt(i, 1);
                double error  = (double) model.getValueAt(i, 2);
                double min    = (double) model.getValueAt(i, 3);
                double max    = (double) model.getValueAt(i, 4);
                boolean fixed = (boolean) model.getValueAt(i, 5);
                // Push updated values into your fit manager
                //
                if(li<2 || li>3) {
                    fitMgr.getFitParams(li).setValue(name, value);
                    fitMgr.getFitParams(li).setError(name, error);
                    fitMgr.getFitParams(li).setLimits(name, min, max);
                    if(fixed) {
                        fitMgr.getFitParams(li).fix(name);
                        fixedParsList.add(i);
                    }
                } else {
                    fitMgr.getFitParams2d(li-2).setValue(name, value);
                    fitMgr.getFitParams2d(li-2).setError(name, error);
                    fitMgr.getFitParams2d(li-2).setLimits(name, min, max);
                    if(fixed) {
                        fitMgr.getFitParams2d(li-2).fix(name);
                        fixedParsList.add(i);
                    }
                }
                //
                
                System.out.printf("Param %s: val=%.3f ± %.3f, range=[%.2f, %.2f], fixed=%b%n",
                        name, value, error, min, max, fixed);
            }
            int[] fixedPars;
            if(!fixedParsList.isEmpty()) {
                fixedPars = new int[fixedParsList.size()];
                for(int ii =0; ii<fixedParsList.size(); ii++) {
                    fixedPars[ii]=fixedParsList.get(ii);
                }
            } else {
                fixedPars= new int[]{-1};
            }
            // Run fit and update convergence info
            
            boolean converged = this.reRunFit(li, fixedPars); 
            for (int i = 0; i < model.getRowCount(); i++) {
                String name   = (String) model.getValueAt(i, 0);
                boolean fixed = (boolean) model.getValueAt(i, 5);
                if(fixed) {
                    if(li<2 || li>3) {
                        fitMgr.getFitParams(li).release(name);
                    } else {
                        fitMgr.getFitParams2d(li-2).release(name);
                    }
                }
            }
            if (converged) {
                convergenceInfo.setText("✅ Fit converged successfully.");
            } else {
                convergenceInfo.setText("❌ Fit did NOT converge. Try adjusting ranges or fixing parameters.");
            }
        });

        JPanel panel = new JPanel(new BorderLayout());
        panel.add(scrollPane, BorderLayout.CENTER);

        JPanel bottomPanel = new JPanel(new BorderLayout());
        bottomPanel.add(convergenceInfo, BorderLayout.CENTER);
        bottomPanel.add(rerunButton, BorderLayout.EAST);

        panel.add(bottomPanel, BorderLayout.SOUTH);
        
        JFrame frame = new JFrame("Fit Parameters (Layer " + (li+1) + ")");
        frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        frame.setSize(700, 400);
        frame.setLocationRelativeTo(null);
        frame.add(panel);
        frame.setVisible(true);
        
    }
    public void plot(int sl) {
        // Plot the profiles and the fits for each superlayer
        GraphErrors nprofile = histogramMgr.getNormalizedProfile(sl);
        GraphErrors profile = histogramMgr.getProfile(sl);
        
        getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").clear();
        if(sl==2 || sl==3) {
            getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").divide(2, 2);
        } else {
            getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").divide(2, 1);
        }
        getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").cd(0);
        getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").draw(nprofile);  // Draw the profile
        if(sl<2 || sl>3) {
            FitLine fitLine = fitMgr.getFitLine(sl);
            getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").draw(fitLine, "same");  // Draw the fit line
        }
        
        if(sl==2 || sl==3) {
            getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").cd(2);
            for(int bbin = 0; bbin<Bsq.length; bbin++) {
                H1F hb = histogramMgr.getBHistogram(sl-2, bbin);
                if(hb.getIntegral()>0) {
                    GraphErrors bprofile = histogramMgr.getBProfile(sl-2, bbin);
                    getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").draw(bprofile, "same");
                    getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").draw(fitLines2d.get(new Coordinate(sl-2, bbin)), "same");
                    getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").draw(fitLines2dcolor.get(new Coordinate(sl-2, bbin)), "same");
                }
            }
        }
            
        getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").cd(1);
        getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").draw(profile);  // Draw the profile
        if(sl==2 || sl==3) {
            getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").cd(3);
            for(int bbin = 0; bbin<Bsq.length; bbin++) {
                H1F hb = histogramMgr.getBHistogram(sl-2, bbin);
                if(hb.getIntegral()>0) {
                    GraphErrors bprofile = histogramMgr.getBProfileU(sl-2, bbin);
                    getAnalysisCanvas().getCanvas("Inefficiency vs TrackDoca").draw(bprofile, "same");
                }
            }
        }
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
    private long lastClickTime = 0;                 // store last click time
    private static final int DOUBLE_CLICK_DELAY = 300; // ms threshold
    @Override
    public void constantsEvent(CalibrationConstants cc, int col, int row) {
        long currentTime = System.currentTimeMillis();
        String str_sector    = (String) cc.getValueAt(row, 0);
        String str_layer     = (String) cc.getValueAt(row, 1);
        String str_component = (String) cc.getValueAt(row, 2);
        int sector    = Integer.parseInt(str_sector);
        int layer     = Integer.parseInt(str_layer);
        int component = Integer.parseInt(str_component);
        System.out.println(str_sector + " " + str_layer + " " + str_component);
        IndexedList<DataGroup> group = this.getDataGroup();

        if (currentTime - lastClickTime <= DOUBLE_CLICK_DELAY) {
            // Handle double click
            System.out.println("Double click detected!");
            this.openFitParsPanel(layer-1);
            // Put your double-click behavior here
            // e.g., open a new window, show details, etc.
        } else {
            // Handle single click
            
            if (group.hasItem(sector, layer, component)) {
                this.plot(layer - 1);
            } else {
                System.out.println(" ERROR: cannot find the data group");
            }
        }

        lastClickTime = currentTime; // update last click time
    }
    
}
