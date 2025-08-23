/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.dc.caltdccuts;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import org.clas.detector.clas12calibration.dc.analysis.TDCParamsPanel;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.field;
import org.clas.detector.clas12calibration.viewer.AnalysisMonitor;
import org.jlab.detector.calib.utils.CalibrationConstants;
import org.jlab.detector.calib.utils.ConstantsManager;
import org.jlab.groot.data.DataLine;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.utils.groups.IndexedList;
import org.jlab.utils.groups.IndexedTable;

/**
 *
 * @author ziegler
 */
public class TDCCuts extends AnalysisMonitor {
    public int runNumber;
    private String[] analTabs = new String[]{"Raw TDC", "TDC vs wire"};
    public boolean analysisDone = false;
    public TDCCuts(String name, ConstantsManager ccdb) {
        super(name, ccdb);
        this.setAnalysisTabNames(analTabs[0], analTabs[1]);
        this.init(false, "TDC upper cut");
       
    }
    
    public static Map<Coordinate, H1F> TDCHis          = new HashMap<Coordinate, H1F>(); 
    public static Map<Coordinate, H1F> TDCHisCoarse    = new HashMap<Coordinate, H1F>(); 
    public static Map<Coordinate, H2F> TDCvsWire       = new HashMap<Coordinate, H2F>(); 
    int nwir  = 112;
    int reg = 3;
    public static final double[] tLow4T0Fits  = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    public static final double[] tHigh4T0Fits  = {500.0, 1500.0, 1200.0}; 
    public static final double[] fitbound = {170, 280, 400};
    public static final double buffer = 50;
    public static final int histRebinFac=5;
    
    @Override
    public void createHistos() { 
        // initialize canvas and create histograms
        this.setNumberOfEvents(0);
        DataGroup hgrps = new DataGroup(6,7);
        String hNm;
        String hTtl;
        int ijk=-1;
        for (int i = 0; i < reg; i++) {
            TDCvsWire.put(new Coordinate(i), new H2F("h2", "TDC vs wire", 112, 0.5, 112.5, (int)tHigh4T0Fits[i]/histRebinFac, tLow4T0Fits[i], tHigh4T0Fits[i]) );
            TDCHis.put(new Coordinate(i,0), new H1F("h", (int)tHigh4T0Fits[i]/histRebinFac/10, tLow4T0Fits[i], tHigh4T0Fits[i])); 
                
            for (int j = 0; j < this.getWireMax(i); j++) {
                hNm = String.format("tdcS%dS%d", i + 1, j + 1);
                TDCHisCoarse.put(new Coordinate(i,j), new H1F(hNm, (int)tHigh4T0Fits[i]/histRebinFac/10, tLow4T0Fits[i], tHigh4T0Fits[i])); 
                hgrps.addDataSet(TDCHisCoarse.get(new Coordinate(i, j)), 0);
                this.getDataGroup().add(hgrps, 0, i+1, j+1);
            }
        }
        
        
        this.getDataGroup().add(hgrps,0,0,0);
        for (int i = 0; i < reg; i++) {
            for (int j = 0; j < this.getWireMax(i); j++) {
                this.getCalib().addEntry(0,i+1,j+1);
            }
        }    
        this.getCalib().fireTableDataChanged();
    }
    int count = 0;
    public static Map<Integer,DataLine[]> cuts = new HashMap<>();
    public static Map<Integer, double[]> cutParams = new HashMap<>();
    private static double[] getR2Parameters(double X1, double Y1, double X2, double Y2, double wiremax,
            boolean oldfcn) {
        double Ca= wiremax;
        if(!oldfcn) Ca/=56;
        double C = -56*(Y2-Y1)/(X2-X1);
        double Cb = (X2*Y1-X1*Y2)/(X2-X1);
        double M = -Ca*C+Cb; 
        return new double[] {M, C};
    }
    public static void readTDCCuts(IndexedTable tdccuts,
            boolean oldfcn) {
        double Bscale=1;
        if(field!=0) Bscale =field*field;
        for(int r =0; r<3; r++) {
            int region = r+1;
            if(region==1 || region==3) {
                double timeCutMin = tdccuts.getIntValue("MinEdge", 0, region, 0);
                double timeCutMax = tdccuts.getIntValue("MaxEdge", 0, region, 0);
                cutParams.put(r, new double[]{timeCutMin, timeCutMax});
                DataLine rmin = new DataLine(1, timeCutMin,  112, timeCutMin);
                DataLine rmax = new DataLine(1, timeCutMax,  112, timeCutMax);
                cuts.put(r, new DataLine[]{rmin,rmax});
            } else {
                double timeCutLC1 = tdccuts.getIntValue("LinearCoeff", 0, region, 1);
                double timeCutMin1 = tdccuts.getIntValue("MinEdge", 0, region, 1);
                double timeCutMax1 = tdccuts.getIntValue("MaxEdge", 0, region, 1);
                double timeCutLC56 = tdccuts.getIntValue("LinearCoeff", 0, region, 56);
                double timeCutMin56 = tdccuts.getIntValue("MinEdge", 0, region, 56);
                double timeCutMax56 = tdccuts.getIntValue("MaxEdge", 0, region, 56);
                //double[] newPars1=getR2Parameters(1, 1300, 55, 600, 56);
                //timeCutMax1=newPars1[0];
                //timeCutLC1=newPars1[1];
                
                //double[] newPars56=getR2Parameters(56, 600, 112, 400, 112);
                //timeCutMax56=newPars56[0];
                //timeCutLC56=newPars56[1];
                
                cutParams.put(r, new double[]{timeCutMin1, timeCutMax1,timeCutLC1,
                                                        timeCutMin56, timeCutMax56,timeCutLC56});
                
                
                double R2Xi56[] = new double[] {56, 112};
                double R2Yi56[] = new double[2] ;
                for(int xi=0; xi<2; xi++) {
                    if(oldfcn) {
                        R2Yi56[xi] = timeCutMax56 + timeCutLC56 * (double) (112 - (double)R2Xi56[xi] / (double) 56) * Bscale;
                    } else {
                        R2Yi56[xi] = timeCutMax56 + timeCutLC56 * (double) (112 - (double)R2Xi56[xi] )/ (double) 56 * Bscale;
                    } 
                }
                double R2Xi1[] = new double[] {1, 55};
                double R2Yi1[] = new double[2] ;
                for(int xi=0; xi<2; xi++) {
                    if(oldfcn) {
                        R2Yi1[xi] = timeCutMax1 + timeCutLC1 * (double) (56 - (double) R2Xi1[xi] / (double) 56) * Bscale ;
                    } else {
                        R2Yi1[xi] = timeCutMax1 + timeCutLC1 * (double) (56 - (double) R2Xi1[xi] )/ (double) 56 * Bscale ;
                    }
                    
                }
                DataLine r2w56U = new DataLine(R2Xi56[0], R2Yi56[0], R2Xi56[1], R2Yi56[1]);
                DataLine r2w1U  = new DataLine(R2Xi1[0], R2Yi1[0], R2Xi1[1], R2Yi1[1]);
                DataLine r2w56L = new DataLine(56, timeCutMin56, 112, timeCutMin56);
                DataLine r2w1L = new DataLine(1, timeCutMin1, 56, timeCutMin1);
                cuts.put(r, new DataLine[]{r2w1L,r2w56L,r2w1U,r2w56U});
            }
        }
        for(Integer i : cuts.keySet()) {
            for(int j=0; j<cuts.get(i).length; j++) {
                cuts.get(i)[j].setLineColor(5);
                cuts.get(i)[j].setLineStyle(4);
                cuts.get(i)[j].setLineWidth(3);
            }
        }
    }
    
    @Override
    public void processEvent(DataEvent event) {
        
        if (!event.hasBank("RUN::config")) {
            return ;
        }
        
        DataBank bank = event.getBank("RUN::config");
        int newRun = bank.getInt("run", 0);
        if (newRun == 0) {
           return ;
        } else {
           count++;
        }
        
        if(count==1) {
            
            runNumber = newRun; 
            field = event.getBank("RUN::config").getFloat("torus",0);
            //IndexedTable tab = TDCViewer.ccdb.getConstants(newRun, 
            //        "/calibration/dc/time_corrections/tdctimingcuts");
            //readTDCCuts(tab, oldfcn);
        }
        
        // get hits property
        if(event.hasBank("DC::tdc") 
                && !(event.hasBank("TimeBasedTrkg::TBHits") )
                && !(event.hasBank("HitBasedTrkg::HBHits") )
                ) {
            DataBank bnkHits = event.getBank("DC::tdc");

            for (int j = 0; j < bnkHits.rows(); j++) {

                //int sec = bnkHits.getInt("region", j);
                int sl = (bnkHits.getInt("layer", j) - 1)/6 + 1;//  goes from 1 to 36 in data
                int tdc = bnkHits.getInt("TDC", j);
                int wire = bnkHits.getInt("component", j);
                //int sec = bnkHits.getInt("sector", j);
                int region = (sl + 1) / 2;
                TDCvsWire.get(new Coordinate(region-1)).fill(wire, (float) tdc);
                this.TDCHis.get(new Coordinate(region-1,0))
                        .fill((float)tdc); 
                if(sl<3 || sl>4) {
                    this.TDCHisCoarse.get(new Coordinate(region-1,0))
                        .fill((float)tdc); 
                } else {
                    this.TDCHisCoarse.get(new Coordinate(region-1,wire-1))
                        .fill((float)tdc);
                }
            } 
        }
        if(event.hasBank("TimeBasedTrkg::TBHits")) {
            DataBank bnkHits = event.getBank("TimeBasedTrkg::TBHits");

            for (int j = 0; j < bnkHits.rows(); j++) {

                int sl = bnkHits.getInt("superlayer", j);
                int tdc = bnkHits.getInt("TDC", j);
                int wire = bnkHits.getInt("wire", j);
                //int sec = bnkHits.getInt("sector", j);
                int region = (sl + 1) / 2;
                TDCvsWire.get(new Coordinate(region-1)).fill(wire, (float) tdc);
                this.TDCHis.get(new Coordinate(region-1,0))
                        .fill((float)tdc); 
                if(sl<3 || sl>4) {
                    this.TDCHisCoarse.get(new Coordinate(region-1,0))
                        .fill((float)tdc); 
                } else {
                    this.TDCHisCoarse.get(new Coordinate(region-1,wire-1))
                        .fill((float)tdc);
                }
            } 
        }
    }
    public void FillCCDBTable(int i, int j, PrintWriter writer) {
        H1F h = this.TDCHisCoarse.get(new Coordinate(i, j));
        double cbound = getCorrectedBound(i,j);
        double obound = getLowBound(i, j);
        
        if(i==0 || i==2) {
            highbound[i]=cbound;
            for(int k = 0; k<112; k++) {
                writer.printf("%d\t %d\t %d\t %.0f\t %.0f\n",
                              0, i+1, k+1, obound, cbound);
            }
        } else {
            writer.printf("%d\t %d\t %d\t %.0f\t %.0f\n",
                              0, i+1, j+1, obound, cbound);
        }
         
    }
    public void Plot(int i , int j) {
        this.getAnalysisCanvas().getCanvas(analTabs[0]).cd(0);
        this.getAnalysisCanvas().getCanvas(analTabs[0]).clear();
        this.getAnalysisCanvas().getCanvas(analTabs[0])
                .draw(this.TDCHisCoarse.get(new Coordinate(i, j)));
//        double[] bounds = new double[2];
//        if(cuts.containsKey(i)) {
//            if(i!=1) {
//                bounds[0] = cuts.get(i)[0].getOriginY();
//                bounds[1] = cuts.get(i)[1].getOriginY();
//            } else {
//                bounds[0] = cuts.get(i)[0].getOriginY();
//                
//                if(j<56) {
//                    bounds[1] = getBound(j, cuts.get(i)[2].getOriginX(),cuts.get(i)[2].getOriginY(),
//                                   cuts.get(i)[2].getEndX(),cuts.get(i)[2].getEndY());
//                } else {
//                    bounds[1] = getBound(j, cuts.get(i)[3].getOriginX(),cuts.get(i)[3].getOriginY(), 
//                                   cuts.get(i)[3].getEndX(),cuts.get(i)[3].getEndY());
//                }
//            }
//        }
        //H1F h = this.TDCHisCoarse.get(new Coordinate(i, j));
        double max = this.TDCHisCoarse.get(new Coordinate(i, j)).getMax()/4; 
        double cbound = getCorrectedBound(i,j);
        DataLine lhc = new DataLine(cbound,0, cbound, max);
        lhc.setLineColor(2);
        this.getAnalysisCanvas().getCanvas(analTabs[0])
                .draw(lhc);
        double obound = getLowBound(i, j);
         DataLine ll = new DataLine(obound,0, obound, max);
        ll.setLineColor(2);
        this.getAnalysisCanvas().getCanvas(analTabs[0])
                .draw(ll);
        
        this.getCalib().setDoubleValue(cbound, "TDC upper cut", 0, i+1, j+1);
        this.getCalib().fireTableDataChanged();
        
    }
    static double getBound(int j, double originX, double originY, double endX, double endY) {
        double sl = (endY-originY)/(endX-originX);
        double i = originY-sl*originX;
        
        return (double) j *sl +i;
    }
    
    
    public static int R1NSigs=4;
    private double getCorrectedBound(int i, int ji) {
        if(i==0) 
            return CutFinderTool.findR1UpperEdge(this.TDCHisCoarse.get(new Coordinate(i, ji)), R1NSigs);
        int j = ji;
        if(i==1 && j<2) j=2;
        H1F h = this.TDCHisCoarse.get(new Coordinate(i, j));
        double binWidth = h.getDataX(1)-h.getDataX(0);
        
        int g1 = this.getLargestGradiant(h,binWidth);
        return h.getDataX(g1+3)+binWidth/2;
    }
    
    private void setLowBound(int i) {
        H1F h = this.TDCHis.get(new Coordinate(i, 0));
        lowbound[i]= CutFinderTool.findLowerEdge(h, 3);
    }
    
    double[] lowbound = new double[3];
    private double getLowBound(int i, int j) {
        return lowbound[i];
    }
    double[] highbound = new double[3];
    @Override
    public void plotHistos() {
        this.getAnalysisCanvas().getCanvas(analTabs[0]).setGridX(false);
        this.getAnalysisCanvas().getCanvas(analTabs[0]).setGridY(false);
        this.getAnalysisCanvas().getCanvas(analTabs[0]).update();

        
    }
    @Override
    public void timerUpdate() {
    }
    @Override
    public void analysis() {
        String filePath = "Files/ccdb_TDC_run" + this.runNumber  + ".txt";
        File file = new File(filePath);
        if (file.exists()) {
            // Delete the file
            if (file.delete()) {
                System.out.println("File deleted successfully.");
            } else {
                System.out.println("Failed to delete the file.");
            }
        } else {
            System.out.println("File does not exist.");
        }
        try (PrintWriter writer = new PrintWriter(new FileWriter(filePath))) {
            // Write the header
            writer.printf("#& sector\t region\t component\t min\t max\n");
            
            for (int i = 0; i < reg; i++) {
                this.setLowBound(i);
                for (int j = 0; j < this.getWireMax(i); j++) {
                    this.FillCCDBTable(i, j, writer);
                }
            }
            writer.close();
        } catch (IOException ex) { 
            Logger.getLogger(TDCCuts.class.getName()).log(Level.SEVERE, null, ex);
        }
        GraphErrors g = new GraphErrors();
        g.setMarkerSize(4);
        for (int i = 0; i < reg; i++) {
            for (int j = 0; j < this.getWireMax(i); j++) {
                this.Plot(i, j);
                if(i==1) {
                    double cbound = getCorrectedBound(i, j);
                    g.addPoint(j+1, cbound, 0, 0);
                }
            }
        }
       System.out.println("Output successfully written to: " + filePath); 
       this.getAnalysisCanvas().getCanvas(analTabs[1]).divide(3, 1);
      
       
       for(int i =0; i<3; i++) {
           this.getAnalysisCanvas().getCanvas(analTabs[1]).getPad(i).getAxisZ().setLog(true);
           this.getAnalysisCanvas().getCanvas(analTabs[1]).cd(i);
           this.getAnalysisCanvas().getCanvas(analTabs[1]).draw(TDCvsWire.get(new Coordinate(i)));
           //for(int j=0; j<cuts.get(i).length; j++) { System.out.println(cuts.get(i)[j].getEndY());
            //    this.getAnalysisCanvas().getCanvas(analTabs[1]).draw(cuts.get(i)[j]);
           //}
        }
        double minY = 0;
        double maxY = tHigh4T0Fits[1];
        DataLine l = new DataLine(56, minY, 56, maxY);
        l.setLineWidth(1);
        l.setLineColor(1);
        l.setLineStyle(1);
        l.setArrowSizeOrigin(10);
        DataLine l1 = new DataLine(minY,lowbound[1], maxY, lowbound[1] );
        l1.setLineWidth(2);
        l1.setLineColor(1);
        l1.setLineStyle(4);
        this.getAnalysisCanvas().getCanvas(analTabs[1]).cd(1);
        this.getAnalysisCanvas().getCanvas(analTabs[1]).draw(l);
        this.getAnalysisCanvas().getCanvas(analTabs[1]).draw(l1);
        this.getAnalysisCanvas().getCanvas(analTabs[1]).draw(g, "same");
        
        int[] c = new int[]{0,2};
        for(int k =0; k<c.length; k++) {
            int kk = c[k];
            DataLine l0 = new DataLine(minY,lowbound[kk], maxY, lowbound[kk] );
            l0.setLineWidth(2);
            l0.setLineColor(1);
            l0.setLineStyle(4);
            DataLine h0 = new DataLine(minY,highbound[kk], maxY, highbound[kk] );
            l0.setLineWidth(2);
            l0.setLineColor(1);
            l0.setLineStyle(4);
            this.getAnalysisCanvas().getCanvas(analTabs[1]).cd(kk);
            this.getAnalysisCanvas().getCanvas(analTabs[1]).draw(l0);
            this.getAnalysisCanvas().getCanvas(analTabs[1]).draw(h0);
        }
        this.getCalib().fireTableDataChanged();  
       
       analysisDone=true;
    }
    @Override
    public void constantsEvent(CalibrationConstants cc, int col, int row) {
        String str_sl    = (String) cc.getValueAt(row, 1);
        String str_wir     = (String) cc.getValueAt(row, 2);
        System.out.println(str_sl + " " + str_wir + " " );
        IndexedList<DataGroup> group = this.getDataGroup();
        int region    = Integer.parseInt(str_sl);
        int wire     = Integer.parseInt(str_wir);

       if(group.hasItem(0,region,wire)==true){
           this.Plot(region-1, wire-1);
           
       } else {
           System.out.println(" ERROR: can not find the data group");
       }
   
    }

    private int getWireMax(int i) {
        if(i==1) {
            return nwir;
        } else {
            return 1;
        }
    }
    public void updateOutputArea(List<String[]> dataRows, TDCParamsPanel pm) {
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("%-8s %-8s %-10s %-8s %-8s %-12s\n", 
            "sector", "layer", "component", "MinEdge", "MaxEdge", "LinearCoeff"));
        for (String[] row : dataRows) {
            sb.append(String.format("%-8s %-8s %-10s %-8s %-8s %-12s\n",
                row[0], row[1], row[2], row[3], row[4], row[5]));
        }
        pm.panel.outputArea.setText(sb.toString());
    }
    public void computeBounds(Map<Coordinate, ArrayList<Double>> pars) {
       double Bscale=1;
        if(field!=0) Bscale =field*field;
        
        double R1X1l = 1;
        double R1Y1l = pars.get(new Coordinate(0,0)).get(0);
        double R1X2l = 112;
        double R1Y2l = pars.get(new Coordinate(0,2)).get(0);
        double R1X1u = 1;
        double R1Y1u = pars.get(new Coordinate(1,0)).get(0);
        double R1X2u = 112;
        double R1Y2u = pars.get(new Coordinate(1,2)).get(0);
        
        DataLine r1L  = new DataLine(R1X1l, R1Y1l, R1X2l, R1Y2l);
        DataLine r1U  = new DataLine(R1X1u, R1Y1u, R1X2u, R1Y2u);
        
        cuts.put(0, new DataLine[]{r1L,r1U});
        
        double R3X1l = 1;
        double R3Y1l = pars.get(new Coordinate(4,0)).get(0);
        double R3X2l = 112;
        double R3Y2l = pars.get(new Coordinate(4,2)).get(0);
        double R3X1u = 1;
        double R3Y1u = pars.get(new Coordinate(5,0)).get(0);
        double R3X2u = 112;
        double R3Y2u = pars.get(new Coordinate(5,2)).get(0);
        
        DataLine r3L  = new DataLine(R3X1l, R3Y1l, R3X2l, R3Y2l);
        DataLine r3U  = new DataLine(R3X1u, R3Y1u, R3X2u, R3Y2u);
        
        cuts.put(2, new DataLine[]{r3L,r3U});
        
        double R2X1l = 1;
        double R2Y1l = pars.get(new Coordinate(2,0)).get(0);
        double R2X2l = 112;
        double R2Y2l = pars.get(new Coordinate(2,2)).get(0);
        double R2X1u1 = 1;
        double R2Y1u1 = pars.get(new Coordinate(3,0)).get(0);
        double R2X2u1 = 56;
        double R2Y2u1 = pars.get(new Coordinate(3,1)).get(0);
        double R2X1u2 = 56;
        double R2Y1u2 = pars.get(new Coordinate(3,1)).get(0);
        double R2X2u2 = 112;
        double R2Y2u2 = pars.get(new Coordinate(3,2)).get(0);
        
       
        
        double[] newPars1=getR2Parameters(R2X1u1, R2Y1u1, R2X2u1, R2Y2u1, 56, oldfcn);
        double[] newPars56=getR2Parameters(R2X1u2, R2Y1u2, R2X2u2, R2Y2u2, 112, oldfcn);
         
        double timeCutMax1=newPars1[0];
        double timeCutLC1=newPars1[1];
                
        double timeCutMax56=newPars56[0];
        double timeCutLC56=newPars56[1];
        
        cutParams.put(0, new double[]{R1Y1l, R1Y1u});
        cutParams.put(2, new double[]{R3Y1l, R3Y1u});
        cutParams.put(1, new double[]{R2Y1l, timeCutMax1,timeCutLC1,
                                                R2Y1l, timeCutMax56,timeCutLC56});
        double R2Xi56[] = new double[] {56, 112};
        double R2Yi56[] = new double[2] ;
        for(int xi=0; xi<2; xi++) {
         if(oldfcn) {
                R2Yi56[xi] = timeCutMax56 + timeCutLC56 * (double) (112 - (double)R2Xi56[xi] / (double) 56) * Bscale;
            } else {
                R2Yi56[xi] = timeCutMax56 + timeCutLC56 * (double) (112 - (double)R2Xi56[xi] )/ (double) 56 * Bscale;
            } 
        }
        double R2Xi1[] = new double[] {1, 55};
        double R2Yi1[] = new double[2] ;
        for(int xi=0; xi<2; xi++) {
            if(oldfcn) {
                R2Yi1[xi] = timeCutMax1 + timeCutLC1 * (double) (56 - (double) R2Xi1[xi] / (double) 56) * Bscale ;
            } else {
                R2Yi1[xi] = timeCutMax1 + timeCutLC1 * (double) (56 - (double) R2Xi1[xi] )/ (double) 56 * Bscale ;
            } 
        }
        DataLine r2w56U = new DataLine(R2Xi56[0], R2Yi56[0], R2Xi56[1], R2Yi56[1]);
        DataLine r2w1U  = new DataLine(R2Xi1[0], R2Yi1[0], R2Xi1[1], R2Yi1[1]);
        DataLine r2w56L = new DataLine(56, R2Y1l, 112, R2Y1l);
        DataLine r2w1L = new DataLine(1, R2Y1l, 56, R2Y1l);
        cuts.put(1, new DataLine[]{r2w1L,r2w56L,r2w1U,r2w56U});
            
        for(Integer i : cuts.keySet()) {
            for(int j=0; j<cuts.get(i).length; j++) {
                cuts.get(i)[j].setLineColor(7);
                cuts.get(i)[j].setLineStyle(4);
                cuts.get(i)[j].setLineWidth(3);
            }
        }
        
        for(int i =0; i<3; i++) {
           if(isLog) {
               this.getAnalysisCanvas().getCanvas(analTabs[1]).getPad(i).getAxisZ().setLog(true);
           } else {
               this.getAnalysisCanvas().getCanvas(analTabs[1]).getPad(i).getAxisZ().setLog(false);
           }
           this.getAnalysisCanvas().getCanvas(analTabs[1]).getPad(i).clear();
           this.getAnalysisCanvas().getCanvas(analTabs[1]).cd(i);
           this.getAnalysisCanvas().getCanvas(analTabs[1]).draw(TDCvsWire.get(new Coordinate(i)));
           for(int j=0; j<cuts.get(i).length; j++) { 
                this.getAnalysisCanvas().getCanvas(analTabs[1]).draw(cuts.get(i)[j]);
           }
        }
        
        double minY = 0;
        double maxY = tHigh4T0Fits[1];
        DataLine l = new DataLine(56, minY, 56, maxY);
        l.setLineWidth(1);
        l.setLineColor(1);
        l.setLineStyle(1);
        l.setArrowSizeOrigin(8);
        this.getAnalysisCanvas().getCanvas(analTabs[1]).cd(1);
        this.getAnalysisCanvas().getCanvas(analTabs[1]).draw(l);
        this.getCalib().fireTableDataChanged(); 
    }
    public boolean isLog=true;
    public boolean oldfcn = true;

    private int getLargestGradiant(H1F h, double binWidth) {
    
        return this.getLargestGradiant(h, binWidth, h.getDataSize(0));
    }
    
    private int getLargestGradiant(H1F h, double binWidth, int maxBin) {
        int maxGradiantBin = 0;
        double mgrd=Double.NEGATIVE_INFINITY;
        int min = CutFinderTool.findPlateau(h)[1]; 
        for(int ii =min; ii<maxBin-1; ii++) {
            double bingrd = (h.getDataY(ii)-h.getDataY(ii+1))/binWidth;
            if(bingrd>mgrd) {
                mgrd=bingrd;
                maxGradiantBin=ii+1;
            }
        }
        return maxGradiantBin;
    }
    
    private int getLowEdgeBin(H1F h, double binWidth, int maxBin) {
        List<Integer> bins = new ArrayList<>();
        int maxGradiantBin = 0;
        double mgrd=Double.NEGATIVE_INFINITY;
        
        for(int ii =0; ii<maxBin; ii++) {
            double bingrd = Math.abs(h.getDataY(ii)-h.getDataY(ii+1))/binWidth;
            if(bingrd>mgrd) {
                mgrd=bingrd;
                maxGradiantBin=ii+1;
                //bins.add(ii+1);
            }
        }
        bins.add(maxGradiantBin); System.out.println(maxGradiantBin);
        maxGradiantBin = 0;
        mgrd=Double.NEGATIVE_INFINITY;
        
        for(int ii =0; ii<bins.get(0)-1; ii++) {
            double bingrd = Math.abs(h.getDataY(ii)-h.getDataY(ii+1))/binWidth;
            if(bingrd>mgrd) {
                mgrd=bingrd;
                maxGradiantBin=ii+1;
                //bins.add(ii+1);
            }
        }
        bins.add(maxGradiantBin);
        bins.sort(Comparator.naturalOrder());
        
        return bins.get(0);
    }
    
}
