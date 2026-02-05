/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.dc.mctuning.analysis.docasmear;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap; 
import java.util.Map;
import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import org.clas.detector.clas12calibration.dc.analysis.CreateDirectory;
import org.clas.detector.clas12calibration.dc.mctuning.viewer.AnalysisMonitor;
import org.freehep.math.minuit.FCNBase;
import org.freehep.math.minuit.FunctionMinimum;
import org.freehep.math.minuit.MnMigrad;
import org.freehep.math.minuit.MnScan;
import org.freehep.math.minuit.MnUserParameters;
import org.jlab.detector.calib.utils.CalibrationConstants;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent; 
import org.jlab.detector.calib.utils.ConstantsManager;
import org.jlab.groot.data.GraphErrors;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.utils.groups.IndexedList;
import org.jlab.utils.system.ClasUtilsFile;
import org.clas.detector.clas12calibration.dc.mctuning.analysis.FitPanel;
import org.clas.detector.clas12calibration.dc.mctuning.viewer.DocaSmearAnalViewer;
import org.clas.detector.clas12calibration.viewer.Driver;
import org.jlab.detector.geant4.v2.DCGeant4Factory;
import org.jlab.groot.fitter.DataFitter;

import org.clas.detector.clas12calibration.dc.t2d.TableLoader;
import org.jlab.rec.dc.Constants;
/**
 *
 * @author ziegler
 */
public class DocaSmearAnal extends AnalysisMonitor{
    private SchemaFactory schemaFactory = new SchemaFactory();
    PrintWriter pw = null;
    File outfile = null;
    private int runNumber;
    private FitPanel fp;
    private String fcString = "fc3";
    public DocaSmearAnal(String name, ConstantsManager ccdb) throws FileNotFoundException {
        super(name, ccdb);
        this.setAnalysisTabNames("Time Resi vs TrackDoca","Time Resi vs TrackDoca Graphs", "Time Smearing vs TrackDoca", "Time Smearing vs TrackDoca Graphs");
        this.init(false, "p0:p1:p2:p3:p4");
        CreateDirectory.create("Files");
        outfile = new File("Files/docasmearConstants.txt");
        pw = new PrintWriter(outfile);
        pw.printf("#& superlayer p0 p1 p2 p3 p4\n");
        
        String dir = ClasUtilsFile.getResourceDir("CLAS12DIR", "etc/bankdefs/hipo4");
        schemaFactory.initFromDirectory(dir);
       
        if(schemaFactory.hasSchema("TimeBasedTrkg::TBHits")) {
            System.out.println(" BANK FOUND........");
        } else {
            System.out.println(" BANK NOT FOUND........");
        }
        
        
    }
    private Map<Coordinate, H2F> timeResVsTrkDoca                      = new HashMap<Coordinate, H2F>();
    private Map<Coordinate, H2F> timeSmearingVsTrkDoca                    = new HashMap<Coordinate, H2F>();
    private Map<Coordinate, GraphErrors> timeResVsTrkDocaProf          = new HashMap<Coordinate, GraphErrors>();
    private Map<Coordinate, GraphErrors> timeSmearingVsTrkDocaProf          = new HashMap<Coordinate, GraphErrors>();
    private Map<Coordinate, FitFunction> timResVsTrkDocaFit            = new HashMap<Coordinate, FitFunction>();
    private Map<Coordinate, MnUserParameters> timeResVsTrkDocaFitPars  = new HashMap<Coordinate, MnUserParameters>();
    public  Map<Coordinate, FitLine> timeResVsTrkDocaFits              = new HashMap<Coordinate, FitLine>();
    private Map<Coordinate, H1F> parsVsBeta                            = new HashMap<Coordinate, H1F>();
    public static Map<Coordinate, H1F> Beta                            = new HashMap<Coordinate, H1F>();
    int nsl = 6;

//    public static double[] betaValues = new double[] {0.65, 0.75, 0.85, 0.95};
//    double betaBinHalfWidth = (betaValues[1]-betaValues[0])*0.5;
//    public static int betaBins = betaValues.length;
//    
     public static double[] betaValues = new double[] {0.5};
    double betaBinHalfWidth = 0.5;
    public static int betaBins = betaValues.length;
    
    @Override
    public void createHistos() {
        // initialize canvas and create histograms
        this.setNumberOfEvents(0);
        
        DataGroup tr = new DataGroup(6,6);
        
        for (int i = 0; i < nsl; i++) {
            for (int j = 0; j < betaBins; j++) {
                DataGroup trkdvst = new DataGroup(1,1);
                DataGroup dvst = new DataGroup(1,1);
                timeResVsTrkDocaFitPars.put(new Coordinate(i,j), new MnUserParameters());
                timeResVsTrkDoca.put(new Coordinate(i,j), new H2F("timeResVsTrkDoca" + (i + 1)*1000+(j+1)+26, "superlayer" + (i + 1)
                            + ", beta ("+(betaValues[j]-betaBinHalfWidth)+", "+(betaValues[j]+betaBinHalfWidth)+")"
                             , 10, 0, 1.0, 40, -1.0, 1.0));
                timeSmearingVsTrkDoca.put(new Coordinate(i,j), new H2F("timeSmearVsTrkDoca" + (i + 1)*1000+(j+1)+26, "superlayer" + (i + 1)
                            + ", beta ("+(betaValues[j]-betaBinHalfWidth)+", "+(betaValues[j]+betaBinHalfWidth)+")"
                             , 10, 0, 1.0, 40, -50.0, 50.0));             
                timeResVsTrkDocaProf.put(new Coordinate(i,j), new GraphErrors());
                timeResVsTrkDocaProf.get(new Coordinate(i,j)).setMarkerColor(1);
                tr.addDataSet(timeResVsTrkDocaProf.get(new Coordinate(i,j)), 0);
                timeResVsTrkDocaFits.put(new Coordinate(i,j), new FitLine());
                timeResVsTrkDocaProf.get(new Coordinate(i,j)).setTitleX("norm. doca");
                timeResVsTrkDocaProf.get(new Coordinate(i,j)).setTitleY("residual (cm)");
                this.getDataGroup().add(tr, 0, i+1, j+1);

                timeSmearingVsTrkDocaProf.put(new Coordinate(i,j), new GraphErrors());
                timeSmearingVsTrkDocaProf.get(new Coordinate(i,j)).setMarkerColor(1);
                tr.addDataSet(timeSmearingVsTrkDocaProf.get(new Coordinate(i,j)), 0);
                timeSmearingVsTrkDocaProf.get(new Coordinate(i,j)).setTitleX("norm. doca");
                timeSmearingVsTrkDocaProf.get(new Coordinate(i,j)).setTitleY("time smearing (ns)");
                this.getDataGroup().add(tr, 0, i+1, j+1);
            }
        }
        
        this.getDataGroup().add(tr, 0,0,0);
        for (int i = 0; i < nsl; i++) {
             for (int j = 0; j < betaBins; j++) {
                Beta.put(new Coordinate(i,j), new H1F("h "+0, "superlayer" + (i + 1), 
                                                             50, betaValues[0]-betaBinHalfWidth, 1.0));
             }
            for(int k = 0; k<5; k++) { //pars
                parsVsBeta.put(new Coordinate(i,k), new H1F("par "+k*1000+i, "superlayer" + (i + 1), 
                        betaValues.length, betaValues[0]-betaBinHalfWidth, betaValues[betaValues.length-1]+betaBinHalfWidth));
                parsVsBeta.get(new Coordinate(i,k)).setTitleX("beta");
                parsVsBeta.get(new Coordinate(i,k)).setTitleY("fit parameter "+k);
            }
        }
        
        
        for (int i = 0; i < nsl; i++) {
            for (int j = 0; j < betaBins; j++) {
                this.getCalib().addEntry(0,i+1,j+1);
                //blank out
                this.getCalib().setDoubleValue((double)999, "p0", 0, i+1, j+1);
                this.getCalib().setDoubleValue((double)999, "p1", 0, i+1, j+1);
                this.getCalib().setDoubleValue((double)999, "p2", 0, i+1, j+1);
                this.getCalib().setDoubleValue((double)999, "p3", 0, i+1, j+1);
                this.getCalib().setDoubleValue((double)999, "p4", 0, i+1, j+1);
                this.getCalib().setDoubleValue((double)999, "p5", 0, i+1, j+1);
            }
        }
        
        this.getCalib().fireTableDataChanged();
    }
    private void updateTable(int i, int j) {
        this.getCalib().setDoubleValue(timeResVsTrkDocaFitPars.get(new Coordinate(i,j)).value(0), "p0", 0, i+1, j+1);
        this.getCalib().setDoubleValue(timeResVsTrkDocaFitPars.get(new Coordinate(i,j)).value(1), "p1", 0, i+1, j+1);
        this.getCalib().setDoubleValue(timeResVsTrkDocaFitPars.get(new Coordinate(i,j)).value(2), "p2", 0, i+1, j+1);
        this.getCalib().setDoubleValue(timeResVsTrkDocaFitPars.get(new Coordinate(i,j)).value(3), "p3", 0, i+1, j+1);
        this.getCalib().setDoubleValue(timeResVsTrkDocaFitPars.get(new Coordinate(i,j)).value(4), "p4", 0, i+1, j+1);
    }    
    @Override
    public void plotHistos() {
        String[] Names = {"Time Resi vs TrackDoca","Time Resi vs TrackDoca Graphs", "Time Smearing vs TrackDoca Graphs"};
        
        for(int s = 0; s<Names.length; s++) {
            this.getAnalysisCanvas().getCanvas(Names[s]).setGridX(false);
            this.getAnalysisCanvas().getCanvas(Names[s]).setGridY(false);
        
            int NumPads = 
            this.getAnalysisCanvas().getCanvas(Names[s]).getCanvasPads().size();
            for (int n = 0; n < NumPads; n++) {
                this.getAnalysisCanvas().getCanvas(Names[s]).getPad(n).getAxisZ().setLog(true);
            }
        }
        
        
        
        this.getAnalysisCanvas().getCanvas("Time Resi vs TrackDoca").update();
        this.getAnalysisCanvas().getCanvas("Time Resi vs TrackDoca Graphs").update();
        this.getAnalysisCanvas().getCanvas("Time Smearing vs TrackDoca").update();
        this.getAnalysisCanvas().getCanvas("Time Smearing vs TrackDoca Graphs").update();
    }
    
    @Override
    public void timerUpdate() {
    }
    
    @Override
    public void analysis() {
        loadFitPars();
        
        for (int i = 0; i < this.nsl; i++) {
            for (int j = 0; j < this.betaBins; j++) {
                this.filltrkDocavsTGraphs(i,j);
                this.runFit(i, j);
            }
        }
        plotFits() ;
        for (int i = 0; i < this.nsl; i++) {
            for (int j = 0; j < this.betaBins; j++) {
                this.Plot(i,j);
            }
        }
//        int ik =0;
//        for (int i = 0; i < this.nsl; i++) {
//            for (int k = 0; k < 5; k++) {
//                this.getAnalysisCanvas().getCanvas("Beta Dependence").cd(ik);
//                GraphErrors ge = parsVsBeta.get(new Coordinate(i,k)).getGraph();
//                ge.addPoint(0, 0, 0, 1);
//                ge.setTitleX("beta");
//                ge.setTitleY("fit parameter "+k);
//                this.getAnalysisCanvas().getCanvas("Beta Dependence").draw(
//                        ge, "E");
//                ik++;
//            }
//        }
    }
    public void plotFits() {
        pw.close();
        File file2 = new File("");
        file2 = outfile;
        DateFormat df = new SimpleDateFormat("MM-dd-yyyy_hh.mm.ss_aa");
        String fileName = "Files/docasmear_run" + this.runNumber + "time_" 
                + df.format(new Date())+ ".txt";
        file2.renameTo(new File(fileName));
        int ij =0;
        int ip =0;
        for (int i = 0; i < this.nsl; i++) {

            for (int j = 0; j < this.betaBins; j++) {

                if(timeResVsTrkDocaProf.get(new Coordinate(i, j)).getVectorX().size()>0) {
                    this.updateTable(i,j);
                    timeResVsTrkDocaFits.put(new Coordinate(i,j), new FitLine("f"+""+i+""+j+"0", i, j, this.fcString,
                    timeResVsTrkDocaFitPars.get(new Coordinate(i,j))));
                    timeResVsTrkDocaFits.get(new Coordinate(i, j)).setLineStyle(4);
                    timeResVsTrkDocaFits.get(new Coordinate(i, j)).setLineWidth(5);
                    timeResVsTrkDocaFits.get(new Coordinate(i, j)).setLineColor(8);
                }

                ij++;
            }
            this.getCalib().fireTableDataChanged();  
        }
    }
    private int maxIter = 2;
    double[][] fixPars = new double[6][10];
    
    private MnScan  scanner = null;
    private MnMigrad migrad = null;
    
    public int NbRunFit = 0;
    public void runFit(int i, int j) {
        // i = superlayer - 1;
        System.out.println(" **************** ");
        System.out.println(" RUNNING THE FITS "+fcString);
        System.out.println(" **************** ");
        timResVsTrkDocaFit.put(new Coordinate(i, j), 
                new FitFunction(i, j, fcString, (Map<Coordinate, GraphErrors>) timeResVsTrkDocaProf));
        
        scanner = new MnScan((FCNBase) timResVsTrkDocaFit.get(new Coordinate(i,j)), 
        timeResVsTrkDocaFitPars.get(new Coordinate(i,j)),2);
        if(this.fcString.equalsIgnoreCase("fc2")==true) {
//            timeResVsTrkDocaFitPars.get(new Coordinate(i,j)).setLimits(0, 6.9-0.0001,6.9+0.0001);
//            timeResVsTrkDocaFitPars.get(new Coordinate(i,j)).setLimits(1, 0-1, 0.01+1);
//            timeResVsTrkDocaFitPars.get(new Coordinate(i,j)).setLimits(2, -0.001-1, 0.01+1);
//            timeResVsTrkDocaFitPars.get(new Coordinate(i,j)).setLimits(3, -0.0001-1, 0.006+1);
//            timeResVsTrkDocaFitPars.get(new Coordinate(i,j)).setLimits(4, 0.6-1, 1.0);
        } else {
        //    timeResVsTrkDocaFitPars.get(new Coordinate(i,j)).setLimits(0, 0.9, 1.1);
        //    timeResVsTrkDocaFitPars.get(new Coordinate(i,j)).setLimits(1, 0, 0.01);
       //     timeResVsTrkDocaFitPars.get(new Coordinate(i,j)).setLimits(2, 0.03, 0.05);
        //    timeResVsTrkDocaFitPars.get(new Coordinate(i,j)).setLimits(3, -0.002, 0.00);
        //    timeResVsTrkDocaFitPars.get(new Coordinate(i,j)).setLimits(4, 0.7, 0.9);
        }
        if(this.fcString.equalsIgnoreCase("fc3")==true){
            double a0 = this.getA0(timeResVsTrkDocaProf.get(new Coordinate(i,j)));
            timeResVsTrkDocaFitPars.get(new Coordinate(i,j)).setValue(0, a0);
            timeResVsTrkDocaFitPars.get(new Coordinate(i,j)).fix(0);
            //timeResVsTrkDocaFitPars.get(new Coordinate(i,j)).fix(4);
        }
	//for (int p = 0; p < 1; p++) {
        //    scanner.fix(p); 
        //}
        //for(int p = 1; p<5; p++) {
        //    scanner.fix(p);
        //}
        System.out.println(" Ready to minimize..... ");
        FunctionMinimum scanmin = scanner.minimize();
        for(int pi = 0; pi<5; pi++) 
                System.out.println("scan par["+pi+"]="+timeResVsTrkDocaFitPars.get(new Coordinate(i,j)).value(pi));
           
        //if(scanmin.isValid())
            timeResVsTrkDocaFitPars.put(new Coordinate(i,j),scanmin.userParameters());
        
        migrad = new MnMigrad((FCNBase) timResVsTrkDocaFit.get(new Coordinate(i,j)), 
                timeResVsTrkDocaFitPars.get(new Coordinate(i,j)),1);
        migrad.setCheckAnalyticalDerivatives(true);
        
        FunctionMinimum min ;
        
        for(int it = 0; it<maxIter; it++) {
            
            min = migrad.minimize();
            System.err.println("****************************************************");
            System.err.println("*   FIT RESULTS  FOR SUPERLAYER  "+(i+1)+" at iteration "+(it+1)+"  *");
            System.err.println("****************************************************");  
            for(int pi = 0; pi<5; pi++) 
                System.out.println("par["+pi+"]="+(timeResVsTrkDocaFitPars.get(new Coordinate(i,j)).value(pi)-min.userParameters().value(pi)));
            
            //if(min.isValid()) {
                timeResVsTrkDocaFitPars.put(new Coordinate(i,j),min.userParameters());  
            //}
            System.err.println(min);
                 
        }
        
//        for(int isec = 0; isec < 6; isec++) {
//           
//            pw.printf("%d\t %d\t %d\t %.6f\t %d\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %d\t %.6f\t %.6f\t %d\n",
//                (isec+1), (i+1), 
//                timeResVsTrkDocaFitPars.get(new Coordinate(i)).value(0),
//                timeResVsTrkDocaFitPars.get(new Coordinate(i)).value(1),
//                timeResVsTrkDocaFitPars.get(new Coordinate(i)).value(2),
//                timeResVsTrkDocaFitPars.get(new Coordinate(i)).value(3),
//                timeResVsTrkDocaFitPars.get(new Coordinate(i)).value(4),
//                0);
//        }
        //for(int k = 0; k<1; k++) {
        //    timeResVsTrkDocaFitPars.get(new Coordinate(i,j)).release(k);
        //}
        if(this.fcString.equalsIgnoreCase("fc2")==true) {
            for(int k = 0; k<5; k++) {
                timeResVsTrkDocaFitPars.get(new Coordinate(i,j)).removeLimits(k);
            }
        }
        if(this.fcString.equalsIgnoreCase("fc3")==true){
            timeResVsTrkDocaFitPars.get(new Coordinate(i,j)).release(0);
            //timeResVsTrkDocaFitPars.get(new Coordinate(i,j)).release(4);
        }
        for(int k = 0; k<5; k++) {
            parsVsBeta.get(new Coordinate(i,k)).setBinContent(j, timeResVsTrkDocaFitPars.get(new Coordinate(i,j)).value(k));
            parsVsBeta.get(new Coordinate(i,k)).setBinError(j, timeResVsTrkDocaFitPars.get(new Coordinate(i,j)).error(k));
            //if(k>0)
            //    timeResVsTrkDocaFitPars.get(new Coordinate(i,j)).release(k);
        } 
        
    }
    
    int counter = 0;
    private int iterationNum = 1;
   
    
    private int getBetaBin(int sli, double b) {
        
        int v = 0;
        for(int i = 0; i<betaValues.length; i++) {
            if(Math.abs(b-betaValues[i])<this.betaBinHalfWidth)
                v = i;
        } 
        if(b==1.0) v = betaValues.length-1;
        Beta.get(new Coordinate(sli,v)).fill(b);
        return v;
    }

    
    private int MINENTRIES = 10;
    F1D f1 = new F1D("f1","[amp]*gaus(x,[mean],[sigma])", -0.5, 0.5);
    F1D f2 = new F1D("f2","[amp]*gaus(x,[mean],[sigma])", -40.0, 40.0);
    // F1D f1 = new F1D("f1","[amp]*gaus(x,[mean],[sigma])+[p0]", -0.5, 0.5);
    int[] BinMin = {5,5,5,5,5,5};
    int[] BinMax = {32,32,32,32,32,32};
    private void filltrkDocavsTGraphs(int i, int j) {
        
        if(timeResVsTrkDocaProf.get(new Coordinate(i, j))!=null) {
            
            timeResVsTrkDocaProf.get(new Coordinate(i, j)).reset();
            H2F h2 = timeResVsTrkDoca.get(new Coordinate(i, j));
            ArrayList<H1F> hslice = h2.getSlicesX();
            
            for(int si=0; si<hslice.size(); si++) {
                double amp   = hslice.get(si).getBinContent(hslice.get(si).getMaximumBin());
                
                if(amp<this.MINENTRIES 
                    //    ||  si<BinMin[i] || si>BinMax[i] 
                        ) {
                    
                } else {
                    double x = h2.getXAxis().getBinCenter(si);
                    double y = hslice.get(si).getMean();
                    double sigma = hslice.get(si).getRMS();
                    f1.setParameter(0, amp);
                    f1.setParameter(1, y);
                    f1.setParameter(2, sigma);
                    //f1.setParameter(3, 0);
                    DataFitter.fit(f1, hslice.get(si), "Q"); //No options uses error for sigma 
                    
                    if(f1.parameter(1).error()>0 && Math.abs(f1.parameter(2).error())<0.1) {
                        timeResVsTrkDocaProf.get(new Coordinate(i, j)).
                                addPoint(x, Math.abs(f1.getParameter(2)), 0, Math.abs(f1.parameter(2).error()));
                    } 
                    timeResVsTrkDocaProf.get(new Coordinate(i, j)).setError(0, 0, 0.01);
                }
            }
            
            
             //timeResVsTrkDocaProf.get(new Coordinate(i, j)).
              //                  addPoint(0, 0.07, 0, 0.1);
        }
        
        if(timeSmearingVsTrkDocaProf.get(new Coordinate(i, j))!=null) {
            
            timeSmearingVsTrkDocaProf.get(new Coordinate(i, j)).reset();
            H2F h2 = timeSmearingVsTrkDoca.get(new Coordinate(i, j));
            ArrayList<H1F> hslice = h2.getSlicesX();
            
            for(int si=0; si<hslice.size(); si++) {
                double amp   = hslice.get(si).getBinContent(hslice.get(si).getMaximumBin());
                double x = h2.getXAxis().getBinCenter(si);
                double y = hslice.get(si).getMean();
                double sigma = hslice.get(si).getRMS();
                int N=0;
                for(int k =0; k<10; k++) {
                    N+=hslice.get(si).getBinContent(k);
                }
                f2.setParameter(0, amp);
                f2.setParameter(1, y);
                f2.setParameter(2, sigma);
                //f1.setParameter(3, 0);
                DataFitter.fit(f2, hslice.get(si), "Q"); //No options uses error for sigma 

                if(f2.parameter(1).error()>0 ) {
                    timeSmearingVsTrkDocaProf.get(new Coordinate(i, j)).
                            addPoint(x, Math.abs(f2.getParameter(2)), 0, Math.abs(f2.parameter(2).error()));
                } else {
                    timeSmearingVsTrkDocaProf.get(new Coordinate(i, j)).
                                addPoint(x, sigma, 0, sigma/(Math.sqrt(2*N)));
                }
                
                timeSmearingVsTrkDocaProf.get(new Coordinate(i, j)).setError(0, 0, 10.0);
                     
            }
        }
        
    }
    public static double[][] v0 = new double[6][6];
    int count = 0;
    public static int polarity =-1;
    public DCGeant4Factory DcDetector;
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
            //Constants.getInstance().initialize("DCCAL");
            Driver.init();
            polarity = (int)Math.signum(event.getBank("RUN::config").getFloat("torus",0));
            runNumber = newRun;
            if(count==1) {
                Driver.init();
                TableLoader.FillT0Tables(newRun, "default");
                TableLoader.Fill(DocaSmearAnalViewer.ccdb.getConstants(newRun, "/calibration/dc/v2/t2d_pressure"),
                        DocaSmearAnalViewer.ccdb.getConstants(newRun, "/calibration/dc/v2/ref_pressure"),
                        DocaSmearAnalViewer.ccdb.getConstants(newRun, "/hall/weather/pressure"));  
                runNumber = newRun;
            }
            
            for(int s = 0; s<6; s++ ){ // loop over sectors
                for(int r = 0; r<6; r++ ){ //loop over slys
                    // Fill constants
                    v0[s][r] = org.jlab.rec.dc.timetodistance.TableLoader.v0[s][r];
                }
            }
        }
        if(!event.hasBank("TimeBasedTrkg::TBHits")) {
            return;
        } 
        // get segment property
        
        DataBank bnkHits = event.getBank("TimeBasedTrkg::TBHits");
        
        for (int i = 0; i < bnkHits.rows(); i++) {
            int sector = bnkHits.getInt("sector", i);    
            int superlayer = bnkHits.getInt("superlayer", i);
            int layer = bnkHits.getInt("layer", i);
            double beta = bnkHits.getFloat("beta", i);    
            double doca = bnkHits.getFloat("trkDoca", i);   
            double B = bnkHits.getFloat("B", i);
            double alpha = bnkHits.getFloat("Alpha", i);
            if (bnkHits.getByte("trkID", i) >0 
                    && beta > betaValues[0]-this.betaBinHalfWidth && beta<1.3
                    && bnkHits.getFloat("TFlight", i)>0 
                    && Math.abs(bnkHits.getFloat("fitResidual", i))<0.075)
            {
                if(beta>1)
                    beta=1;
                int betaBin = this.getBetaBin(superlayer-1, beta);
                //System.out.println(beta+" --> beta bin : "+betaBin);
                double x1 = DcDetector.getWireMidpoint(sector-1, superlayer-1, layer-1, 1).x;
                double x0 = DcDetector.getWireMidpoint(sector-1, superlayer-1, layer-1, 0).x;
                double deltax = Math.abs(x1-x0)/2;
                timeResVsTrkDoca.get(new Coordinate(superlayer - 1, betaBin))
                                .fill(Math.abs(bnkHits.getFloat("trkDoca", i))/deltax, 
                                        bnkHits.getFloat("timeResidual", i));
                double v0 =org.jlab.rec.dc.timetodistance.TableLoader.v0[sector-1][superlayer-1];
                double vmid = org.jlab.rec.dc.timetodistance.TableLoader.vmid[sector-1][superlayer-1];
                double R = org.jlab.rec.dc.timetodistance.TableLoader.FracDmaxAtMinVel[sector-1][superlayer-1];
                double tmax = org.jlab.rec.dc.timetodistance.TableLoader.Tmax[sector-1][superlayer-1];
                double distbeta = org.jlab.rec.dc.timetodistance.TableLoader.distbeta[sector-1][i];
                double delBf = org.jlab.rec.dc.timetodistance.TableLoader.delta_bfield_coefficient[sector-1][superlayer-1];
                double bB1 = org.jlab.rec.dc.timetodistance.TableLoader.b1[sector-1][superlayer-1];
                double bB2 = org.jlab.rec.dc.timetodistance.TableLoader.b2[sector-1][superlayer-1];
                double bB3 = org.jlab.rec.dc.timetodistance.TableLoader.b3[sector-1][superlayer-1];
                double bB4 = org.jlab.rec.dc.timetodistance.TableLoader.b4[sector-1][superlayer-1];
                double dmax = 2.*Constants.getInstance().wpdist[superlayer-1];
                //T2DDerivative(double x, double alpha, double bfield, double v_0, double vm, double R, 
                //double tmax, double dmax, double delBf, double Bb1, double Bb2, double Bb3, double Bb4, double beta, double distbeta, int superlayer)
                double derivT2D = org.clas.detector.clas12calibration.dc.calt2d.FitFunction.T2DDerivative(doca, alpha, B, v0, vmid, R, tmax, dmax, delBf, bB1, bB2, bB3, bB4, beta, distbeta, superlayer);
                timeSmearingVsTrkDoca.get(new Coordinate(superlayer - 1, betaBin))
                                .fill(Math.abs(bnkHits.getFloat("trkDoca", i))/deltax, 
                                        bnkHits.getFloat("timeResidual", i)*derivT2D);
            }
        }
    }
    
    private String[] parNames = {"p0", "p1", "p2", "p3", "p4"};
    private double[] errs = {0.001,0.001,0.01,0.01,0.01};
    
    public void loadFitPars() {
        for (int i = 0; i < this.nsl; i++) {
            for (int j = 0; j < betaValues.length; j++) {
                double[] pars = new double[parNames.length];
                //a1 = 0.05 cm
                //a2 = 0.02 cm
                //a3 = 0.012 cm
                //a4 = 0.1 cm
                //scale_factor is nominally 1
                //V0 (cm/ns) is gotten from CCDB
                if(this.fcString.equalsIgnoreCase("fc3")==true){
                    pars[0] = 3.511850e-02;
                    pars[1] = -2.528271e-02;
                    pars[2] = 1.873393e-03;
                    pars[3] =  3.005418e-02;
                    pars[4] = 0.0;
                } else {
                    pars[0] = 3.5;
                    pars[1] = 0.002041;
                    pars[2] = -0.000382;
                    pars[3] = 0.028848;
                    pars[4] = 0.861841;
                }
                timeResVsTrkDocaFitPars.put(new Coordinate(i,j), new MnUserParameters());
                for(int p = 0; p < parNames.length; p++) {
                    timeResVsTrkDocaFitPars.get(new Coordinate(i,j)).add(parNames[p], pars[p], errs[p]);
                }
                for(int pi = 0; pi<5; pi++) 
                    System.out.println("loaded par["+pi+"]="+(timeResVsTrkDocaFitPars.get(new Coordinate(i,j)).value(pi)));
           
            }   
        }
        // Fit panel
        //fp = new FitPanel();
        //fp.openFitPanel("fit panel", timeResVsTrkDocaFitPars);
    }
    
    public void Plot(int i , int j) {
        
        this.getAnalysisCanvas().getCanvas("Time Resi vs TrackDoca").cd(0);
        this.getAnalysisCanvas().getCanvas("Time Resi vs TrackDoca").draw(timeResVsTrkDoca.get(new Coordinate(i, j)));
        this.getAnalysisCanvas().getCanvas("Time Smearing vs TrackDoca").cd(0);
        this.getAnalysisCanvas().getCanvas("Time Smearing vs TrackDoca").draw(timeSmearingVsTrkDoca.get(new Coordinate(i, j)));
        this.getAnalysisCanvas().getCanvas("Time Resi vs TrackDoca Graphs").cd(0);
        //this.getAnalysisCanvas().getCanvas("Time Resi vs TrackDoca Graphs").draw(timeResVsTrkDoca.get(new Coordinate(i, j)));
        if(timeResVsTrkDocaProf.get(new Coordinate(i, j)).getVectorX().size()>0){
                    this.getAnalysisCanvas().getCanvas("Time Resi vs TrackDoca Graphs").
                            draw(timeResVsTrkDocaProf.get(new Coordinate(i, j)), "E");
                    this.getAnalysisCanvas().getCanvas("Time Resi vs TrackDoca Graphs").
                            draw(timeResVsTrkDocaFits.get(new Coordinate(i, j)), "same");
        }
        
        this.getAnalysisCanvas().getCanvas("Time Smearing vs TrackDoca Graphs").cd(0);
        //this.getAnalysisCanvas().getCanvas("Time Resi vs TrackDoca Graphs").draw(timeResVsTrkDoca.get(new Coordinate(i, j)));
        if(timeSmearingVsTrkDocaProf.get(new Coordinate(i, j)).getVectorX().size()>0){
                    this.getAnalysisCanvas().getCanvas("Time Smearing vs TrackDoca Graphs").
                            draw(timeSmearingVsTrkDocaProf.get(new Coordinate(i, j)), "E");
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
           this.Plot(layer-1, component-1);
       } else {
           System.out.println(" ERROR: can not find the data group");
       }
       
   
    }

    private double getA0(GraphErrors g) {
        double x1 = g.getDataX(1);
        double x2 = g.getDataX(2);
        double y1 = g.getDataY(1);
        double y2 = g.getDataY(2);
        double s = (y2-y1)/(x2-x1);
        return y1-s*x1;
    }
}

