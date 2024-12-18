/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.dc.calt2d;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import org.clas.detector.clas12calibration.dc.analysis.FitPanel;
import org.clas.detector.clas12calibration.dc.calt2d.FitUtility.MinuitPar;
import org.clas.detector.clas12calibration.dc.t2d.TableLoader;
import org.clas.detector.clas12calibration.viewer.AnalysisMonitor;
import org.clas.detector.clas12calibration.viewer.Driver;
import org.clas.detector.clas12calibration.viewer.T2DViewer;
import static org.clas.detector.clas12calibration.viewer.T2DViewer.ccdb;
import static org.clas.detector.clas12calibration.viewer.T2DViewer.voice;
import org.jlab.detector.calib.utils.CalibrationConstants;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent; 
import org.jlab.detector.calib.utils.ConstantsManager;
import org.jlab.groot.base.TColorPalette;
import org.jlab.groot.data.GraphErrors;
import org.jlab.io.hipo.HipoDataEvent;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.io.hipo.HipoDataSync;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.rec.dc.hit.FittedHit;
import org.jlab.rec.dc.timetodistance.TimeToDistanceEstimator;
import org.jlab.utils.groups.IndexedList;
import org.jlab.utils.system.ClasUtilsFile;
/**
 *
 * @author ziegler
 */
public class T2DCalib extends AnalysisMonitor{

    public static double DeltaTimeCut = 50;
    public HipoDataSync calwriter = null;
    public HipoDataSync writer = null;
    private HipoDataEvent calhipoEvent = null;
    private HipoDataEvent hipoEvent = null;
    private SchemaFactory schemaFactory = new SchemaFactory();
    public FitPanel fp;
    PrintWriter pw = null;
    PrintWriter pw2 = null;
    public PrintWriter pw3 = null;
    File outfile = null;
    private int runNumber;
    public static FcnUtility util = new FcnUtility();
    private int numberprocessedevents;
    private static double betaAve = 1;
    private T2DFitter t2df;
    public static int minSec=6;
    public static int maxSec=7;
    public static boolean vocal= false;
    
    public T2DCalib(String name, ConstantsManager ccdb) throws FileNotFoundException {
        super(name, ccdb);
        this.setAnalysisTabNames("TrackDoca vs T","TrackDoca vs T Graphs","TrackDoca vs T Fit Resi", 
                "CalcDoca vs T","Time Residuals", "Parameters", "Fit Function");
        this.init(false, "v0:vmid:R:tmax:distbeta:delBf:b1:b2:b3:b4");
        
        String dir = ClasUtilsFile.getResourceDir("CLAS12DIR", "etc/bankdefs/hipo4");
        //init BBin Centers
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < this.alphaBins; j++) {
                for (int k = 0; k < this.BBins; k++) {
                    BfieldValuesUpd[i][j][k] = BfieldValues[k];
                }
            }
        }
        //init ABin Centers
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < this.alphaBins; j++) {
                for (int k = 0; k <= this.BBins; k++) {
                    AlphaValuesUpd[i][j][k] = AlphaValues[j];
                }
            }
        }
        
    }
    private Map<Coordinate, H2F> Tvstrkdocas                    = new HashMap<Coordinate, H2F>();
    private Map<Coordinate, H2F> Tvscalcdocas                   = new HashMap<Coordinate, H2F>();
    private Map<Coordinate, H2F> Tresvstrkdocas                 = new HashMap<Coordinate, H2F>();
    private Map<Coordinate, GraphErrors> TvstrkdocasProf        = new HashMap<Coordinate, GraphErrors>();
    private final Map<Coordinate, GraphErrors> TvstrkdocasInit  = new HashMap<Coordinate, GraphErrors>();
    private Map<Coordinate, FitFunction> TvstrkdocasFit         = new HashMap<Coordinate, FitFunction>();
    public static Map<Coordinate, MinuitPar> TvstrkdocasFitPars = new HashMap<Coordinate, MinuitPar>();
    public  Map<Coordinate, FitLine> TvstrkdocasFits            = new HashMap<Coordinate, FitLine>();
    private Map<Coordinate, H1F> timeResi                       = new HashMap<Coordinate, H1F>();
    private Map<Coordinate, H1F> timeResiFromFile               = new HashMap<Coordinate, H1F>();
    private Map<Coordinate, H1F> timeResiNew                    = new HashMap<Coordinate, H1F>();
    private Map<Coordinate, H1F> timeResiB                      = new HashMap<Coordinate, H1F>();
    private Map<Coordinate, H1F> fitResi                        = new HashMap<Coordinate, H1F>();
    private Map<Coordinate, H1F> B                              = new HashMap<Coordinate, H1F>(); //histogram to get B values centroids
    private Map<Coordinate, H1F> A                              = new HashMap<Coordinate, H1F>(); //histogram to get Alpha values centroids
    private Map<Coordinate, H1F> BAlphaBins                     = new HashMap<Coordinate, H1F>();
    private Map<Coordinate, H1F> ParsVsIter                     = new HashMap<Coordinate, H1F>();
    
    Map<Integer, ArrayList<Integer>> segMapTBHits = new HashMap<Integer, ArrayList<Integer>>();
    Map<Integer, ArrayList<HitUtility.CalHit>> segMap2TBHits = new HashMap<Integer, ArrayList<HitUtility.CalHit>>();
    Map<Integer, ArrayList<FittedHit>> segMap = new HashMap<Integer, ArrayList<FittedHit>>();
    Map<Integer, SegmentProperty> segPropMap = new HashMap<Integer, SegmentProperty>();
    
    
    public boolean betaLoaded=false;
    public static double avebeta = 0;
    public static int betacnt =0;
    

    public static double[] BfieldValues = new double[]{0.707106781,1.224744871,1.58113883,1.87082869,2.121320344,2.34520788,2.549509757,2.738612788};   
    public static double[] AlphaValues = new double[]{-40,-36,-32,-28,-24,-20,-16,-12,-8,-4,0,4,8,12,16,20,24,28,32,36,40};
    public static double AlphaBinHalfWidth = 2;
    public static int alphaBins = AlphaValues.length; 
    public static int BBins = BfieldValues.length;
    //update middle of the B bins
    public static double[][][] BfieldValuesUpd = new double[2][alphaBins][BBins];
    public static double[][][] AlphaValuesUpd = new double[6][alphaBins][BBins+1];
    public static int nbinx[] = new int[6];
    public static int nbiny[] = new int[6];
    public static  double docaBinWidth = 0.025;
    public static  double timeBinWidth = 1.0;
    public static  double maxx[] = new double[]{0.85,0.97,1.34,1.38,1.95,2.0};
    public static  double maxy[] = new double[]{400,410,1200,1500,1000,1200};
    
    public static boolean fitBySector = false;
    
    
    @Override
    public void createHistos() {
        // initialize canvas and create histograms
        this.setNumberOfEvents(0);
        
        HistoUtility.createHistos(timeResi, timeResiFromFile, timeResiNew, fitResi, 
                TvstrkdocasFitPars, TvstrkdocasProf, TvstrkdocasInit, Tvstrkdocas, 
                Tresvstrkdocas, Tvscalcdocas, TvstrkdocasFits, A, B, BAlphaBins, timeResiB,
                this.getDataGroup());
       
         
        for(int s = 0; s<6; s++) {
            for (int i0 = 0; i0 < 6; i0++) {   
                for (int j = 0; j < alphaBins; j++) {
                    this.getCalib().addEntry(s+1,i0+1,j+1);
                    //blank out
                    this.getCalib().setDoubleValue((double)999, "v0", s+1, i0+1, j+1);
                    this.getCalib().setDoubleValue((double)999, "vmid", s+1, i0+1, j+1);
                    this.getCalib().setDoubleValue((double)999, "R", s+1, i0+1, j+1);
                    this.getCalib().setDoubleValue((double)999, "tmax", s+1, i0+1, j+1);
                    this.getCalib().setDoubleValue((double)999, "distbeta", s+1, i0+1, j+1);
                    this.getCalib().setDoubleValue((double)999, "delBf", s+1, i0+1, j+1);
                    this.getCalib().setDoubleValue((double)999, "b1", s+1, i0+1, j+1);
                    this.getCalib().setDoubleValue((double)999, "b2", s+1, i0+1, j+1);
                    this.getCalib().setDoubleValue((double)999, "b3", s+1, i0+1, j+1);
                    this.getCalib().setDoubleValue((double)999, "b4", s+1, i0+1, j+1);
                }
            }
        }
        this.getCalib().fireTableDataChanged();
    }
    

    private void resetTable(int secidx, int slidx, int j) {
        this.getCalib().setDoubleValue(999., "v0", secidx+1, slidx+1, j+1);
        this.getCalib().setDoubleValue(999., "vmid", secidx+1, slidx+1, j+1);
        this.getCalib().setDoubleValue(999., "R", secidx+1, slidx+1, j+1);
        this.getCalib().setDoubleValue(999., "tmax", secidx+1, slidx+1, j+1);
        this.getCalib().setDoubleValue(999., "distbeta", secidx+1, slidx+1, j+1);
        this.getCalib().setDoubleValue(999., "delBf", secidx+1, slidx+1, j+1);
        this.getCalib().setDoubleValue(999., "b1", secidx+1, slidx+1, j+1);
        this.getCalib().setDoubleValue(999., "b2", secidx+1, slidx+1, j+1);
        this.getCalib().setDoubleValue(999., "b3", secidx+1, slidx+1, j+1);
        this.getCalib().setDoubleValue(999., "b4", secidx+1, slidx+1, j+1);
    }    
    @Override
    public void plotHistos() {
        String[] Names = {"TrackDoca vs T","TrackDoca vs T Graphs","TrackDoca vs T Fit Resi","CalcDoca vs T","Time Residuals","Parameters",
                            "Fit Function"};
        for(int s = 0; s<3; s++) {
            this.getAnalysisCanvas().getCanvas(Names[s]).setGridX(false);
            this.getAnalysisCanvas().getCanvas(Names[s]).setGridY(false);
            int NumPads = 
            this.getAnalysisCanvas().getCanvas(Names[s]).getCanvasPads().size();
            for (int n = 0; n < NumPads; n++) {
                this.getAnalysisCanvas().getCanvas(Names[s]).getPad(n).getAxisZ().setLog(true);
            }
        }
        this.getAnalysisCanvas().getCanvas(Names[2]).getPad(0).getAxisZ().setLog(true);
        for(int s = 3; s<6; s++) {
            this.getAnalysisCanvas().getCanvas(Names[s]).setGridX(false);
            this.getAnalysisCanvas().getCanvas(Names[s]).setGridY(false);
        }
        this.getAnalysisCanvas().getCanvas(Names[4]).divide(6, 2);
        this.getAnalysisCanvas().getCanvas(Names[5]).divide(6, 6);
        
        this.getAnalysisCanvas().getCanvas(Names[6]).divide(4, 3);
        
        this.getAnalysisCanvas().getCanvas( "TrackDoca vs T Fit Resi").getPad().getAxisY().setRange(-150, 150);
        this.getAnalysisCanvas().getCanvas( "TrackDoca vs T Fit Resi").getPad().setPalette(TColorPalette.PaletteName.kRainBow);
        this.getAnalysisCanvas().getCanvas("TrackDoca vs T").update();
        this.getAnalysisCanvas().getCanvas("TrackDoca vs T Graphs").update();
        this.getAnalysisCanvas().getCanvas( "TrackDoca vs T Fit Resi").update();
        this.getAnalysisCanvas().getCanvas("CalcDoca vs T").update();
        this.getAnalysisCanvas().getCanvas("Time Residuals").update();
        this.getAnalysisCanvas().getCanvas("Parameters").update();
        this.getAnalysisCanvas().getCanvas("Fit Function").update();
    }
    @Override
    public void timerUpdate() {
    }
    
    @Override
    public void analysis() {
        
        try {
            calwriter.close();
            writer.close();
            CalUtility.UpdateAlphaBinCenters(A);
            CalUtility.UpdateBBinCenters(B, BAlphaBins);
            
            for(int s=0; s<7;s++) {
                for (int i = 0; i < 6; i++) {
                    for (int j = 0; j < this.alphaBins; j++) {
                        GraphUtility.filltrkDocavsTGraphs(s,i,j,Tvstrkdocas, TvstrkdocasProf, useBProf);
                        System.out.println("FilledTGraph "+s+" "+i+" "+j);
                    }
                    //runFit(i);
                }
            }
            
            reLoadFitPars();
            System.out.println("FitPars reloaded!");
            //fp.refit();
            //pw.close();
            this.getAnalysisCanvas().getCanvas("Time Residuals").divide(6, 3);
            //
            for(int s=0; s<7;s++) {
                for (int i = 0; i < 6; i++) {
                    this.runInitFit(s,i);
                    
                }
            }
            for (int i = 0; i < 6; i++) {
                this.getAnalysisCanvas().getCanvas("Time Residuals").cd(i);
                PlotUtility.fitTimeResPlot(timeResiFromFile.get(new Coordinate(i)), 
                        this.getAnalysisCanvas().getCanvas("Time Residuals"));
                
            }
                   
            //fp.setGreenFitButton();
            this.plotFits(true);
            
            System.out.println("PLOTS WITH CCDB CONSTANTS DONE");
            this.plotHistos();
            
            for(int s=0; s<6;s++) {
                for (int i = 0; i < 6; i++) {
                    for (int j = 0; j < this.alphaBins; j++) {
                        this.Plot(s, i, j); 
                    }
                }
            }
            eventProcessingDone = true;
            System.out.println("ANALYSIS Done ....");
            
            if(T2DCalib.vocal==true) voice.speak("Event  processing done!");
        } catch (FileNotFoundException ex) {
            Logger.getLogger(T2DCalib.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void plotFits(boolean fitted) throws FileNotFoundException {
        if(fitted==true) {
            DateFormat df = new SimpleDateFormat("MM-dd-yyyy_hh.mm.ss_aa");
            String fileName = "Files/ccdb_run" + this.runNumber + "time_" 
                    + df.format(new Date())+ "iteration_"+this.iterationNum  + ".txt";
           
            String fileName2 = "Files/parameteranderror_run" + this.runNumber + "time_" 
                    + df.format(new Date())+ "iteration_"+this.iterationNum  + ".txt";
           
            
            
            pw = new PrintWriter(fileName);
            pw.printf("#& sector superlayer component v0_a0 v0_a1 v0_a2 vmid_a0 vmid_a1 vmid_a2 tmax_a0 tmax_a1 tmax_a2 distbeta_a0 distbeta_a1 distbeta_a2 delta_bfield_a0 delta_bfield_a1 delta_bfield_a2 b1_a0 b1_a1 b1_a2 b2_a0 b2_a1 b2_a2 b3_a0 b3_a1 b3_a2 b4_a0 b4_a1 b4_a2 c1_a0 c1_a1\n");
            pw2 = new PrintWriter(fileName2);
            pw2.printf("#& sector superlayer component v0 +/-v0 tmax +/-tmax vmid +/-vmid delta_bf +/-delta_bf distbeta +/-distbeta \n");
            
            PlotUtility.plotFits(pw, pw2, NbRunFit,ParsVsIter,TvstrkdocasFitPars,TvstrkdocasFits, TvstrkdocasProf, this.getCalib(),
                    this.getAnalysisCanvas().getCanvas("Parameters"));
            
        }
    }
    int maxIter = 10;
    
    double[][] fixPars = new double[6][10];
    public boolean eventProcessingDone=false;
    
    public int NbRunFit = 0;
    
    public void runInitFit(int s, int i) {
        TvstrkdocasFit.put(new Coordinate(s,i), 
                new FitFunction(s,i, (Map<Coordinate, GraphErrors>) TvstrkdocasProf));
        
    }
    public void initFitParsToFile() throws FileNotFoundException {
        DateFormat df = new SimpleDateFormat("MM-dd-yyyy_hh.mm.ss_aa");
        String fileName3 = "Files/minuit_run" + this.runNumber + "time_" 
                    + df.format(new Date())+ "iteration_"+this.iterationNum  + ".txt";
        pw3 = new PrintWriter(fileName3);
    }
    
    public static double v0limits[][]         = new double[][]{{0.002, 0.009},{0.002, 0.009},{0.002, 0.009},{0.002, 0.009},{0.002, 0.009},{0.002, 0.009}};     //  limits for each superlayer
    public static double vmidlimits[][]       = new double[][]{{0.002, 0.009},{0.002, 0.009},{0.001, 0.009},{0.001, 0.009},{0.001, 0.009},{0.001, 0.009}};     
    public static double Rlimits[][]          = new double[][]{{0.55,0.75},{0.55,0.75},{0.55,0.75},{0.55,0.75},{0.55,0.75},{0.55,0.75}};                            
    public static double distbetalimits[][]   = new double[][]{{0.02, 0.1},{0.02, 0.1},{0.02, 0.1},{0.02, 0.1},{0.02, 0.1},{0.02, 0.1}};                       
    public static double delBflimits[][]      = new double[][]{{0.10, 0.4},{0.10, 0.4},{0.10, 0.4},{0.10, 0.4},{0.10, 0.4},{0.10, 0.4}};                       //  limits for each superlayer
    
    public static double limits[][][] = new double[][][]{v0limits, vmidlimits, Rlimits, distbetalimits, delBflimits};
    boolean useFixedBoundsMethod = true;
    double[][] chi2FitNum = new double[100][6];
    int[] fitNum =new int[6];
    
    public static boolean runParallel = false; //no thread safe; DO NOT USE YET!!!
    
    FitUtility fitUtil = new FitUtility();
    public void runParamScan(boolean fixFit[][][],Map<Coordinate, MinuitPar> TvstrkdocasFitPars) {
        if(T2DCalib.vocal==true) voice.speak("PARAMETER SCAN STARTED");
        for(int s=0; s<7;s++) {
            for (int i = 0; i < 6; i++) {
                this.runInitFit(s,i);

            }
        }
        
        fitUtil.initParsForFit(TvstrkdocasFitPars, fixFit);
        if(runParallel) {
            if(T2DCalib.vocal==true) voice.speak("Starting parallel processing");
            fitUtil.runParamScanParallel(fixFit, TvstrkdocasFitPars, TvstrkdocasFit, TvstrkdocasProf);
        } else {
            fitUtil.runParamScan(fixFit, TvstrkdocasFitPars, TvstrkdocasFit, TvstrkdocasProf);
        }
        fitUtil.releaseParsAfterFit(TvstrkdocasFitPars);
        
        if(T2DCalib.vocal==true) voice.speak("PARAMETER SCAN DONE");
    }
    
    public void runFit(boolean fixFit[][][], Map<Coordinate, MinuitPar> TvstrkdocasFitPars) {
        if(T2DCalib.vocal==true) voice.speak("PARAMETER FIT STARTED");
        for(int s=0; s<7;s++) {
            for (int i = 0; i < 6; i++) {
                this.runInitFit(s,i);
            }
        }
        fitUtil.initParsForFit(TvstrkdocasFitPars, fixFit);
        if(runParallel) {
            if(T2DCalib.vocal==true) voice.speak("Starting parallel processing");
            fitUtil.runFitParallel(fixFit, TvstrkdocasFitPars, TvstrkdocasFit, TvstrkdocasProf);
        } else {
            fitUtil.runFit(fixFit, TvstrkdocasFitPars, TvstrkdocasFit, TvstrkdocasProf);
        }
        fitUtil.releaseParsAfterFit(TvstrkdocasFitPars);
        if(T2DCalib.vocal==true) voice.speak("PARAMETER FIT DONE");
        NbRunFit++;
    }
    
    
    public boolean useBProf = false;
    int counter = 0;
    public static int iterationNum = 0;
    public  HipoDataSource calreader = new HipoDataSource();
    public  HipoDataSource reader = new HipoDataSource();
    private TimeToDistanceEstimator t2d = new TimeToDistanceEstimator();
    public void reCook() {
        if(T2DCalib.vocal==true) voice.speak("Reprocessing the segment fits");
        iterationNum++;
        fp.setRedFitButton();
        //reset histos to refill
        CalUtility.reCook(timeResi, timeResiNew, timeResiB, A, B, BAlphaBins, Tvstrkdocas, Tvscalcdocas, Tresvstrkdocas, 
                calreader, hits, calhits, ParsVsIter, TvstrkdocasFits, TvstrkdocasProf, TvstrkdocasFitPars, useBProf);
        fp.setGreenFitButton();
        if(T2DCalib.vocal==true) voice.speak("Reprocessing done");
    }
    
    public void rePlotResi() {
        PlotUtility.rePlotResi(this.getAnalysisCanvas(), timeResiFromFile, timeResi, timeResiNew);
        
        System.out.println("REPROCESSING DONE!");
        reader.close();
        this.getCalib().fireTableDataChanged();  
    }
      
    
    int count = 0;
    public static int polarity =-1;
    public static double field =1;
    public List<FittedHit> calhits = new ArrayList<>();
    public List<FittedHit> hits = new ArrayList<>();
    Map<Integer, FittedHit> calhitmap = new HashMap<>();
    List<FittedHit> calhitlist = new ArrayList<>();
    Map<Integer, FittedHit> hitmap = new HashMap<>();
    List<FittedHit> hitlist = new ArrayList<>();
    Refit rf = new Refit();
    public static boolean refitSegs = false;
    @Override
    public void processEvent(DataEvent event) {
        
        hitmap.clear();
        calhitmap.clear();
        hitlist.clear();
        calhitlist.clear();
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
         
        //if(count>20000) return;
        if(count==1) {
            if(T2DCalib.vocal==true) voice.speak("Event processing started!");
            HipoReader read = new HipoReader();
            read.open(T2DViewer.theFile);
            schemaFactory = read.getSchemaFactory();
        
            if(schemaFactory.hasSchema("TimeBasedTrkg::TBHits")) {
                System.out.println(" BANK FOUND........");
            } else {
                System.out.println(" BANK NOT FOUND........");
            }
            calwriter = new HipoDataSync(schemaFactory);
            calwriter.setCompressionType(2);
            writer = new HipoDataSync(schemaFactory);
            writer.setCompressionType(2);
            calhipoEvent = (HipoDataEvent) calwriter.createEvent();
            hipoEvent = (HipoDataEvent) writer.createEvent();
            CalUtility.checkFile("TestCalOutPut.hipo");
            CalUtility.checkFile("TestOutPut.hipo");
            calwriter.open("TestCalOutPut.hipo");
            calwriter.writeEvent(calhipoEvent);
            writer.open("TestOutPut.hipo");
            writer.writeEvent(hipoEvent);
            System.out.println("FILES READY......");
            //Constants.getInstance().initialize("DCCAL");
            Driver.init();
            System.out.println("Driver initialized");
            String newVar = String.valueOf(T2DViewer.calVariation.getSelectedItem());
            System.out.println("***************************************************************** VARIATION *"+newVar);
            ccdb.setVariation(newVar);
            TableLoader.t2dc=this;
            TableLoader.FillT0Tables(newRun, newVar);
            
            TableLoader.Fill(T2DViewer.ccdb.getConstants(newRun, "/calibration/dc/time_to_distance/t2d_pressure"),
                    T2DViewer.ccdb.getConstants(newRun, "/calibration/dc/time_to_distance/ref_pressure"),
                    T2DViewer.ccdb.getConstants(newRun, "/hall/weather/pressure"));  
            System.out.println("Filling fit pars");
            try { 
                this.loadFitPars();
            } catch (FileNotFoundException ex) { 
                Logger.getLogger(T2DCalib.class.getName()).log(Level.SEVERE, null, ex);
            }
            polarity = (int)Math.signum(event.getBank("RUN::config").getFloat("torus",0));
            field = event.getBank("RUN::config").getFloat("torus",0);
            runNumber = newRun;
             
            numberprocessedevents = Integer.parseInt(T2DViewer.enternofevents.getText());
            if (numberprocessedevents==-1)
            {
            	System.out.println("All events will be processed!!!!");
            }
            else {
            	System.out.println(numberprocessedevents + " events will be processed!!!!");
            }
        }
        if(!event.hasBank("TimeBasedTrkg::TBHits")) {
            return;
        } 
     
        if (count > numberprocessedevents && numberprocessedevents!=-1) {
        	return;
        }
        
        
        // get segment property     
        DataBank bnkHits = event.getBank("TimeBasedTrkg::TBHits");
        bnkHits=HitUtility.FixDoubleHitsTFlight(event, bnkHits, this, segMap);
        
        HitUtility.getSegProperty(bnkHits,segPropMap,segMapTBHits);
        HistoUtility.fill(event, bnkHits, hitmap, calhitmap, Tvstrkdocas, Tvscalcdocas, TvstrkdocasFits,
                Tresvstrkdocas, timeResiFromFile, A, B, segPropMap);
        
        //hitmap.forEach((k,v) -> hitlist.add(v));
        //calhitmap.forEach((k,v) -> calhitlist.add(v));
        
        List<FittedHit> hits = new ArrayList<>();
        List<FittedHit> calhits = new ArrayList<>();
        hitmap.forEach((k,v) -> hits.add(v));
        calhitmap.forEach((k,v) -> calhits.add(v));
        Set<Integer> hitBankRows=HitUtility.getSegNumLayer(hits, segMap2TBHits);
        for(FittedHit h : hits) {
            if(hitBankRows.contains(h.get_Id())) {
                hitlist.add(h);
            }
        }
        Set<Integer> calhitBankRows=HitUtility.getSegNumLayer(calhits, segMap2TBHits);
        for(FittedHit h : calhits) {
            if(calhitBankRows.contains(h.get_Id())) {
                calhitlist.add(h);
            }
        }
        if(!calhitlist.isEmpty()) {
            calhipoEvent = (HipoDataEvent) calwriter.createEvent();
            calhipoEvent.appendBank(HitUtility.fillTBHitsBank(event, calhitlist));
            calhipoEvent.appendBank(event.getBank("RUN::config"));
            calwriter.writeEvent(calhipoEvent);
        }
        if(!hitlist.isEmpty()) {
            hipoEvent = (HipoDataEvent) writer.createEvent();
            hipoEvent.appendBank(HitUtility.fillTBHitsBank(event, hitlist));
            hipoEvent.appendBank(event.getBank("RUN::config"));
            writer.writeEvent(hipoEvent);
        }
        if(betacnt!=0) {
            setBetaAve(avebeta/(double)betacnt);
            avebeta=0;
            betacnt=0;
        }
        
        calhitlist.clear();
        hitlist.clear();
       
    }
    
    private double[][][] resetPars = new double[7][6][11];
    public static String[] parNames = {"v0", "vmid", "R", "tmax", "distbeta", "delBf", 
        "b1", "b2", "b3", "b4", "dmax"};
    
    public static double[] errs = {0.00001,0.00001,0.01,0.1,0.01,0.001,0.00001,0.00001,0.00001,0.00001,0.0000001};
    //private double[] errs = {0.00001,0.00001,0.001,0.1,0.001,0.001,0.00001,0.00001,0.00001,0.00001,0.0000001};
  //private double[] errs = {0.000001,0.000001,0.00001,0.01,0.0001,0.00001,0.00001,0.00001,0.00001,0.00001,0.0000001};
    
    public void resetPars() {
        fitUtil.resetPars(TvstrkdocasFitPars, ParsVsIter, parNames,resetPars);
        fp.openFitPanel("fit panel", TvstrkdocasFitPars);
        reLoadFitPars();
    }
    
    public void loadFitPars() throws FileNotFoundException {
        
        fitUtil.loadFitPars(Tvstrkdocas, TvstrkdocasFitPars, TvstrkdocasFits, ParsVsIter, resetPars, parNames, this.getCalib());
        
        // Fit panel
        fp = new FitPanel(this);
        fp.openFitPanel("fit panel", TvstrkdocasFitPars);
        System.out.println("FIT PANEL OK");
    }
    private void reLoadFitPars() {
        fitUtil.reLoadFitPars(ParsVsIter,TvstrkdocasFitPars);
    }
    
    public void Plot(int s, int i , int j) {
        PlotUtility.Plot(s, i, j,Tvstrkdocas, Tvscalcdocas, Tresvstrkdocas, TvstrkdocasProf, 
            TvstrkdocasFits,  TvstrkdocasInit, this.getAnalysisCanvas(),  iterationNum);
        
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
           this.Plot(sector-1, layer-1, component-1);
       } else {
           System.out.println(" ERROR: can not find the data group");
       }
       
   
    }

    /**
     * @return the betaAve
     */
    public static double getBetaAve() {
        return betaAve;
    }

    /**
     * @param aBetaAve the betaAve to set
     */
    public static void setBetaAve(double aBetaAve) {
        betaAve = aBetaAve;
    }
    
    public void resetPars(int s, int i, boolean[][] fixFit) {
        // reset
        FitUtility.releasePar(10, TvstrkdocasFitPars.get(new Coordinate(s, i)));
        for (int p = 0; p < 10; p++) {
            if(fixFit[p][i]==true) {
                FitUtility.releasePar(p, TvstrkdocasFitPars.get(new Coordinate(s, i)));
            }
        }
    }

    
}

