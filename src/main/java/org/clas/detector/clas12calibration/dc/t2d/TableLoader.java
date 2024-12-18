package org.clas.detector.clas12calibration.dc.t2d;

import java.util.HashMap;
import java.util.Map; 
import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import org.clas.detector.clas12calibration.dc.calt2d.FitLine;
import org.clas.detector.clas12calibration.dc.calt2d.T2DCalib;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.BBins;
import org.clas.detector.clas12calibration.dc.calt2d.FcnUtility;
import org.clas.detector.clas12calibration.dc.calt2d.FitUtility.MinuitPar;
import org.clas.detector.clas12calibration.viewer.T2DViewer;
import org.freehep.math.minuit.MnUserParameters;
import org.jlab.detector.calib.utils.DatabaseConstantProvider;
import org.jlab.rec.dc.Constants;
import org.jlab.rec.dc.timetodistance.TimeToDistanceEstimator;
import org.jlab.utils.groups.IndexedTable;


public class TableLoader {

    public TableLoader() {
            // TODO Auto-generated constructor stub
    }
    public static T2DCalib t2dc;
    static final protected int nBinsT=2000;
    //public static double[][][][][] DISTFROMTIME = new double[6][6][6][6][850]; // sector slyr alpha Bfield time bins
    
    static boolean T2DLOADED = false;
    static boolean T0LOADED = false;
       
    private static double[][][][] T0 ;
    private static double[][][][] T0ERR ;
    public static synchronized void FillT0Tables(int run, String variation) {
        if (T0LOADED) return;
        System.out.println(" T0 TABLE FILLED..... for Run "+run+" with VARIATION "+variation);
        DatabaseConstantProvider dbprovider = new DatabaseConstantProvider(run, variation);
        dbprovider.loadTable("/calibration/dc/time_corrections/T0Corrections");
        //disconnect from database. Important to do this after loading tables.
        dbprovider.disconnect();
        // T0-subtraction
        double[][][][] T0 ;
        double[][][][] T0ERR ;
        //T0s
        T0 = new double[6][6][7][6]; //nSec*nSL*nSlots*nCables
        T0ERR = new double[6][6][7][6]; //nSec*nSL*nSlots*nCables
        for (int i = 0; i < dbprovider.length("/calibration/dc/time_corrections/T0Corrections/Sector"); i++) {
            int iSec = dbprovider.getInteger("/calibration/dc/time_corrections/T0Corrections/Sector", i);
            int iSly = dbprovider.getInteger("/calibration/dc/time_corrections/T0Corrections/Superlayer", i);
            int iSlot = dbprovider.getInteger("/calibration/dc/time_corrections/T0Corrections/Slot", i);
            int iCab = dbprovider.getInteger("/calibration/dc/time_corrections/T0Corrections/Cable", i);
            double t0 = dbprovider.getDouble("/calibration/dc/time_corrections/T0Corrections/T0Correction", i);
            double t0Error = dbprovider.getDouble("/calibration/dc/time_corrections/T0Corrections/T0Error", i);
            T0[iSec - 1][iSly - 1][iSlot - 1][iCab - 1] = t0; 
            T0ERR[iSec - 1][iSly - 1][iSlot - 1][iCab - 1] = t0Error;
            setT0(T0);
            setT0Err(T0ERR);
            //System.out.println("T0 = "+t0);
        }
        T0LOADED = true;
    }

    /**
     * @return the T0
     */
    public static double[][][][] getT0() {
        return T0;
    }

    /**
     * @param aT0 the T0 to set
     */
    public static void setT0(double[][][][] aT0) {
        T0 = aT0;
    }

    /**
     * @return the T0ERR
     */
    public static double[][][][] getT0Err() {
        return T0ERR;
    }

    /**
     * @param aT0ERR the T0ERR to set
     */
    public static void setT0Err(double[][][][] aT0ERR) {
        T0ERR = aT0ERR;
    }
    
   

    static Map<Coordinate, FitLine> TvsDB = new HashMap<Coordinate, FitLine>();
    static Map<Coordinate, FitLine> TvsDBr = new HashMap<Coordinate, FitLine>();
    //FitLine(String name, int i, int j, int k, MnUserParameters pars)
    static String[] parNames = {"v0", "vmid", "R", "tmax", "distbeta", "delBf", 
        "b1", "b2", "b3", "b4", "dmax"};
    static double[] errs = {0.001,0.001,0.01,1.0,0.01,0.001,0.001,0.001,0.001,0.001,0.00001};
    
    private static synchronized void fitsinit(){
        for (int s = 0; s < 6; s++) {
            for (int i = 0; i < 6; i++) {
                double[] pars = new double[11];
                pars[0] = org.jlab.rec.dc.timetodistance.TableLoader.v0[s][i];
                pars[1] = org.jlab.rec.dc.timetodistance.TableLoader.vmid[s][i];
                pars[2] = org.jlab.rec.dc.timetodistance.TableLoader.FracDmaxAtMinVel[s][i];
                pars[3] = org.jlab.rec.dc.timetodistance.TableLoader.Tmax[s][i];
                pars[4] = org.jlab.rec.dc.timetodistance.TableLoader.distbeta[s][i];
                pars[5] = org.jlab.rec.dc.timetodistance.TableLoader.delta_bfield_coefficient[s][i];
                pars[6] = org.jlab.rec.dc.timetodistance.TableLoader.b1[s][i];
                pars[7] = org.jlab.rec.dc.timetodistance.TableLoader.b2[s][i];
                pars[8] = org.jlab.rec.dc.timetodistance.TableLoader.b3[s][i];
                pars[9] = org.jlab.rec.dc.timetodistance.TableLoader.b4[s][i];
                pars[10] = 2.*Constants.getInstance().wpdist[i];//fix dmax
                MinuitPar mnp = new MinuitPar();

                for(int p = 0; p < 10; p++) {
                    mnp.add(parNames[p], pars[p], errs[p]);
                }
                mnp.add(parNames[10], pars[10], errs[10]);
                for (int j = 0; j < T2DCalib.alphaBins; j++) {
                    for (int k = 0; k < BBins+1; k++) {
                            TvsDB.put(new Coordinate(s, i,j,k), new FitLine("f_"+"s"+s+"i"+i+"j"+j+"k"+k, i, j, k, 
                                    mnp));
                            TvsDB.get(new Coordinate(s, i, j, k)).setLineWidth(2);
                            TvsDB.get(new Coordinate(s, i, j, k)).setLineColor(k+1);
                            TvsDB.get(new Coordinate(s, i, j, k)).setRange(0, pars[10]);   
                            TvsDB.get(new Coordinate(s, i, j, k)).useMidBfieldBin=true;
                            TvsDB.get(new Coordinate(s, i, j, k)).useMidAlphaBin=true;
                            TvsDBr.put(new Coordinate(s, i,j,k), new FitLine("f_"+"s"+s+"i"+i+"j"+j+"k"+k, i, j, k, 
                                   mnp));
                            TvsDBr.get(new Coordinate(s, i, j, k)).setLineWidth(2);
                            TvsDBr.get(new Coordinate(s, i, j, k)).setLineColor(k+1);
                            TvsDBr.get(new Coordinate(s, i, j, k)).setRange(0, pars[10]); 
                            TvsDBr.get(new Coordinate(s, i, j, k)).useMidBfieldBin=true;
                            TvsDBr.get(new Coordinate(s, i, j, k)).useMidAlphaBin=true;
                    }
                }
            }
        }
        System.out.println("Fit Function reloaded...");
    }
    private static synchronized void reset(){
        TvsDBr.clear();
        for (int s = 0; s < 6; s++) {
            for (int i = 0; i < 6; i++) {
                double[] pars = new double[11];
                pars[0] = org.jlab.rec.dc.timetodistance.TableLoader.v0[s][i];
                pars[1] = org.jlab.rec.dc.timetodistance.TableLoader.vmid[s][i];
                pars[2] = org.jlab.rec.dc.timetodistance.TableLoader.FracDmaxAtMinVel[s][i];
                pars[3] = org.jlab.rec.dc.timetodistance.TableLoader.Tmax[s][i];
                pars[4] = org.jlab.rec.dc.timetodistance.TableLoader.distbeta[s][i];
                pars[5] = org.jlab.rec.dc.timetodistance.TableLoader.delta_bfield_coefficient[s][i];
                pars[6] = org.jlab.rec.dc.timetodistance.TableLoader.b1[s][i];
                pars[7] = org.jlab.rec.dc.timetodistance.TableLoader.b2[s][i];
                pars[8] = org.jlab.rec.dc.timetodistance.TableLoader.b3[s][i];
                pars[9] = org.jlab.rec.dc.timetodistance.TableLoader.b4[s][i];
                pars[10] = 2.*Constants.getInstance().wpdist[i];//fix dmax
                MinuitPar mnp = new MinuitPar();

                for(int p = 0; p < 10; p++) {
                    mnp.add(parNames[p], pars[p], errs[p]);
                }
                mnp.add(parNames[10], pars[10], errs[10]);
                for (int j = 0; j < T2DCalib.alphaBins; j++) {
                    for (int k = 0; k < BBins+1; k++) {
                        TvsDBr.put(new Coordinate(s, i,j,k), new FitLine("f_"+"i"+i+"i"+s+"j"+j+"k"+k, i, j, k, 
                                   mnp));
                        TvsDBr.get(new Coordinate(s, i, j, k)).useMidBfieldBin=true;
                        TvsDBr.get(new Coordinate(s, i, j, k)).useMidAlphaBin=true;
                        TvsDBr.get(new Coordinate(s, i, j, k)).setRange(0, pars[10]);
                        TvsDBr.get(new Coordinate(s, i, j, k)).setLineWidth(2);
                        TvsDBr.get(new Coordinate(s, i, j, k)).setLineColor(k+1);
                    }
                }
            }
        }
        
    }
    
    
    public static synchronized void Fill(IndexedTable t2dPressure, IndexedTable t2dPressRef, IndexedTable pressure) {
        System.out.println("****** T2D TABLE");	
        if (T2DLOADED) return;
        System.out.println("Filling T2D TABLE");
        org.jlab.rec.dc.timetodistance.TableLoader.Fill(t2dPressure, t2dPressRef, pressure);
        if(Boolean.parseBoolean(T2DViewer.updatedBConstants.getText())) {
            for(int s =0; s<6; s++) {
                for(int l = 2; l<4; l++) {
                    org.jlab.rec.dc.timetodistance.TableLoader.delta_bfield_coefficient[s][l]=0.14;
                    org.jlab.rec.dc.timetodistance.TableLoader.b1[s][l]=1.0;
                    org.jlab.rec.dc.timetodistance.TableLoader.b3[s][l]=11.5;
                }
            }
           
        }
        
        
        fillT2DGraphs();
        System.out.println(" T2D TABLE FILLED.....");
        //testBeq1();
        //test();
        T2DLOADED = true;
     }
    
    
    
    public static synchronized void ReFill() {
        //reset
        reset();
        System.out.println("RESET DONE....");
        org.jlab.rec.dc.timetodistance.TableLoader.DISTFROMTIME = new double[6][6][org.jlab.rec.dc.timetodistance.TableLoader.maxBinIdxB+1][org.jlab.rec.dc.timetodistance.TableLoader.maxBinIdxAlpha+1][org.jlab.rec.dc.timetodistance.TableLoader.betaValues.length][nBinsT]; // sector slyr alpha Bfield time bins [s][r][ibfield][icosalpha][tbin]      
        org.jlab.rec.dc.timetodistance.TableLoader.FillTable();
        if(t2dc.NbRunFit>0) {
            System.out.println("REMAKING T2D FUNCTION PLOTS....");
            refillT2DGraphs();
            System.out.println("T2D FUNCTION PLOTS DONE ");
        }
     }
    public static void fillT2DGraphs() {
        fitsinit();
        t2dc.getAnalysisCanvas().getCanvas("Fit Function").divide(4, 3);
        t2dc.getAnalysisCanvas().getCanvas("Fit Function").setAxisLabelSize(9);
        int s=0;
        for (int i = 0; i < 6; i++) {
           int cd1=0;
           if(i/2==0) {
               cd1=i;
           }
           if(i/2==1){
               cd1=i+2;
           }
           if(i/2==2){
               cd1=i+4;
           }

            for (int j = 0; j < T2DCalib.alphaBins; j++) {
                if(i<2 || i>3) {
                    t2dc.getAnalysisCanvas().getCanvas("Fit Function").cd(cd1);
                    t2dc.getAnalysisCanvas().getCanvas("Fit Function").draw(TvsDB.get(new Coordinate(s,i,j,8)), "same");
                } else {
                    t2dc.getAnalysisCanvas().getCanvas("Fit Function").cd(cd1);
                    t2dc.getAnalysisCanvas().getCanvas("Fit Function").draw(TvsDB.get(new Coordinate(s, i,j,0)), "same");
                    for (int k = 1; k < BBins; k++) {
                        t2dc.getAnalysisCanvas().getCanvas("Fit Function").draw(TvsDB.get(new Coordinate(s, i,j,k)), "same");
                    } 
                }
            }
        }
        t2dc.getAnalysisCanvas().getCanvas("Fit Function").getPad(0).getAxisY().setTitle("time (ns)");
        t2dc.getAnalysisCanvas().getCanvas("Fit Function").getPad(4).getAxisY().setTitle("time (ns)");
        t2dc.getAnalysisCanvas().getCanvas("Fit Function").getPad(8).getAxisY().setTitle("time (ns)");
        t2dc.getAnalysisCanvas().getCanvas("Fit Function").getPad(8).getAxisX().setTitle("doca (cm)");
        t2dc.getAnalysisCanvas().getCanvas("Fit Function").getPad(9).getAxisX().setTitle("doca (cm)");
        t2dc.getAnalysisCanvas().getCanvas("Fit Function").getPad(0).setTitle("with ccdb parameters");
        t2dc.getAnalysisCanvas().getCanvas("Fit Function").getPad(1).setTitle("with ccdb parameters");
        t2dc.getAnalysisCanvas().getCanvas("Fit Function").update();
    }
     public static void refillT2DGraphs() {
        reset();
        t2dc.getAnalysisCanvas().getCanvas("Fit Function").clear();
        t2dc.getAnalysisCanvas().getCanvas("Fit Function").divide(4, 3);
        t2dc.getAnalysisCanvas().getCanvas("Fit Function").setAxisLabelSize(9);
        //TColorPalette palette = new TColorPalette();
        //palette.setBackgroundColor(Color.yellow);
        System.out.println("RESETTING Fit Function Panel");
        int s = 0; 
        for (int i = 0; i < 6; i++) {
           int cd1=0;
           int cd2=0;
           if(i/2==0) {
               cd1=i;
               cd2=cd1+2;
           }
           if(i/2==1){
               cd1=i+2;
               cd2=cd1+2;
           }
           if(i/2==2){
               cd1=i+4;
               cd2=cd1+2;
           }
            for (int j = 0; j < T2DCalib.alphaBins; j++) {
                if(i<2 || i>3) {
                    t2dc.getAnalysisCanvas().getCanvas("Fit Function").cd(cd1);
                    t2dc.getAnalysisCanvas().getCanvas("Fit Function").draw(TvsDB.get(new Coordinate(s, i,j,8)), "same");
                    t2dc.getAnalysisCanvas().getCanvas("Fit Function").cd(cd2);
                    t2dc.getAnalysisCanvas().getCanvas("Fit Function").draw(TvsDBr.get(new Coordinate(s, i,j,8)), "same");
                } else {
                    t2dc.getAnalysisCanvas().getCanvas("Fit Function").cd(cd1);
                    t2dc.getAnalysisCanvas().getCanvas("Fit Function").draw(TvsDB.get(new Coordinate(s, i,j,0)), "same");
                    for (int k = 1; k < BBins; k++) {
                        t2dc.getAnalysisCanvas().getCanvas("Fit Function").draw(TvsDB.get(new Coordinate(s, i,j,k)), "same");
                    } 
                    t2dc.getAnalysisCanvas().getCanvas("Fit Function").cd(cd2);
                    t2dc.getAnalysisCanvas().getCanvas("Fit Function").draw(TvsDBr.get(new Coordinate(s, i,j,0)), "same");
                    for (int k = 1; k < BBins; k++) {
                        t2dc.getAnalysisCanvas().getCanvas("Fit Function").draw(TvsDBr.get(new Coordinate(s, i,j,k)), "same");
                    } 
                }
            }
        } 
        
        t2dc.getAnalysisCanvas().getCanvas("Fit Function").getPad(0).getAxisY().setTitle("time (ns)");
        t2dc.getAnalysisCanvas().getCanvas("Fit Function").getPad(4).getAxisY().setTitle("time (ns)");
        t2dc.getAnalysisCanvas().getCanvas("Fit Function").getPad(8).getAxisY().setTitle("time (ns)");
        t2dc.getAnalysisCanvas().getCanvas("Fit Function").getPad(8).getAxisX().setTitle("doca (cm)");
        t2dc.getAnalysisCanvas().getCanvas("Fit Function").getPad(9).getAxisX().setTitle("doca (cm)");
        t2dc.getAnalysisCanvas().getCanvas("Fit Function").getPad(10).getAxisX().setTitle("doca (cm)");
        t2dc.getAnalysisCanvas().getCanvas("Fit Function").getPad(11).getAxisX().setTitle("doca (cm)");
        t2dc.getAnalysisCanvas().getCanvas("Fit Function").getPad(0).setTitle("with ccdb parameters");
        t2dc.getAnalysisCanvas().getCanvas("Fit Function").getPad(1).setTitle("with ccdb parameters");
        t2dc.getAnalysisCanvas().getCanvas("Fit Function").getPad(2).setTitle("with fit parameters");
        t2dc.getAnalysisCanvas().getCanvas("Fit Function").getPad(3).setTitle("with fit parameters");
        t2dc.getAnalysisCanvas().getCanvas("Fit Function").update();
    }
    public static void test() {
        TimeToDistanceEstimator t2de = new TimeToDistanceEstimator();
        for(int r = 3; r<4; r++ ){ //loop over slys
            for(int ibfield =1; ibfield<2; ibfield++) {
                for(int icosalpha =org.jlab.rec.dc.timetodistance.TableLoader.maxBinIdxAlpha; icosalpha<org.jlab.rec.dc.timetodistance.TableLoader.maxBinIdxAlpha+1; icosalpha++) {
                    double cos30minusalpha = Math.cos(Math.toRadians(30.)) + (double) (icosalpha)*(1. - Math.cos(Math.toRadians(30.)))/5.;
                    double alpha = -(Math.toDegrees(Math.acos(cos30minusalpha)) - 30);
                    for(int ibeta=4; ibeta<5; ibeta++) {
                        for(int tbin = 0; tbin<org.jlab.rec.dc.timetodistance.TableLoader.maxTBin; tbin++) {
                            double time = 2*tbin; 
                            double doca = org.jlab.rec.dc.timetodistance.TableLoader.DISTFROMTIME[0][r][ibfield][icosalpha][ibeta][tbin];
                            double calctime = org.jlab.rec.dc.timetodistance.TableLoader.calc_Time( doca,  alpha, 1, 1, r+1) ;
                            double deltatime_beta = util.getDeltaTimeBeta(doca,org.jlab.rec.dc.timetodistance.TableLoader.betaValues[ibeta],org.jlab.rec.dc.timetodistance.TableLoader.distbeta[0][r],org.jlab.rec.dc.timetodistance.TableLoader.v0[0][r]);
                            calctime+=deltatime_beta;
                            double calcdoca = t2de.interpolateOnGrid(1, alpha, 1, time,  0, r);
                            if(time-calctime <4)System.out.println("alpha "+alpha+" time "+time+" time calculated from doca in T2B "+(float)calctime+" doca in T2B "+doca+" calculated doca from interp "+calcdoca);
                            if(time-calctime <4)System.out.println("interpolation for time range[" +(time -1)+","+(time+1)+"] calculated doca from interp "+(float)t2de.interpolateOnGrid(1, alpha, 1, time-1,  0, r)+", "+
                                    t2de.interpolateOnGrid(0, alpha, 1, time+1,  0, r));
                            
                        }
                    }
                }
            }
        }
    }
    /**
     * 
     * @param x distance to wire in cm
     * @param alpha local angle in deg
     * @param bfield B field value a x in T
     * @param sector sector  
     * @param superlayer superlayer 
     * @return returns time (ns) when given inputs of distance x (cm), local angle alpha (degrees) and magnitude of bfield (Tesla).  
     */
    private static FcnUtility util = new FcnUtility();
    
    private static int getAlphaBinT2DC(double alpha) {
        int v = -1;
        for(int i = 0; i<T2DCalib.AlphaValues.length; i++) {
            
            if(Math.abs(alpha-T2DCalib.AlphaValues[i])<T2DCalib.AlphaBinHalfWidth)
                v = i;
        } 
        
        return v;
    }
    
    
    
    
    public static void main(String[] args) {
//        System.setProperty("CLAS12DIR", "/Users/ziegler/BASE/Tracking/HighLumi/coatjava/coatjava/");
//        ConstantsManager ccdb = new ConstantsManager();
//        ccdb.init(Arrays.asList(new String[]{
//            "/geometry/dc/superlayer",
//            "/calibration/dc/time_to_distance/t2d_pressure", 
//            "/hall/weather/pressure",
//            "/calibration/dc/time_to_distance/ref_pressure",
//            "/calibration/dc/time_jitter"}));
//        int newRun=5700;
//        String var = "default";
//        Driver.init();
//        System.out.println("* VARIATION *"+var);
//        T2DViewer.ccdb.setVariation(var);
//        TableLoader.Fill(ccdb.getConstants(newRun, "/calibration/dc/time_to_distance/t2d_pressure"),
//               ccdb.getConstants(newRun, "/calibration/dc/time_to_distance/ref_pressure"),
//                ccdb.getConstants(newRun, "/calibration/dc/time_to_distance/ref_pressure"));  
//        
//        for (int i = 0; i < 6; i++) {
//            TvsDCan.get(new Coordinate(i)).cd(0);
//            TvsDBCan.get(new Coordinate(i)).cd(0);
//            for (int j = 0; j < maxBinIdxAlpha+1; j++) {
//                if(i<2 || i>3) {
//                    TvsDCan.get(new Coordinate(i)).draw(TvsD.get(new Coordinate(i,j,8)), "same");
//                    TvsDBCan.get(new Coordinate(i)).draw(TvsD.get(new Coordinate(i,j,8)), "same");
//                } else {
//                    for (int k = 1; k < BBins; k++) {
//                        TvsDCan.get(new Coordinate(i)).draw(TvsD.get(new Coordinate(i,j,k)), "same");
//                        TvsDBCan.get(new Coordinate(i)).draw(TvsD.get(new Coordinate(i,j,k)), "same");
//                    } 
//                }
//            }
//        }
//                            
    }
}
