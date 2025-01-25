/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.clas.detector.clas12calibration.dc.calt2d;

import java.io.PrintWriter;
import java.util.Map;
import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import org.clas.detector.clas12calibration.dc.calt2d.FitUtility.MinuitPar;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.BBins;
import static org.clas.detector.clas12calibration.dc.calt2d.T2DCalib.alphaBins;
import org.jlab.detector.calib.utils.CalibrationConstants;
import org.jlab.groot.base.TColorPalette;
import org.jlab.groot.data.DataLine;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;
import org.jlab.groot.math.F1D;

/**
 *
 * @author ziegler
 */
public class PlotUtility {
    
   public static void Plot(int s, int i, int j, Map<Coordinate, H2F> Tvstrkdocas, Map<Coordinate, H2F> Tvscalcdocas, 
            Map<Coordinate, H2F> Tresvstrkdocas, Map<Coordinate, GraphErrors> TvstrkdocasProf, 
            Map<Coordinate, FitLine> TvstrkdocasFits, Map<Coordinate, GraphErrors> TvstrkdocasInit, 
            EmbeddedCanvasTabbed analysisCanvas, int iterationNum) {
        DataLine l = new DataLine(0, 0, 2.5, 0);
        l.setLineStyle(2);
        l.setLineColor(2);
        if(i<2 || i>3) { // regions 1 and 3 --> no b-field
            if(TvstrkdocasProf.get(new Coordinate(s, i, j, BBins)).getVectorX().size()>0) {
                analysisCanvas.getCanvas("TrackDoca vs T").cd(0);
                analysisCanvas.getCanvas("TrackDoca vs T").draw(Tvstrkdocas.get(new Coordinate(s, i, j, BBins)));
                analysisCanvas.getCanvas("TrackDoca vs T Graphs").cd(0);
                analysisCanvas.getCanvas("TrackDoca vs T Graphs").draw(Tvstrkdocas.get(new Coordinate(s, i, j, BBins)));
                analysisCanvas.getCanvas("CalcDoca vs T").cd(0);
                analysisCanvas.getCanvas("CalcDoca vs T").draw(Tvscalcdocas.get(new Coordinate(s, i, j, BBins)));
                analysisCanvas.getCanvas("TrackDoca vs T Graphs").
                        draw(TvstrkdocasProf.get(new Coordinate(s, i, j, BBins)), "same");
                //TvstrkdocasFits.get(new Coordinate(i, j, BBins)).setRange(0, maxx[i]);
                analysisCanvas.getCanvas("TrackDoca vs T Graphs").
                        draw(TvstrkdocasFits.get(new Coordinate(s, i, j, BBins)), "same");
                analysisCanvas.getCanvas("CalcDoca vs T").
                        draw(TvstrkdocasFits.get(new Coordinate(s, i, j, BBins)), "same");
                //Resi canvas
                //analysisCanvas.getCanvas("TrackDoca vs T Fit Resi").cd(0);
                
                GraphErrors g1 = new GraphErrors();
                GraphErrors g2 = new GraphErrors();
                g2.setTitle(TvstrkdocasInit.get(new Coordinate(s, i, j, BBins)).getTitle());
                g1.copy(TvstrkdocasProf.get(new Coordinate(s, i, j, BBins)));
                for(int ip =0; ip<g1.getVectorX().getSize(); ip++) {
                    if(g1.getDataEY(ip)!=0) {
                        double yf = TvstrkdocasFits.get(new Coordinate(s, i, j, BBins)).evaluate(g1.getDataX(ip));
                        double y = g1.getDataY(ip);
                        g2.addPoint(g1.getDataX(ip), y-yf, 0, g1.getDataEY(ip));
                        if(iterationNum==0) {
                            TvstrkdocasInit.get(new Coordinate(s, i, j, BBins)).addPoint(g1.getDataX(ip), y-yf, 0, 0);
                        }
                    }       
                }
                analysisCanvas.getCanvas("TrackDoca vs T Fit Resi").clear();
                analysisCanvas.getCanvas("TrackDoca vs T Fit Resi").cd(0);
                analysisCanvas.getCanvas("TrackDoca vs T Fit Resi").getPad().setPalette(TColorPalette.PaletteName.kRainBow);
                analysisCanvas.getCanvas("TrackDoca vs T Fit Resi").draw(Tresvstrkdocas.get(new Coordinate(s, i, j, BBins)));
                analysisCanvas.getCanvas("TrackDoca vs T Fit Resi").draw(g2, "Esame");
                analysisCanvas.getCanvas("TrackDoca vs T Fit Resi").draw(TvstrkdocasInit.get(new Coordinate(s, i, j, BBins)), "Esame");                   
                analysisCanvas.getCanvas("TrackDoca vs T Fit Resi").draw(l);
            }
        } else {   
            //plot the profiles for the various B-field components
            
            analysisCanvas.getCanvas("TrackDoca vs T").cd(0);
            analysisCanvas.getCanvas("TrackDoca vs T").draw(Tvstrkdocas.get(new Coordinate(s, i, j, BBins)));
            analysisCanvas.getCanvas("TrackDoca vs T Graphs").cd(0);
            analysisCanvas.getCanvas("TrackDoca vs T Graphs").draw(Tvstrkdocas.get(new Coordinate(s, i, j, BBins)));
            analysisCanvas.getCanvas("CalcDoca vs T").cd(0);
            analysisCanvas.getCanvas("CalcDoca vs T").draw(Tvscalcdocas.get(new Coordinate(s, i, j, BBins)));    
            int maxBbin=0;
            for(int k = 0; k < BBins; k++) {
                if(TvstrkdocasProf.get(new Coordinate(s, i, j, k)).getVectorX().size()>0){
                    maxBbin++;
                    analysisCanvas.getCanvas("TrackDoca vs T Graphs").
                            draw(TvstrkdocasProf.get(new Coordinate(s, i, j, k)), "same");
                    //TvstrkdocasFits.get(new Coordinate(i, j, k)).setRange(0, maxx[i]);
                    analysisCanvas.getCanvas("TrackDoca vs T Graphs").
                            draw(TvstrkdocasFits.get(new Coordinate(s, i, j, k)), "same");
                    analysisCanvas.getCanvas("CalcDoca vs T").
                    draw(TvstrkdocasFits.get(new Coordinate(s, i, j, k)), "same");
                }
            }
           
            
            if(TvstrkdocasProf.get(new Coordinate(s, i, j, 0)).getVectorX().size()>0) {
            //Resi canvas
                GraphErrors g1 = new GraphErrors();
                GraphErrors g2 = new GraphErrors();
                g2.setTitle(TvstrkdocasInit.get(new Coordinate(s, i, j, 0)).getTitle());
                g1.copy(TvstrkdocasProf.get(new Coordinate(s, i, j, 0)));
                for(int ip =0; ip<g1.getVectorX().getSize(); ip++) {
                    if(g1.getDataEY(ip)!=0) {
                        double yf = TvstrkdocasFits.get(new Coordinate(s, i, j, 0)).evaluate(g1.getDataX(ip));
                        double y = g1.getDataY(ip);
                        g2.addPoint(g1.getDataX(ip), y-yf, 0, g1.getDataEY(ip));
                        if(iterationNum==0) {
                            TvstrkdocasInit.get(new Coordinate(s, i, j, 0)).addPoint(g1.getDataX(ip), y-yf, 0, 0);
                        }
                    }       
                }
                analysisCanvas.getCanvas("TrackDoca vs T Fit Resi").clear();
                analysisCanvas.getCanvas("TrackDoca vs T Fit Resi").cd(0);
                analysisCanvas.getCanvas("TrackDoca vs T Fit Resi").getPad().setPalette(TColorPalette.PaletteName.kRainBow);
                analysisCanvas.getCanvas("TrackDoca vs T Fit Resi").draw(Tresvstrkdocas.get(new Coordinate(s, i, j, 0)));
                analysisCanvas.getCanvas("TrackDoca vs T Fit Resi").draw(g2, "Esame");
                analysisCanvas.getCanvas("TrackDoca vs T Fit Resi").draw(TvstrkdocasInit.get(new Coordinate(s, i, j, 0)), "Esame");                   
                analysisCanvas.getCanvas("TrackDoca vs T Fit Resi").draw(l);
            }
        }   
    }

    public static void rePlotResi(EmbeddedCanvasTabbed analysisCanvas, Map<Coordinate, H1F> timeResiFromFile, 
            Map<Coordinate, H1F> timeResi, Map<Coordinate, H1F> timeResiNew) {
        analysisCanvas.getCanvas("Time Residuals").clear();
        analysisCanvas.getCanvas("Time Residuals").divide(6, 3);
        for(int s=0; s<6; s++) {
            String st = "sector ";
            st+=s+1;
            analysisCanvas.getCanvas(st).clear();
            analysisCanvas.getCanvas(st).divide(6, 3);
        }
        //
        for(int i = 0; i<6; i++) {
            analysisCanvas.getCanvas("Time Residuals").cd(i);
            fitTimeResPlot(timeResiFromFile.get(new Coordinate(6,i)), 
                    analysisCanvas.getCanvas("Time Residuals"));
            analysisCanvas.getCanvas("Time Residuals").cd(i+6);
            fitTimeResPlot(timeResi.get(new Coordinate(6,i)), 
                    analysisCanvas.getCanvas("Time Residuals"));
            analysisCanvas.getCanvas("Time Residuals").cd(i+12);
            fitTimeResPlot(timeResiNew.get(new Coordinate(6,i)), 
                    analysisCanvas.getCanvas("Time Residuals"));
        }
        //sectors
        for(int s =0; s<6; s++) {
            String st = "sector ";
            st+=s+1;
            for(int i = 0; i<6; i++) {
                analysisCanvas.getCanvas(st).cd(i);
                fitTimeResPlot(timeResiFromFile.get(new Coordinate(s,i)), 
                        analysisCanvas.getCanvas(st));
                analysisCanvas.getCanvas(st).cd(i+6);
                fitTimeResPlot(timeResi.get(new Coordinate(s,i)), 
                        analysisCanvas.getCanvas(st));
                analysisCanvas.getCanvas(st).cd(i+12);
                fitTimeResPlot(timeResiNew.get(new Coordinate(s,i)), 
                        analysisCanvas.getCanvas(st));
            }
        }
    
        
    }
    
    
    public static void fitTimeResPlot(H1F h1, EmbeddedCanvas canvasRes) {
        if (h1==null) return;
        F1D gaus1Func = new F1D("gaus1Func", "[amp]*gaus(x,[mean],[sigma])", -0.3, 0.3); 
        gaus1Func.setParameter(0, h1.getMax());
        gaus1Func.setParameter(1, 0.0);
        gaus1Func.setParameter(2, 0.05);
        DataFitter.fit(gaus1Func, h1, "Q");
        
        int binMax = h1.getMaximumBin(); 
        double hm = h1.getDataY(binMax)/3; 
        double xlo = 0; 
        double xhi =0; 
        for(int b = 0; b<binMax-1; b++) {
            if(h1.getDataY(b)<=hm && h1.getDataY(b+1) >hm) {
                xlo = h1.getDataX(b);
                break;
            }
                
        }
        
        for(int b = binMax; b<99; b++) {
            if(h1.getDataY(b+1)<=hm && h1.getDataY(b) >hm) {
                xhi = h1.getDataX(b);
                break;
            }
                
        }
        gaus1Func.setRange(xlo, xhi);
        DataFitter.fit(gaus1Func, h1, "Q");
        
        //refit using a double gaussian 
        F1D gausFunc = new F1D("gausFunc", "[amp]*gaus(x,[mean],[sigma])+0.25*[amp]*gaus(x,[mean],[sigma2])", -0.3, 0.3); 
        gausFunc.setLineColor(4);
        gausFunc.setLineStyle(1);
        gausFunc.setLineWidth(2);
        gausFunc.setParameter(0, gaus1Func.getParameter(0));
        gausFunc.setParameter(1, gaus1Func.getParameter(1));
        gausFunc.setParameter(2, gaus1Func.getParameter(2)*0.75);
        gausFunc.setParameter(3, gaus1Func.getParameter(2));
        gausFunc.setOptStat(11100);
        h1.setOptStat(11); //only number of entries
        //canvasRes.clear();
        
        DataFitter.fit(gausFunc, h1, "Q");
        //gausFunc.setOptStat(101100);
        //gausFunc.setOptStat(101100); //mean and both sigmas
        
        
        int effSig = (int) Math.round(10000*Math.sqrt(gausFunc.getParameter(2)*gausFunc.getParameter(2)+
                                            0.25*0.25*gausFunc.getParameter(3)*gausFunc.getParameter(3))/
                                            Math.sqrt(1+0.25*0.25));
        String t = "Eff. Sig. "+effSig + "microns";
        h1.setTitle(t);
        canvasRes.draw(h1, "E1");
        //canvasRes.draw(gausFunc, "same");
    }

    public static void plotFits(PrintWriter pw, PrintWriter pw2, int NbRunFit, Map<Coordinate, H1F> ParsVsIter, 
            Map<Coordinate, MinuitPar> TvstrkdocasFitPars, 
            Map<Coordinate, FitLine> TvstrkdocasFits, Map<Coordinate, GraphErrors> TvstrkdocasProf, CalibrationConstants calib,
            EmbeddedCanvas canvas) {
        int ij=0;
        int ip =0;
            NbRunFit++;
            for(int s=0; s<6;s++) {
                for (int i = 0; i < 6; i++) {
               
                    for(int p = 0; p<6; p++) {
                        ParsVsIter.get(new Coordinate(i,p)).setBinContent(NbRunFit, TvstrkdocasFitPars.get(new Coordinate(0,i)).value(p));
                        ParsVsIter.get(new Coordinate(i,p)).setBinError(NbRunFit, TvstrkdocasFitPars.get(new Coordinate(0,i)).error(p));
                        ParsVsIter.get(new Coordinate(i,p)).setOptStat(0);
                        canvas.cd(ip);
                        
                        double min = ParsVsIter.get(new Coordinate(i,p)).getMin()-2*ParsVsIter.get(new Coordinate(i,p)).getBinError(1);
                        double max = ParsVsIter.get(new Coordinate(i,p)).getMax()+2*ParsVsIter.get(new Coordinate(i,p)).getBinError(1);
                        if(Math.abs(min)<1.e-06 && Math.abs(max)<1.e-06) {
                            min = -0.1;
                            max = 0.1;
                        }

                        canvas.draw(ParsVsIter.get(new Coordinate(i,p)));

                        ip++;
                    }
                }
            }    
            
            for(int s=0; s<6;s++) {
                for (int i = 0; i < 6; i++) {
                    for (int j = 0; j < alphaBins; j++) {
                        if(i<2 || i>3) {
                            if(TvstrkdocasProf.get(new Coordinate(s,i, j, BBins)).getVectorX().size()>0) {
                                FitUtility.updateTable(s, i, j, calib, TvstrkdocasFitPars);
                                TvstrkdocasFits.put(new Coordinate(s,i,j,BBins), new FitLine("f"+""+s+""+i+""+j+"0",s, i, j, BBins, 
                                TvstrkdocasFitPars.get(new Coordinate(s,i))));
                                TvstrkdocasFits.get(new Coordinate(s,i, j, BBins)).setLineStyle(4);
                                TvstrkdocasFits.get(new Coordinate(s,i, j, BBins)).setLineWidth(5);
                                TvstrkdocasFits.get(new Coordinate(s,i, j, BBins)).setLineColor(8);
                            } else {
                                //this.resetTable(i,j);
                            }

                        } else {
                            for(int k = 0; k < BBins; k++) { 
                                if(TvstrkdocasProf.get(new Coordinate(s,i, j, k)).getVectorX().size()>0){
                                    FitUtility.updateTable(s, i, j, calib, TvstrkdocasFitPars);
                                    TvstrkdocasFits.put(new Coordinate(s,i,j,k), new FitLine("f"+""+s+""+i+""+j+""+k, s, i, j, k, 
                                    TvstrkdocasFitPars.get(new Coordinate(s,i))));
                                    TvstrkdocasFits.get(new Coordinate(s, i, j, k)).setLineStyle(4);
                                    TvstrkdocasFits.get(new Coordinate(s,i, j, k)).setLineWidth(5);
                                    TvstrkdocasFits.get(new Coordinate(s,i, j, k)).setLineColor(k+1);
                                } else {
                                    //this.resetTable(i,j);
                                }
                            }
                        }
                        ij++;
                    }
                }
            }
            calib.fireTableDataChanged();    
            for(int s = 0; s < 6; s++) {
                for(int i = 0; i<6; i++) {
                    pw.printf("%d\t %d\t %d\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\n",
                    (s+1), (i+1), 0,
                    TvstrkdocasFitPars.get(new Coordinate(s,i)).value(0), //vo
                    0.0,
                    0.0,
                    TvstrkdocasFitPars.get(new Coordinate(s,i)).value(1), //vmid
                    0.0,
                    0.0,
                    TvstrkdocasFitPars.get(new Coordinate(s,i)).value(3), //tmax
                    0.0,
                    0.0,
                    TvstrkdocasFitPars.get(new Coordinate(s,i)).value(4), //distbeta
                    0.0,
                    0.0,
                    TvstrkdocasFitPars.get(new Coordinate(s,i)).value(5), //delta_bfield
                    0.0,
                    0.0,
                    TvstrkdocasFitPars.get(new Coordinate(s,i)).value(6), //b1
                    0.0,
                    0.0,
                    TvstrkdocasFitPars.get(new Coordinate(s,i)).value(7), //b2
                    0.0,
                    0.0,
                    TvstrkdocasFitPars.get(new Coordinate(s,i)).value(8), //b3
                    0.0,
                    0.0,
                    TvstrkdocasFitPars.get(new Coordinate(s,i)).value(9), //b4
                    0.0,
                    0.0,
                    TvstrkdocasFitPars.get(new Coordinate(s,i)).value(2), //R
                    0.0,
                    0.0);
                    
                    pw2.printf("%d\t %d\t %d\t %.6f\t +/- %.6f\t %.6f\t +/- %.6f\t %.6f\t +/- %.6f\t %.6f\t +/- %.6f\t %.6f\t +/- %.6f\n",
                            (s+1), (i+1), 0,
                            //v0
                            TvstrkdocasFitPars.get(new Coordinate(s,i)).value(0),
                            TvstrkdocasFitPars.get(new Coordinate(s,i)).error(0),
                            //tmax
                            TvstrkdocasFitPars.get(new Coordinate(s,i)).value(3),
                            TvstrkdocasFitPars.get(new Coordinate(s,i)).error(3),
                            //vmid
                            TvstrkdocasFitPars.get(new Coordinate(s,i)).value(1),
                            TvstrkdocasFitPars.get(new Coordinate(s,i)).error(1),
                            //deltaBf
                            TvstrkdocasFitPars.get(new Coordinate(s,i)).value(5),
                            TvstrkdocasFitPars.get(new Coordinate(s,i)).error(5),
                            //distbeta
                            TvstrkdocasFitPars.get(new Coordinate(s,i)).value(4),
                            TvstrkdocasFitPars.get(new Coordinate(s,i)).error(4));
                            
                    
                }
            }
            pw.close(); 
            pw2.close(); 
            //this.rePlotResi();
    }
}
