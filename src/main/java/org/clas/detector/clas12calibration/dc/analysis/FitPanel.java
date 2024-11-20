/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.dc.analysis;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;
import javax.swing.plaf.metal.MetalButtonUI;
import org.clas.detector.clas12calibration.dc.calt2d.T2DCalib;
import org.freehep.math.minuit.MnUserParameters;


public class FitPanel {
    
    
    private Map<Integer, ArrayList<Double>> pars     = new HashMap<Integer, ArrayList<Double>>();
    private double[]          range     = new double[2];
    private JFrame            frame     = new JFrame();
    private CustomPanel2      panel     = null;            
    private T2DCalib _pM;

    private int getSectorLayer(int i, int s) {
        return s*6+i;
    }
    public FitPanel(T2DCalib pM) {
        this._pM = pM;
        //init pars container
        for(int s = 0; s<6; s++) {
            for(int j = 0; j<6; j++) {
                pars.put(getSectorLayer(j,s), new ArrayList<Double>());
            }
        }
    }
    
    public void openFitPanel(String title, Map<Coordinate, MnUserParameters> TvstrkdocasFitPars){
        frame    = new JFrame();
        panel = new CustomPanel2(TvstrkdocasFitPars);
        frame.setSize(350, 300); 
        frame.setTitle(title);
        frame.add(panel);
        frame.pack();
        frame.setVisible(true);
//        frame.setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
            
    }
    public void setGreenFitButton() {
        panel.fitButton.setBackground(Color.GREEN);
        panel.fitButton.setVisible(true);
        panel.repaint();
    }
    public void setRedFitButton() {
        panel.fitButton.setBackground(Color.RED);
        panel.fitButton.setVisible(true);
        panel.repaint();
    }
    public boolean fitted = false;
    public boolean initscan = false;
    public void refit(Map<Coordinate, MnUserParameters> TvstrkdocasFitPars) throws FileNotFoundException{
        System.out.println("READY TO RUN THE FIT "+initscan);
        if(initscan==false) return;
        if(!this._pM.useBProf) {
            for(int s=0; s<6; s++) {
                for(int j = 2; j<4; j++) {
                    for(int i=5; i<10; i++){ 
                        if(panel.fixFit[i][j][s].isSelected()==false) {
                            panel.fixFit[i][j][s].setSelected(true);
                            System.out.println("REDO THE SEGMENT FITS BEFORE RELEASING THE B-field DEPENDENT PARAMETERS!");
                        }
                    }
                }
            }
        }
        boolean[][][] fixedPars = new boolean[10][6][6];
        for(int s = 0; s<6; s++) {
            for(int j = 0; j<6; j++) {
                pars.get(this.getSectorLayer(j, s)).clear();
            }
        }
        int npar = 10;
        for(int s = 0; s<6; s++) {
            for(int j = 0; j<6; j++) {
                for(int i=0; i<npar; i++){   
                    if(panel.params[i][j][s].getText().isEmpty()){
                        this.pars.get(this.getSectorLayer(j, s)).add(TvstrkdocasFitPars.get(new Coordinate(j)).value(i));
                    }
                    else { 
                        this.pars.get(this.getSectorLayer(j, s)).add(Double.parseDouble(panel.params[i][j][s].getText()));
                    }
                }
            }
        }
        if(!panel.minRange.getText().isEmpty())this.range[0] = Double.parseDouble(panel.minRange.getText());
        else this.range[0] = 0.0;
        if(!panel.maxRange.getText().isEmpty())this.range[1] = Double.parseDouble(panel.maxRange.getText());
        else this.range[1] = 2.0;
        for(int s = 0; s<6; s++) {
            for(int j = 0; j<6; j++) {    
                for(int i=0; i<npar; i++){
                    TvstrkdocasFitPars.get(new Coordinate(j)).setValue(i,this.pars.get(this.getSectorLayer(j, s)).get(i));
                }
            }
        }
        this._pM.initFitParsToFile();
        
        for(int s = 0; s<6; s++) {
            for(int j = 0; j<6; j++) {
                for(int i=0; i<npar; i++){ 
                    if(panel.fixFit[i][j][s].isSelected()==true)
                        fixedPars[i][j][s] = true;
                }
            }
        }
        //Don't allow to fit the B>0 profile if they are not filled
        
        this._pM.runFit(fixedPars);
         
            //this._pM.runParamScan(j, fixedPars);
        for(int s = 0; s<6; s++) {        
            for(int j0 = 0; j0<6; j0++) { 
                int j = j0+6*s;
                for(int p = 0; p<10; p++) {
                    panel.pars[p][j0][s] = TvstrkdocasFitPars.get(new Coordinate(j)).value(p);
                    if(p!=3) {
                        panel.params[p][j0][s].setText(String.format("%.5f",TvstrkdocasFitPars.get(new Coordinate(j)).value(p)));
                    } else {
                        panel.params[p][j0][s].setText(String.format("%.3f",TvstrkdocasFitPars.get(new Coordinate(j)).value(p)));
                    }
                }
                panel.fixFit[2][j0][s].setSelected(true);
                panel.fixFit[4][j0][s].setSelected(true);
            }
        }
        fitted = true;
        this._pM.plotFits(fitted);
        this._pM.pw3.close();
    }
    
    public void parscan(Map<Coordinate, MnUserParameters> TvstrkdocasFitPars) throws FileNotFoundException{
        if(this._pM.eventProcessingDone==false) {
            System.out.println("PATIENCE ... WAIT UNTIL THE EVENT PROCESSING IS DONE....");
            return;
        }
        boolean[][][] fixedPars = new boolean[10][6][6];
        
        int npar = 10;
        for(int s = 0; s<6; s++) {
            for(int j0 = 0; j0<6; j0++) {
                int j = j0+6*s;
                pars.get(j).clear();
                for(int i=0; i<npar; i++){   
                    if(panel.params[i][j0][s].getText().isEmpty()){
                        this.pars.get(j).add(TvstrkdocasFitPars.get(new Coordinate(j)).value(i));
                    }
                    else { 
                        this.pars.get(j).add(Double.parseDouble(panel.params[i][j0][s].getText()));
                    }
                }
            }
        }
        if(!panel.minRange.getText().isEmpty())this.range[0] = Double.parseDouble(panel.minRange.getText());
        else this.range[0] = 0.0;
        if(!panel.maxRange.getText().isEmpty())this.range[1] = Double.parseDouble(panel.maxRange.getText());
        else this.range[1] = 2.0;
        for(int s = 0; s<6; s++) {
            for(int j0 = 0; j0<6; j0++) {
                int j = j0+6*s;    
                for(int i=0; i<npar; i++){
                    TvstrkdocasFitPars.get(new Coordinate(j)).setValue(i,this.pars.get(j).get(i));
                
                    if(panel.fixFit[i][j0][s].isSelected()==true)
                        fixedPars[i][j0][s] = true;
                }
            }
        }
        if(panel.runIndivSectors.isSelected()==true) {
            T2DCalib.minSec=0;
            T2DCalib.maxSec=6;
        } else {
            T2DCalib.minSec=6;
            T2DCalib.maxSec=7;
        }
        this._pM.runParamScan(fixedPars);
        initscan=true;
        for(int s = 0; s<6; s++) {
            for(int j0 = 0; j0<6; j0++) {    
                int j = j0+6*s;
                for(int p = 0; p<10; p++) {
                    panel.pars[p][j0][s] = TvstrkdocasFitPars.get(new Coordinate(j)).value(p);
                    if(p!=3) {
                        panel.params[p][j0][s].setText(String.format("%.5f",TvstrkdocasFitPars.get(new Coordinate(j)).value(p)));
                    } else {
                        panel.params[p][j0][s].setText(String.format("%.3f",TvstrkdocasFitPars.get(new Coordinate(j)).value(p)));
                    }
                }
                panel.fixFit[2][j0][s].setSelected(true);
                panel.fixFit[4][j0][s].setSelected(true);
            }
        }
        fitted = true;
        this._pM.plotFits(fitted);
        
    }
    
    public void plotResiduals() {
        this._pM.rePlotResi();
    }
    public void reCook() {
        if(panel.useBprof.isSelected()==true)
            this._pM.useBProf=true;
        this._pM.reCook();
    }
    public void reset() {
        this._pM.resetPars();
        setGreenFitButton();
        panel.updateUI();
        
    }
    public void savePars() throws FileNotFoundException {
        this._pM.plotFits(true);
    }
    private final class CustomPanel2 extends JPanel {
        JLabel label;
        JPanel panel;
        
        JTabbedPane tpane;
        JCheckBox useBprof;
        JCheckBox runIndivSectors;
    	JTextField minRange = new JTextField(5);
	JTextField maxRange = new JTextField(5);
	JTextField[][][] params = new JTextField[10][6][6];
        JCheckBox[][][]  fixFit = new JCheckBox[10][6][6];;
        JPanel[] secpanels = new JPanel[6];   
        
        JButton   scanButton = null;
        JButton   fitButton = null;
        JButton   resetButton = null;
        JButton   saveToFileButton = null;
        JButton   resButton = null;
        JButton   reCookButton = null;
        
        String[] parNames = new String[] {"v0", "vmid", "R", "tmax", "distbeta", "delBf", "b1", "b2", "b3", "b4"};
        double[][][] pars = new double[10][6][6];
        
        private void addSectorPanel(int s) {
            secpanels[s] = addPanel(s);
        }
        private JPanel addPanel(int s) {
            int npar = 10;
            JPanel spanel;            
            spanel = new JPanel(new GridLayout(npar+1, 6));
            for (int i = 0; i < 6; i++) {
                String SuperLayer = Integer.toString(i + 1);
                spanel.add(new JLabel(""));
                spanel.add(new JLabel("Superlayer " + SuperLayer));
            }
            spanel.add(new JLabel(""));
            for (int i = 0; i < npar; i++) {  
                JLabel l = new JLabel("      "+parNames[i], JLabel.LEADING);
                spanel.add(l);
                for (int j = 0; j < 6; j++) {
                    fixFit[i][j][s] = new JCheckBox("Fix");
                    // aa is true for parameters "R", "distbeta", "b1", "b2", "b3", and "b4"
                    //boolean aa = i==2 || i==4 || i>5;
                    // aa is true for parameters "b1", "b2", "b3", and "b4"
                    boolean aa = i>5;
                    // bb is true for parameter "tmax" on superlayers 3 and 4 only
                    boolean bb = i==3 && j==6; //(j==2 || j==3);
                    // cc is true for parameter "delBf" on superlayers 1, 2, 5, and 6 only
                    boolean cc = i==5; // && (j==0 || j==1 || j==4 || j==5);
                    if(aa || bb || cc) {
                        fixFit[i][j][s].setSelected(true);
                    } else {
                        fixFit[i][j][s].setSelected(false);
                    }
                    
                    params[i][j][s] = new JTextField(3);
                    if(i!=3) {
                        params[i][j][s].setText(String.format("%.5f", pars[i][j][s]));
                    } else {
                        params[i][j][s].setText(String.format("%.3f", pars[i][j][s]));
                    }
                    spanel.add(params[i][j][s]);
                    spanel.add(fixFit[i][j][s]);
                }
            }
            spanel.setVisible(true);
            return spanel;
        }
        
        public CustomPanel2(Map<Coordinate, MnUserParameters> TvstrkdocasFitPars) { 
            super(new BorderLayout());
            for(int s =0; s<6; s++) {
                for(int i = 0; i < 6; i++) {
                    for(int p = 0; p<10; p++) {
                        pars[p][i][s] = TvstrkdocasFitPars.get(new Coordinate(i+6*s)).value(p);
                    }
                }
            }
            panel = new JPanel(new GridLayout(1, 1));
            panel.setBackground(Color.white);
            tpane = new JTabbedPane();
            for(int s=0; s<6; s++) {
                this.addSectorPanel(s);
                tpane.addTab("Sector "+(s+1), secpanels[s]);
            }
            panel.add(tpane);
            JPanel settings = new JPanel(new GridLayout(1, 6));
            settings.add(new JLabel("    Fit range min"));
            minRange.setText(Double.toString(0));
            settings.add(minRange);
            settings.add(new JLabel("    Fit range max"));
            maxRange.setText(Double.toString(2.0));
            settings.add(maxRange);
            settings.add(new JLabel("                  "));
            useBprof = new JCheckBox(" Use B>0 profiles");
            useBprof.setSelected(false);
            settings.add(useBprof);
            settings.add(new JLabel("                  "));
            runIndivSectors = new JCheckBox(" Fit by sector");
            settings.add(runIndivSectors);
            
            JPanel buttonsPanel = new JPanel();
            buttonsPanel.setLayout(new GridLayout(1, 4));

            saveToFileButton = new JButton("SAVE PARAMETERS");
            saveToFileButton.setUI(new MetalButtonUI());
            saveToFileButton.setBackground(Color.PINK);
            saveToFileButton.setContentAreaFilled(true);
            saveToFileButton.setOpaque(true);
            saveToFileButton.setFont(bBold);
            saveToFileButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    try {
                        savePars();
                    } catch (FileNotFoundException ex) {
                        Logger.getLogger(FitPanel.class.getName()).log(Level.SEVERE, null, ex);
                    }
                    return;
                }

                
            });
            
            resetButton = new JButton("RESET PARAMETERS");
            resetButton.setUI(new MetalButtonUI());
            resetButton.setBackground(Color.CYAN);
            resetButton.setContentAreaFilled(true);
            resetButton.setOpaque(true);
            resetButton.setFont(bBold);
            resetButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    reset();
                    return;
                }
            });
            
            resButton = new JButton("PLOT RESIDUALS");
            resButton.setUI(new MetalButtonUI());
            resButton.setBackground(Color.YELLOW);
            resButton.setContentAreaFilled(false);
            resButton.setOpaque(true);
            resButton.setFont(bBold);
            resButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    plotResiduals();
                    return;
                }
            });
            
            reCookButton = new JButton("REDO SEGMENT FITS");
            reCookButton.setUI(new MetalButtonUI());
            reCookButton.setBackground(Color.ORANGE);
            reCookButton.setContentAreaFilled(false);
            reCookButton.setOpaque(true);
            reCookButton.setFont(bBold);
            reCookButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    //fitButton.setContentAreaFilled(false);
                    //fitButton.setOpaque(true);
                    //setRedFitButton();
                    fitButton.setBackground(Color.RED);
                    
                    System.out.println("******************************************");
                    System.out.println("* REFITTING SEGMENTS WITH NEW PARAMETERS *");
                    System.out.println("******************************************");
                    
                    reCook();
                    System.out.println("******************************************");
                    System.out.println("*      READY TO REFIT THE HISTOGRAMS     *");
                    System.out.println("******************************************");
                    
                    //setGreenFitButton();
                    fitButton.setBackground(Color.GREEN);
                    return;
                }
            });
            
            fitButton = new JButton("FIT TIME TO DISTANCE");
            fitButton.setUI(new MetalButtonUI());
            fitButton.setBackground(Color.RED);
            fitButton.setContentAreaFilled(false);
            fitButton.setOpaque(true);
            fitButton.setFont(bBold);
            fitButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    try {
                        refit(TvstrkdocasFitPars);
                    } catch (FileNotFoundException ex) {
                        Logger.getLogger(FitPanel.class.getName()).log(Level.SEVERE, null, ex);
                    }
                    return;
                }
            });

            scanButton = new JButton("INITIAL PARAMETER SCAN");
            scanButton.setUI(new MetalButtonUI());
            scanButton.setBackground(Color.WHITE);
            scanButton.setBorderPainted(true);
            scanButton.setContentAreaFilled(false);
            scanButton.setOpaque(true);
            scanButton.setFont(bBold);
            scanButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    try {
                        parscan(TvstrkdocasFitPars);
                    } catch (FileNotFoundException ex) {
                        Logger.getLogger(FitPanel.class.getName()).log(Level.SEVERE, null, ex);
                    }
                    return;
                }
            });
            
            buttonsPanel.add(scanButton);
            buttonsPanel.add(fitButton);
            buttonsPanel.add(reCookButton);
            buttonsPanel.add(resButton);
            buttonsPanel.add(resetButton);
            buttonsPanel.add(saveToFileButton);
            
            
            this.add(panel, BorderLayout.PAGE_START);
            this.add(settings, BorderLayout.CENTER);
            this.add(buttonsPanel, BorderLayout.PAGE_END);
            
            label = new JLabel("Click the \"Show it!\" button"
                           + " to bring up the selected dialog.",
                           JLabel.CENTER);
            
        }
        private Font bBold = new Font("Arial", Font.BOLD, 16);
        void setLabel(String newText) {
            label.setText(newText);
        }

    }
}
