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
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.plaf.metal.MetalButtonUI;
import org.clas.detector.clas12calibration.dc.calt2d.T2DCalib;
import org.freehep.math.minuit.MnUserParameters;


public class FitPanel {
    
    
    private Map<Integer, ArrayList<Double>> pars     = new HashMap<Integer, ArrayList<Double>>();
    private double[]          range    = new double[2];
    private JFrame            frame    = new JFrame();
    private CustomPanel2      panel    = null;
    private T2DCalib _pM;

    public FitPanel(T2DCalib pM) {
        this._pM = pM;
        //init pars container
        for(int j = 0; j<6; j++) {
            pars.put(j, new ArrayList<Double>());
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
            for(int j = 2; j<4; j++) {
                for(int i=5; i<10; i++){ 
                    if(panel.fixFit[i][j].isSelected()==false) {
                        panel.fixFit[i][j].setSelected(true);
                        System.out.println("REDO THE SEGMENT FITS BEFORE RELEASING THE B-field DEPENDENT PARAMETERS!");
                    }
                }
            }
        }
        boolean[][] fixedPars = new boolean[10][6];
        for(int j = 0; j<6; j++) {
            pars.get(j).clear();
        }
        int npar = 10;
        for(int j = 0; j<6; j++) {
            for(int i=0; i<npar; i++){   
                if(panel.params[i][j].getText().isEmpty()){
                    this.pars.get(j).add(TvstrkdocasFitPars.get(new Coordinate(j)).value(i));
                }
                else { 
                    this.pars.get(j).add(Double.parseDouble(panel.params[i][j].getText()));
                }
            }
        }
        if(!panel.minRange.getText().isEmpty())this.range[0] = Double.parseDouble(panel.minRange.getText());
        else this.range[0] = 0.0;
        if(!panel.maxRange.getText().isEmpty())this.range[1] = Double.parseDouble(panel.maxRange.getText());
        else this.range[1] = 2.0;
        for(int j = 0; j<6; j++) {    
            for(int i=0; i<npar; i++){
                TvstrkdocasFitPars.get(new Coordinate(j)).setValue(i,this.pars.get(j).get(i));
            }
        }
        this._pM.initFitParsToFile();
       
        for(int j = 0; j<6; j++) {
            for(int i=0; i<npar; i++){ 
                if(panel.fixFit[i][j].isSelected()==true)
                    fixedPars[i][j] = true;
            }
        }
        //Don't allow to fit the B>0 profile if they are not filled
        
        this._pM.runFit(fixedPars);
         
            //this._pM.runParamScan(j, fixedPars);
        for(int j = 0; j<6; j++) {    
            for(int p = 0; p<10; p++) {
                panel.pars[p][j] = TvstrkdocasFitPars.get(new Coordinate(j)).value(p);
                if(p!=3) {
                    panel.params[p][j].setText(String.format("%.5f",TvstrkdocasFitPars.get(new Coordinate(j)).value(p)));
                } else {
                    panel.params[p][j].setText(String.format("%.3f",TvstrkdocasFitPars.get(new Coordinate(j)).value(p)));
                }
            }
            panel.fixFit[2][j].setSelected(true);
            panel.fixFit[4][j].setSelected(true);
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
        boolean[][] fixedPars = new boolean[10][6];
        for(int j = 0; j<6; j++) {
            pars.get(j).clear();
        }
        int npar = 10;
        for(int j = 0; j<6; j++) {
            for(int i=0; i<npar; i++){   
                if(panel.params[i][j].getText().isEmpty()){
                    this.pars.get(j).add(TvstrkdocasFitPars.get(new Coordinate(j)).value(i));
                }
                else { 
                    this.pars.get(j).add(Double.parseDouble(panel.params[i][j].getText()));
                }
            }
        }
        if(!panel.minRange.getText().isEmpty())this.range[0] = Double.parseDouble(panel.minRange.getText());
        else this.range[0] = 0.0;
        if(!panel.maxRange.getText().isEmpty())this.range[1] = Double.parseDouble(panel.maxRange.getText());
        else this.range[1] = 2.0;
        for(int j = 0; j<6; j++) {    
            for(int i=0; i<npar; i++){
                TvstrkdocasFitPars.get(new Coordinate(j)).setValue(i,this.pars.get(j).get(i));
            }
        }
        
        for(int j = 0; j<6; j++) {
            for(int i=0; i<npar; i++){ 
                if(panel.fixFit[i][j].isSelected()==true)
                    fixedPars[i][j] = true;
            }
        }
        this._pM.runParamScan(fixedPars);
        initscan=true;
        for(int j = 0; j<6; j++) {    
            for(int p = 0; p<10; p++) {
                panel.pars[p][j] = TvstrkdocasFitPars.get(new Coordinate(j)).value(p);
                if(p!=3) {
                    panel.params[p][j].setText(String.format("%.5f",TvstrkdocasFitPars.get(new Coordinate(j)).value(p)));
                } else {
                    panel.params[p][j].setText(String.format("%.3f",TvstrkdocasFitPars.get(new Coordinate(j)).value(p)));
                }
            }
            panel.fixFit[2][j].setSelected(true);
            panel.fixFit[4][j].setSelected(true);
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
        JCheckBox useBprof;
    	JTextField minRange = new JTextField(5);
	JTextField maxRange = new JTextField(5);
	JTextField[][] params = new JTextField[10][6];
        JCheckBox[][]    fixFit ;
        
        JButton   scanButton = null;
        JButton   fitButton = null;
        JButton   resetButton = null;
        JButton   saveToFileButton = null;
        JButton   resButton = null;
        JButton   reCookButton = null;
        
        String[] parNames = new String[] {"v0", "vmid", "R", "tmax", "distbeta", "delBf", "b1", "b2", "b3", "b4"};
        double[][] pars = new double[10][6];
        
        public CustomPanel2(Map<Coordinate, MnUserParameters> TvstrkdocasFitPars) { 
            
            super(new BorderLayout());
            for(int i = 0; i < 6; i++) {
                for(int p = 0; p<10; p++) {
                    pars[p][i] = TvstrkdocasFitPars.get(new Coordinate(i)).value(p);
                }
            }
            int npar = 10;
            panel = new JPanel(new GridLayout(npar+2, 6));            
            for (int i = 0; i < 6; i++) {
                String SuperLayer = Integer.toString(i + 1);
                panel.add(new JLabel(""));
                panel.add(new JLabel("Super Layer " + SuperLayer));
            }
            panel.add(new JLabel(""));
            fixFit = new JCheckBox[10][6];
            for (int i = 0; i < npar; i++) {  
                JLabel l = new JLabel("      "+parNames[i], JLabel.LEADING);
                panel.add(l);
                for (int j = 0; j < 6; j++) {
                    fixFit[i][j] = new JCheckBox("Fix");
                    // aa is true for parameters "R", "distbeta", "b1", "b2", "b3", and "b4"
                    //boolean aa = i==2 || i==4 || i>5;
                    // aa is true for parameters "b1", "b2", "b3", and "b4"
                    boolean aa = i>5;
                    // bb is true for parameter "tmax" on superlayers 3 and 4 only
                    boolean bb = i==3 && j==6; //(j==2 || j==3);
                    // cc is true for parameter "delBf" on superlayers 1, 2, 5, and 6 only
                    boolean cc = i==5; // && (j==0 || j==1 || j==4 || j==5);
                    if(aa || bb || cc) {
                        fixFit[i][j].setSelected(true);
                    } else {
                        fixFit[i][j].setSelected(false);
                    }
                    
                    params[i][j] = new JTextField(3);
                    if(i!=3) {
                        params[i][j].setText(String.format("%.5f", pars[i][j]));
                    } else {
                        params[i][j].setText(String.format("%.3f", pars[i][j]));
                    }
                    panel.add(params[i][j]);
                    panel.add(fixFit[i][j]);
                }
            }
            panel.add(new JLabel("    Fit range min"));
            minRange.setText(Double.toString(0));
            panel.add(minRange);
            panel.add(new JLabel("    Fit range max"));
            maxRange.setText(Double.toString(2.0));
            panel.add(maxRange);
            useBprof = new JCheckBox(" Use B>0 profiles");
            useBprof.setSelected(false);
            panel.add(useBprof);
                
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

            this.add(panel, BorderLayout.CENTER);
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
