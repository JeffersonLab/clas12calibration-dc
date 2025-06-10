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
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JFrame;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.JToggleButton;
import javax.swing.plaf.metal.MetalButtonUI;
import org.clas.detector.clas12calibration.dc.caltdccuts.TDCCuts;
import static org.clas.detector.clas12calibration.dc.caltdccuts.TDCCuts.cutParams;
import org.clas.detector.clas12calibration.viewer.TDCViewer;

public class TDCParamsPanel {
    
    
    private Map<Coordinate, ArrayList<Double>> pars     = new HashMap<>();
    private JFrame            frame     = new JFrame();
    public CustomPanel2      panel     = null;            
    private TDCCuts _pM;
    private TDCViewer _vM;
    public TDCParamsPanel(TDCCuts pM, TDCViewer vM) {
        this._vM = vM;
        this._pM = pM;
        //init pars container
        for(int r = 0; r<6; r++) {
            for(int w = 0; w<3; w++) {
                pars.put(new Coordinate(r,w), new ArrayList<>());
            }
        }
    }
    
    public void openFitPanel(String title){
        System.out.println("TRYING>>>>>>>>>>>>>");
        frame    = new JFrame();
        panel = new CustomPanel2();
        frame.setSize(350, 300); 
        frame.setTitle(title);
        frame.add(panel);
        frame.pack();
        frame.setVisible(true);
    }
    public void setGreenFitButton() {
        panel.calcParsButton.setBackground(Color.GREEN);
        panel.calcParsButton.setVisible(true);
        panel.repaint();
    }
    public void setRedFitButton() {
        panel.calcParsButton.setBackground(Color.RED);
        panel.calcParsButton.setVisible(true);
        panel.repaint();
    }
    public boolean fitted = false;
    public boolean initscan = false;
    public void computeBoundPars() {
        if (this._pM.analysisDone == false) return;
        System.out.println(" COMPUTING BOUNDS PARAMETERS ");

        // Clear the current map before re-filling it
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 3; j++) {
                Coordinate coord = new Coordinate(i, j);
                pars.get(coord).clear(); // clear previous entries

                String text = panel.params[i][j].getText();
                if (!text.isEmpty()) {
                    try {
                        double value = Double.parseDouble(text);
                        pars.get(coord).add(value);
                        System.out.println("Parsed value: " + value + " for (" + i + "," + j + ")");
                    } catch (NumberFormatException ex) {
                        System.err.println("Invalid number in (" + i + "," + j + "): " + text);
                    }
                }
            }
        }

        this._pM.computeBounds(this.pars);
        double r2max1rnd = new BigDecimal(cutParams.get(1)[1])
                    .setScale(1, RoundingMode.HALF_UP)
                    .doubleValue();
        double r2LC1rnd = new BigDecimal(cutParams.get(1)[2])
                    .setScale(1, RoundingMode.HALF_UP)
                    .doubleValue(); 
        double r2max56rnd = new BigDecimal(cutParams.get(1)[4])
                    .setScale(1, RoundingMode.HALF_UP)
                    .doubleValue();
        double r2LC56rnd = new BigDecimal(cutParams.get(1)[5])
                    .setScale(1, RoundingMode.HALF_UP)
                    .doubleValue(); 
        
        List<String[]> outputTable = new ArrayList<>();
        outputTable.add(new String[]{"0", "1", "0", Double.toString(cutParams.get(0)[0]), Double.toString(cutParams.get(0)[1]), "0"});
        outputTable.add(new String[]{"0", "2", "1", Double.toString(cutParams.get(1)[0]), Double.toString(r2max1rnd), Double.toString(r2LC1rnd)});
        outputTable.add(new String[]{"0", "2", "56",Double.toString(cutParams.get(1)[3]), Double.toString(r2max56rnd), Double.toString(r2LC56rnd)});
        outputTable.add(new String[]{"0", "3", "0", Double.toString(cutParams.get(2)[0]), Double.toString(cutParams.get(2)[1]), "0"});
        
        
        this._pM.updateOutputArea(outputTable, this);
        
    }
    DateFormat df = new SimpleDateFormat("MM-dd-yyyy_hh.mm.ss_aa");
    public void savePars() {
        double r2max1rnd = new BigDecimal(cutParams.get(1)[1])
                    .setScale(1, RoundingMode.HALF_UP)
                    .doubleValue();
        double r2LC1rnd = new BigDecimal(cutParams.get(1)[2])
                    .setScale(1, RoundingMode.HALF_UP)
                    .doubleValue(); 
        double r2max56rnd = new BigDecimal(cutParams.get(1)[4])
                    .setScale(1, RoundingMode.HALF_UP)
                    .doubleValue();
        double r2LC56rnd = new BigDecimal(cutParams.get(1)[5])
                    .setScale(1, RoundingMode.HALF_UP)
                    .doubleValue(); 
        
        List<String[]> outputTable = new ArrayList<>();
        outputTable.add(new String[]{"0", "1", "0", Double.toString(cutParams.get(0)[0]), Double.toString(cutParams.get(0)[1]), "0"});
        outputTable.add(new String[]{"0", "2", "1", Double.toString(cutParams.get(1)[0]), Double.toString(r2max1rnd), Double.toString(r2LC1rnd)});
        outputTable.add(new String[]{"0", "2", "56",Double.toString(cutParams.get(1)[3]), Double.toString(r2max56rnd), Double.toString(r2LC56rnd)});
        outputTable.add(new String[]{"0", "3", "0", Double.toString(cutParams.get(2)[0]), Double.toString(cutParams.get(2)[1]), "0"});
        
        String filePath = "Files/ccdb_T00Corr_run" + _pM.runNumber + "time_" 
                + df.format(new Date())  + ".txt";
        try (PrintWriter writer = new PrintWriter(new FileWriter(filePath))) {
            // Write the header
            writer.printf("%-8s %-8s %-10s %-8s %-8s %-12s%n", 
                          "& sector", "layer", "component", "MinEdge", "MaxEdge", "LinearCoeff");

            // Write each row
            for (String[] row : outputTable) {
                writer.printf("%-8s %-8s %-10s %-8s %-8s %-12s%n",
                              row[0], row[1], row[2], row[3], row[4], row[5]);
            }
            writer.flush();
            System.out.println("Output successfully written to: " + filePath);
        } catch (IOException ex) {
            Logger.getLogger(TDCParamsPanel.class.getName()).log(Level.SEVERE, null, ex);
        }    
    }
    
   
    public final class CustomPanel2 extends JPanel {
        JLabel label;
        JPanel panel;
        
        JTabbedPane tpane;
        public JTextArea outputArea;
	JTextField[][] params = new JTextField[6][3];
        
        JButton   saveToFileButton = null;
        JButton   calcParsButton = null;
        JToggleButton logScaleToggle = null;
        
        String[] parNames = new String[] {"R1 Y lower bound", "R1 Y upper bound", 
                                          "R2 Y lower bound", "R2 Y upper bound",
                                          "R3 Y lower bound", "R3 Y upper bound"};
        double[][] pars = new double[6][3];
        
        
        private JPanel addPanel() {
            JPanel spanel;            
            spanel = new JPanel(new GridLayout(7, 4));
            spanel.add(new JLabel("coordinates"));
            spanel.add(new JLabel("wire 1"));
            //spanel.add(new JLabel(""));
            spanel.add(new JLabel("wire 56"));
            //spanel.add(new JLabel(""));
            spanel.add(new JLabel("wire 112"));
            
            //spanel.add(new JLabel(""));
            for (int i = 0; i < parNames.length; i++) {  
                JLabel l = new JLabel("      "+parNames[i], JLabel.LEADING);
                spanel.add(l);
                
                for(int j = 0; j <3 ; j++) {
                    params[i][j] = new JTextField(1);
                    if(j==1 && (i<3 || i>3)) {
                        params[i][j].setBackground(Color.DARK_GRAY);
                    }
                    params[i][j].addActionListener(_vM);
                    spanel.add(params[i][j]);
                }
            }
            spanel.setVisible(true);
            return spanel;
        }
        
        public CustomPanel2() { 
            super(new BorderLayout());
                for(int i = 0; i < 6; i++) {
                    for(int p = 0; p<3; p++) {
                        pars[i][p] = 0;
                    
                }
            }
            panel = new JPanel(new GridLayout(1, 1));
            panel.setBackground(Color.white);
            tpane = new JTabbedPane();
            
            panel.add(addPanel());
            
            
            outputArea = new JTextArea(8, 40);
            outputArea.setEditable(false);
            outputArea.setFont(new Font("Monospaced", Font.PLAIN, 12));
            JScrollPane scrollPane = new JScrollPane(outputArea);
            scrollPane.setBorder(BorderFactory.createTitledBorder("ccdb parameters estimate"));

            this.add(scrollPane, BorderLayout.CENTER);
            
            JPanel buttonsPanel = new JPanel();

            saveToFileButton = new JButton("SAVE PARAMETERS");
            saveToFileButton.setUI(new MetalButtonUI());
            saveToFileButton.setBackground(Color.PINK);
            saveToFileButton.setContentAreaFilled(true);
            saveToFileButton.setOpaque(true);
            saveToFileButton.setFont(bBold);
            saveToFileButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    savePars();
                    return;
                }

                
            });
            
            
            calcParsButton = new JButton("COMPUTE BOUNDS PARAMS");
            calcParsButton.setUI(new MetalButtonUI());
            calcParsButton.setBackground(Color.YELLOW);
            calcParsButton.setContentAreaFilled(false);
            calcParsButton.setOpaque(true);
            calcParsButton.setFont(bBold);
            calcParsButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    computeBoundPars();
                    return;
                }
            });
            
            logScaleToggle = new JToggleButton("Log Scale: ON");
            logScaleToggle.setFont(bBold);
            logScaleToggle.setBackground(Color.CYAN);
            logScaleToggle.setFocusPainted(false);
            logScaleToggle.addItemListener(new ItemListener() {
                public void itemStateChanged(ItemEvent e) {
                    boolean isLog = e.getStateChange() == ItemEvent.SELECTED;
                    logScaleToggle.setText("Log Scale: " + (isLog ? "OFF" : "ON"));
                    handleLogScaleToggle(isLog);
                }

                private void handleLogScaleToggle(boolean log) {
                    _pM.isLog=log;
                }
            });
            
            buttonsPanel.add(calcParsButton);
            buttonsPanel.add(logScaleToggle);
            buttonsPanel.add(saveToFileButton);
            
            this.add(panel, BorderLayout.PAGE_START);
            this.add(buttonsPanel, BorderLayout.PAGE_END);
            
            label = new JLabel("Click the \"Show it!\" button"
                           + " to bring up the selected dialog.",
                           JLabel.CENTER);
            
        }
        private Font bBold = new Font("Arial", Font.BOLD, 13);
        void setLabel(String newText) {
            label.setText(newText);
        }
    }
}
