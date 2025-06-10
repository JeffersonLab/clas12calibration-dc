/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.clas.detector.clas12calibration.dc.t2d;

/**
 *
 * @author ziegler
 */
import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.*;
import java.util.List;

public class ParamAveragingGUI {

    private JFrame frame;
    private List<File> runFiles = new ArrayList<>();
    
    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> new ParamAveragingGUI().createAndShowGUI());
    }

    private void createAndShowGUI() {
        frame = new JFrame("R and distbeta Parameters Averaging Tool");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setSize(400, 200);

        JButton runFilesButton = new JButton("Select Run Data Files");
        JButton runButton = new JButton("Run Average");

        JLabel statusLabel = new JLabel("Select files to begin");

        runFilesButton.addActionListener(e -> selectRunFiles(statusLabel));
        runButton.addActionListener(e -> {
            try {
                runFit(statusLabel);
            } catch (IOException ex) {
                JOptionPane.showMessageDialog(frame, "Error: " + ex.getMessage());
            }
        });

        JPanel panel = new JPanel();
        panel.setLayout(new GridLayout(4, 1));
        panel.add(runFilesButton);
        panel.add(runButton);
        panel.add(statusLabel);

        frame.getContentPane().add(panel);
        frame.setVisible(true);
    }

    
    private void selectRunFiles(JLabel statusLabel) {
        JFileChooser fileChooser = new JFileChooser();
        fileChooser.setDialogTitle("Select Run Data Files");
        fileChooser.setMultiSelectionEnabled(true);
        int result = fileChooser.showOpenDialog(frame);
        if (result == JFileChooser.APPROVE_OPTION) {
            runFiles = Arrays.asList(fileChooser.getSelectedFiles());
            statusLabel.setText(runFiles.size() + " run files selected.");
        }
    }
   
    private void runFit(JLabel statusLabel) throws IOException {
        if (runFiles.isEmpty()) {
            JOptionPane.showMessageDialog(frame, "Please select both pressure file and run data files.");
            return;
        }
        
        
        double[] R = new double[3];
        double[] dB = new double[3];
        for (File runFile : runFiles) {
            System.out.println(runFile.getName());
            try (BufferedReader reader = new BufferedReader(new FileReader(runFile))) {
                String line;
                while ((line = reader.readLine()) != null) {
                    if (line.startsWith("#") || line.trim().isEmpty()) continue;

                    String[] tokens = line.trim().split("\\s+");
                    int sector = Integer.parseInt(tokens[0]);
                    int superlayer = Integer.parseInt(tokens[1]);
                    double r = Double.parseDouble(tokens[30]);
                    double distbeta = Double.parseDouble(tokens[12]);
                    int regionIdx = (int) ((superlayer-1)/2);
                    
                    if(sector==1) {
                        R[regionIdx]+=r; System.out.println("superlayer "+superlayer); System.out.println("R "+r);
                        dB[regionIdx]+=distbeta; System.out.println("dB "+distbeta);
                    }
                }
            }
        }
        
        df.setRoundingMode(RoundingMode.UP);
        double N = runFiles.size()*2; 
       
        String st = "";
        for(int i = 0; i < 3; i++) {
            st+="R_"+(i+1)+" = "+df.format(R[i]/N)+"\n";
            st+="distbeta_"+(i+1)+" = "+df.format(dB[i]/N)+"\n";
        }
        
        // Output to console
        System.out.println(st);
        

        statusLabel.setText("Task complete. See console output.");
    }
      private static final DecimalFormat df = new DecimalFormat("0.000");

      

}
