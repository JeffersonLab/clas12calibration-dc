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
import java.util.*;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.math.F1D;
import org.jlab.groot.ui.TCanvas;

public class PressureFitGUI {

    private JFrame frame;
    private File pressureFile;
    private List<File> runFiles;
    
    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> new PressureFitGUI().createAndShowGUI());
    }

    private void createAndShowGUI() {
        frame = new JFrame("Pressure Fitting Tool");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setSize(400, 200);

        JButton pressureButton = new JButton("Select Pressure File");
        JButton runFilesButton = new JButton("Select Run Data Files");
        JButton RDbButton = new JButton("Select R, Distbeta Parameter Files");
        JButton runButton = new JButton("Run Fit");

        JLabel statusLabel = new JLabel("Select files to begin");

        pressureButton.addActionListener(e -> selectPressureFile(statusLabel));
        runFilesButton.addActionListener(e -> selectRunFiles(statusLabel));
        RDbButton.addActionListener(e -> {
            try {
                selectRDbParamsFile(statusLabel);
            } catch (IOException ex) {
                Logger.getLogger(PressureFitGUI.class.getName()).log(Level.SEVERE, null, ex);
            }
        });
        runButton.addActionListener(e -> {
            try {
                runFit(statusLabel);
            } catch (IOException ex) {
                JOptionPane.showMessageDialog(frame, "Error: " + ex.getMessage());
            }
        });

        JPanel panel = new JPanel();
        panel.setLayout(new GridLayout(4, 1));
        panel.add(pressureButton);
        panel.add(runFilesButton);
       panel.add( RDbButton);
        panel.add(runButton);
        panel.add(statusLabel);

        frame.getContentPane().add(panel);
        frame.setVisible(true);
    }

    private void selectPressureFile(JLabel statusLabel) {
        JFileChooser fileChooser = new JFileChooser();
        fileChooser.setDialogTitle("Select Pressure File");
        int result = fileChooser.showOpenDialog(frame);
        if (result == JFileChooser.APPROVE_OPTION) {
            pressureFile = fileChooser.getSelectedFile();
            statusLabel.setText("Pressure file selected: " + pressureFile.getName());
        }
    }

    private void selectRDbParamsFile(JLabel statusLabel) throws FileNotFoundException, IOException {
        File file=null;
        JFileChooser fileChooser = new JFileChooser();
        fileChooser.setDialogTitle("Select R, Db Parameter File");
        int result = fileChooser.showOpenDialog(frame);
        if (result == JFileChooser.APPROVE_OPTION) {
            file = fileChooser.getSelectedFile();
            statusLabel.setText("parameter file selected: " + file.getName());
        }
        try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
            String line;

            while ((line = reader.readLine()) != null) {
                line = line.trim();

                // Skip empty lines
                if (line.isEmpty()) continue;

                String[] parts = line.split("=");
                if (parts.length == 2) {
                    String key = parts[0].trim(); 
                    double value = Double.parseDouble(parts[1].trim());
                    String[] parts2 = key.split("_");
                    String varName = parts2[0].trim();
                    int region = Integer.parseInt(parts2[1].trim());
                    if(varName.equals("R")) {
                        R[region-1]=(float)value; System.out.println("R in region "+region +" = "+R[region-1]);
                    } else {
                        Db[region-1]=(float)value;
                    }
                }
            }
        }
    }
    float[] R = new float[3];
    float[] Db = new float[3];
    double[] vo0 = new double[6];
    double[] vo1 = new double[6];
    double[] vmid0 = new double[6];
    double[] vmid1 = new double[6];
    double[] tmax0 = new double[6];
    double[] tmax1 = new double[6];
    double[] delta_bfield0 = new double[6];
    double[] delta_bfield1 = new double[6];
    double[] b1 = new double[]{1,1,1,1,1,1};
    double[] b2 = new double[]{-2,-2,-2,-2,-2,-2};
    double[] b3 = new double[]{11.5,11.5,11.5,11.5,11.5,11.5};
    double[] b4 = new double[]{-6.5,-6.5,-6.5,-6.5,-6.5,-6.5};
    
    int markerStyle=1;
    int fitColor = 1;
    private void selectRunFiles(JLabel statusLabel) {
        markerStyle++;
        fitColor++;
        runFiles = new ArrayList<>();
        JFileChooser fileChooser = new JFileChooser();
        fileChooser.setDialogTitle("Select Run Data Files");
        fileChooser.setMultiSelectionEnabled(true);
        int result = fileChooser.showOpenDialog(frame);
        if (result == JFileChooser.APPROVE_OPTION) {
            runFiles = Arrays.asList(fileChooser.getSelectedFiles());
            statusLabel.setText(runFiles.size() + " run files selected.");
        }
    }
    F1D f1 = new F1D("f1","[a0]+[a1]*x", -20, 20);
    double min = Double.MAX_VALUE;
    double max = Double.MIN_VALUE;
    int nSec = 1;
    double p0=410;
    private void runFit(JLabel statusLabel) throws IOException {
        if (pressureFile == null || runFiles.isEmpty()) {
            JOptionPane.showMessageDialog(frame, "Please select both pressure file and run data files.");
            return;
        }
        
        Map<Integer, Double> runToPressure = readPressureFile(pressureFile);
        for (double pressure : runToPressure.values()) {
            if (pressure < min) min = pressure-p0;
            if (pressure > max) max = pressure-p0;
        }
        
        Map<Coordinate, Map<String, GraphErrors>> groupedData = new HashMap<>();
        double N = 0;
        double v0A = 0;
        double vmidA = 0;
        double tmaxA = 0;
        double dBA = 0;
        double[][] maxAxy = new double[][]{
                {Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY},
                {Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY},
                {Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY},
                {Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY}
                };
        double[][] minAxy = new double[][]{
                {Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY},
                {Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY},
                {Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY},
                {Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY}
                };
        double[][] miny = new double[4][6];
        for (File runFile : runFiles) {
            String fileName = runFile.getName(); System.out.println("READING "+runFile);
            int runNumber = extractRunNumber(fileName); 
            if (!runToPressure.containsKey(runNumber)) {
                System.err.println("No pressure for run: " + runNumber);
                continue;
            }
            double pressure = runToPressure.get(runNumber)-p0;
            System.out.println("pressure "+pressure+" for run number "+runNumber);
            
            try (BufferedReader reader = new BufferedReader(new FileReader(runFile))) {
                String line;
                while ((line = reader.readLine()) != null) {
                    if (line.startsWith("#") || line.trim().isEmpty()) continue;

                    String[] tokens = line.trim().split("\\s+");
                    
                    int sector = Integer.parseInt(tokens[0]);
                    if(sector==1) {
                        N+=1.0;
                        int superlayer = Integer.parseInt(tokens[1]);
                        double v0 = Double.parseDouble(tokens[3]);
                        double v0E = Double.parseDouble(tokens[5])*100;
                        double tmax = Double.parseDouble(tokens[6]);
                        double tmaxE = Double.parseDouble(tokens[8])*100;
                        double vmid = Double.parseDouble(tokens[9]);
                        double vmidE = Double.parseDouble(tokens[11])*100;  
                        if(vmidE==0) vmidE=0.001;
                        double delta_bf = Double.parseDouble(tokens[12]);
                        double delta_bfE = Double.parseDouble(tokens[14])*100;
                        if(superlayer<3 || superlayer>4) delta_bfE=1.e-9;
                        v0A += v0;
                        vmidA += vmid;
                        tmaxA += tmax;
                        dBA += delta_bf;
                        
                        if(v0+10*v0E>maxAxy[0][superlayer-1]) {
                            maxAxy[0][superlayer-1]=v0+10*v0E;
                        }
                        if(v0-10*v0E<minAxy[0][superlayer-1]) {
                            minAxy[0][superlayer-1]=v0-10*v0E;
                        }
                        if(tmax+10*tmaxE>maxAxy[1][superlayer-1]) {
                            maxAxy[1][superlayer-1]=tmax+10*tmaxE;
                        }
                        if(tmax-10*tmaxE<minAxy[1][superlayer-1]) {
                            minAxy[1][superlayer-1]=tmax-10*tmaxE;
                        }
                        if(vmid+10*vmidE>maxAxy[2][superlayer-1]) {
                            maxAxy[2][superlayer-1]=vmid+10*vmidE;
                        }
                        if(vmid-10*vmidE<minAxy[2][superlayer-1]) {
                            minAxy[2][superlayer-1]=vmid-10*vmidE;
                        }
                        if(delta_bf+10*delta_bfE>maxAxy[3][superlayer-1]) {
                            maxAxy[3][superlayer-1]=delta_bf+10*delta_bfE;
                        }
                        if(delta_bf-10*delta_bfE<minAxy[3][superlayer-1]) {
                            minAxy[3][superlayer-1]=delta_bf-10*delta_bfE;
                        }
                        
                        Coordinate coord = new Coordinate(sector-1, superlayer-1);
                        groupedData.computeIfAbsent(coord, k -> new HashMap<>())
                                   .computeIfAbsent("v0", k -> new GraphErrors())
                                   .addPoint(pressure, v0, 0, v0E);
                        groupedData.computeIfAbsent(coord, k -> new HashMap<>())
                                   .computeIfAbsent("tmax", k -> new GraphErrors())
                                   .addPoint(pressure, tmax, 0, tmaxE);
                        groupedData.computeIfAbsent(coord, k -> new HashMap<>())
                                   .computeIfAbsent("vmid", k -> new GraphErrors())
                                   .addPoint(pressure, vmid, 0, vmidE);
                        //if(superlayer==3 || superlayer==4) {
                            groupedData.computeIfAbsent(coord, k -> new HashMap<>())
                                   .computeIfAbsent("delta_bf", k -> new GraphErrors())
                                   .addPoint(pressure, delta_bf, 0, delta_bfE);
                        //}
                    }
                }
            }
        }
        v0A /= N;
        vmidA /= N;
        tmaxA /= N;
        dBA /= N;
        //Fit
        double[] iPars = new double[] {v0A, tmaxA, vmidA, dBA};
        String[] vars = new String[] {"v0", "tmax", "vmid", "delta_bf"};
        String output = new String();
        
        String ccdbFile = runFiles.getFirst().getName();
        PrintWriter pw = new PrintWriter("PressureDepT2D_"+ccdbFile);
        pw.printf("#& sector superlayer component v0_a0 v0_a1 v0_a2 vmid_a0 vmid_a1 vmid_a2 tmax_a0 tmax_a1 tmax_a2 distbeta_a0 distbeta_a1 distbeta_a2 delta_bfield_a0 delta_bfield_a1 delta_bfield_a2 b1_a0 b1_a1 b1_a2 b2_a0 b2_a1 b2_a2 b3_a0 b3_a1 b3_a2 b4_a0 b4_a1 b4_a2 c1_a0 c1_a1\n");
           
        for(int v = 0; v < vars.length; v++) {
            TCanvas can = new TCanvas(vars[v], 1500, 300);
            can.divide(6, nSec);
            int cindex=0;
            for(int s = 0; s < nSec; s++) {
                for(int l =0; l < 6; l++) {
                    if(v==3 && l<2) continue;
                    if(v==3 && l>3) continue;
                    GraphErrors g = groupedData.get(new Coordinate(s, l )).get(vars[v]);
                    g.setMarkerStyle(this.markerStyle);
                    f1 = new F1D("f1","[a0]+[a1]*x", -20, 20);
                    f1.setParameter(0, iPars[v]);
                    f1.setParameter(1, 0);
                    f1.setOptStat(0111);
                    DataFitter.fit(f1, g, "V");
                    f1.setLineColor(this.fitColor);
                    String output2="sector "+(s+1)+" superlayer "+(l+1)+": ";
                    output2+=vars[v]+"_0 = ";
                    output2+=(float)f1.getParameter(0)+"+/- "+(float)f1.parameter(0).error()+", ";
                    output2+=vars[v]+"_1 = ";
                    output2+=(float)f1.getParameter(1)+"+/- "+(float)f1.parameter(1).error()+"\n ";
                    System.out.println(output2);
                    output+=output2;
                    if(v==0) {
                        vo0[l]=f1.getParameter(0);
                        vo1[l]=f1.getParameter(1);
                    } 
                    if(v==2) {
                        vmid0[l]=f1.getParameter(0);
                        vmid1[l]=f1.getParameter(1);
                    } 
                    if(v==1) {
                        tmax0[l]=f1.getParameter(0);
                        tmax1[l]=f1.getParameter(1);
                    } 
                    if(v==3) {
                        delta_bfield0[l]=f1.getParameter(0);
                        delta_bfield1[l]=f1.getParameter(1);
                    } 
                    
                    can.cd(cindex);
                    
                    can.getPad().setAxisRange(-25, 25, minAxy[v][l], maxAxy[v][l]);
                    can.draw(g);
                    
                    cindex++;
                }        
            }
        }
        for(int s = 0; s < 6; s++) {
            for(int l =0; l < 6; l++) {
                pw.printf("%d\t %d\t %d\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %.6f\n",
                    (s+1), (l+1), 0,
                    vo0[l],
                    vo1[l],
                    0.0,
                    vmid0[l],
                    vmid1[l],
                    0.0,
                    tmax0[l],
                    tmax1[l],
                    0.0,
                    this.Db[(int)(l/2)],
                    0.0,
                    0.0,
                    delta_bfield0[l],
                    delta_bfield1[l],
                    0.0,
                    b1[l],
                    0.0,
                    0.0,
                    b2[l],
                    0.0,
                    0.0,
                    b3[l],
                    0.0,
                    0.0,
                    b4[l],
                    0.0,
                    0.0,
                    this.R[(int)(l/2)],
                    0.0,
                    0.0);
            }
        }
        pw.close();
        // Output to console
        System.out.println(output);
        

        statusLabel.setText("Fitting complete. See console output.");
    }
    
    private int extractRunNumber(String fileName) {
        // Match "run" followed by digits, case-insensitive
        Matcher matcher = Pattern.compile("run(\\d+)", Pattern.CASE_INSENSITIVE).matcher(fileName);
        if (matcher.find()) {
            return Integer.parseInt(matcher.group(1));
        } else {
            throw new IllegalArgumentException("Run number not found in file name: " + fileName);
        }
    }


    private Map<Integer, Double> readPressureFile(File file) throws IOException {
        Map<Integer, Double> map = new HashMap<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#") || line.trim().isEmpty()) continue;
                String[] tokens = line.trim().split("\\s+");
                System.out.println(tokens[0]+"  "+tokens[1]);
                if (tokens.length >= 2) {
                    map.put(Integer.parseInt(tokens[0]), Double.parseDouble(tokens[1]));
                }
            }
        }
        return map;
    }

}
