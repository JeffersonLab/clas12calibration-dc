/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.clas.detector.clas12calibration.dc.t2d;

/**
 *
 * @author ziegler
 */
import java.io.*;
import java.util.*;
import org.clas.detector.clas12calibration.dc.calt2d.T2DCalib;
import static org.clas.detector.clas12calibration.viewer.T2DViewer.voice;

public class ParameterParser {
    public static double R_1;
    public static double distbeta_1;
    public static double R_2;
    public static double distbeta_2;
    public static double R_3;
    public static double distbeta_3;
    
    public static boolean parseRDPars() {
        String fileName = "parameters.txt"; 
        Map<String, Double> parameters = new HashMap<>();

        File file = new File(fileName);
        if (!file.exists()) {
            System.out.println("File does not exist: " + fileName);
            return false;
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
                    parameters.put(key, value);
                }
            }

            // Output values
            R_1 =  parameters.getOrDefault("R_1", null);
            distbeta_1 = parameters.getOrDefault("distbeta_1", null);
            R_2 =  parameters.getOrDefault("R_2", null);
            distbeta_2 =  parameters.getOrDefault("distbeta_2", null);
            R_3 =  parameters.getOrDefault("R_3", null);
            distbeta_3 =  parameters.getOrDefault("distbeta_3", null);
            
            System.out.println("PARAMETERS READ FROM FILE:");
            System.out.printf("%-15s | %-15s%n", "Parameter", "Value");
            System.out.println("-----------------+-----------------");
            System.out.printf("%-15s | %-15s%n", "R_1", R_1);
            System.out.printf("%-15s | %-15s%n", "distbeta_1", distbeta_1);
            System.out.printf("%-15s | %-15s%n", "R_2", R_2);
            System.out.printf("%-15s | %-15s%n", "distbeta_2", distbeta_2);
            System.out.printf("%-15s | %-15s%n", "R_3", R_3);
            System.out.printf("%-15s | %-15s%n", "distbeta_3", distbeta_3);
            
            if(T2DCalib.vocal) voice.speak("SCAN PARAMETERS READ FROM FILE.");

        } catch (IOException | NumberFormatException e) {
            System.err.println("Error reading file: " + e.getMessage());
            return false;
        }
        return true;
    }
}

