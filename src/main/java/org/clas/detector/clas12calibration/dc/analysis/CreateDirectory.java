/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.clas.detector.clas12calibration.dc.analysis;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
/**
 *
 * @author ziegler
 */


public class CreateDirectory {
    public static void create(String dir) {
        String dirPath = "./"+dir;
        Path path = Paths.get(dirPath);

        // Check if the directory exists
        if (!Files.exists(path)) {
            try {
                // Create the directory
                Files.createDirectory(path);
                System.out.println("Directory created: " + path.toAbsolutePath());
            } catch (IOException e) {
                System.err.println("Failed to create directory: " + e.getMessage());
                e.printStackTrace();
            }
        } else {
            System.out.println("Directory already exists: " + path.toAbsolutePath());
        }
    }
}