/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.clas.detector.clas12calibration.dc.mctuning.analysis.wireineff;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Date;
import org.freehep.math.minuit.MnUserParameters;

/**
 *
 * @author veronique
 */

public class OutputManager implements AutoCloseable {
    private final String outputDir;
    private final PrintWriter pw;
    private final File baseFile;

    // Initializes the output manager and sets up the file writer.
    public OutputManager(String outputDir, String baseFilename) throws FileNotFoundException {
        this.outputDir = outputDir;
        this.baseFile = new File(outputDir, baseFilename);
        this.pw = new PrintWriter(baseFile);
        pw.println("#& Superlayer p1 p2 p3 p4 p0");
    }

    // Writes fit parameters to the output file.
    public void writeFitParams(int sl, MnUserParameters params) {
        pw.printf("%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f%n",
                  sl + 1,
                  params.value(1),
                  params.value(2),
                  params.value(3),
                  params.value(4),
                  params.value(0));
    }

    @Override
    public void close() {
        pw.close();
    }

    // Renames the output file with a timestamp and run number.
    public void renameOutput(int runNumber, int iteration) {
        String timestamp = new SimpleDateFormat("MM-dd-yyyy_hh.mm.ss_aa").format(new Date());
        String newName = String.format("wireineff_run%d_time_%s_iteration_%d.txt",
                                       runNumber, timestamp, iteration);
        File newFile = new File(outputDir, newName);
        baseFile.renameTo(newFile);
    }
}
