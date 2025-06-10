package org.clas.detector.clas12calibration.dc.fitting;

import org.freehep.math.minuit.FunctionMinimum;
import org.freehep.math.minuit.MnMigrad;
import org.freehep.math.minuit.MnScan;
import org.freehep.math.minuit.MnUserParameters;
import org.jlab.groot.data.IDataSet;
import org.jlab.groot.fitter.FitterFunction;
import org.jlab.groot.math.Func1D;
import org.jlab.groot.math.UserParameter;

/**
 * DataFitter performs a fit on a GROOT function and dataset.
 */
public class DataFitter {

    public static boolean PRINT_OUTPUT = true;

    public DataFitter() {
        // Default constructor
    }

    /**
     * Fits a GROOT Func1D function to a dataset using Minuit-based optimization.
     *
     * @param func    The function to fit (e.g. Gaussian, Polynomial)
     * @param data    The dataset (e.g. histogram)
     * @param options Fitting options: "Q" for quiet, "V" for verbose
     */
    public static void fit(Func1D func, IDataSet data, String options) {
        PRINT_OUTPUT = !options.contains("Q");

        if (PRINT_OUTPUT) {
            System.out.println(">>> Starting fit..."+func.getExpression());
        }

        try {
            FitterFunction fitter = new FitterFunction(func, data, options);
            
            int npars = fitter.getFunction().getNPars();
            MnUserParameters upar = new MnUserParameters();

            for (int i = 0; i < npars; i++) {
                UserParameter par = fitter.getFunction().parameter(i);
                double step = (par.getStep() > 1e-10) ? par.getStep() : 0.0001;
                upar.add(par.name(), par.value(), step);
                upar.setLimits(par.name(), par.min(), par.max());
                if (step < 1e-10) {
                    upar.fix(par.name());
                }
            }

            // Optional scan to improve starting point
            MnScan scanner = new MnScan(fitter, upar);
            FunctionMinimum scanMin = scanner.minimize();
            MnUserParameters scanParams = scanMin.userParameters();
            for (int i = 0; i < npars; i++) {
                UserParameter par = fitter.getFunction().parameter(i);
                par.setValue(scanParams.value(par.name()));
                if(Double.isNaN(scanParams.error(par.name()))) {
                    par.setError(upar.error(par.name()));
                } else {
                    par.setError(scanParams.error(par.name()));
                }
            }
            if (PRINT_OUTPUT) {
                System.out.println(">>> Scanning parameter space...");
                System.out.println(">>> Scan minimum: " + scanMin);
            }

            // Perform main fit
            MnMigrad migrad = new MnMigrad(fitter, upar);
            FunctionMinimum min = migrad.minimize();

            if (min.isValid()) {
                func.setFitValid(true);
                MnUserParameters fitParams = min.userParameters();

                for (int i = 0; i < npars; i++) {
                    UserParameter par = fitter.getFunction().parameter(i);
                    par.setValue(fitParams.value(par.name()));
                    if(Double.isNaN(fitParams.error(par.name()))){
                        par.setError(upar.error(par.name()));
                    } else {
                        par.setError(fitParams.error(par.name()));
                    }
                }

                if (PRINT_OUTPUT) {
                    System.out.println(">>> Fit successful.");
                    System.out.println(">>> Final parameters:");
                    for (int i = 0; i < npars; i++) {
                        UserParameter par = fitter.getFunction().parameter(i);
                        System.out.printf("  - %s: %.6f ± %.6f\n", par.name(), par.value(), par.error());
                    }
                    System.out.println(fitter.getBenchmarkString());
                }
            }
        } catch (Exception e) {
            func.setFitValid(false);
            if (PRINT_OUTPUT) {
                System.err.println("!!! Fit failed: " + e.getClass().getSimpleName() + ": " + e.getMessage());
            }
        }
    }
}
