/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.clas.detector.clas12calibration.dc.mctuning.analysis.wireineff;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.clas.detector.clas12calibration.dc.mctuning.viewer.WireIneffAnalViewer;
import org.jlab.detector.geant4.v2.DCGeant4Factory;
import org.jlab.rec.dc.cluster.Cluster;
import org.jlab.rec.dc.cluster.ClusterFitter;
import org.jlab.rec.dc.cluster.FittedCluster;
import org.jlab.rec.dc.hit.FittedHit;
import org.jlab.rec.dc.segment.Segment;
//import org.jlab.rec.dc.trajectory.SegmentTrajectory;
import org.clas.detector.clas12calibration.dc.calt2d.HitUtility;
import static org.clas.detector.clas12calibration.dc.mctuning.analysis.wireineff.EventProcessor.swim;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.rec.dc.Constants;
/**
 *
 * @author veronique
 */



public class SegmentAnalysis {

    private static FittedHit getHit(DataBank bnkHits, int i) {
  
        FittedHit hit = HitUtility.getHit(bnkHits, i,WireIneffAnalViewer.dcDetector);
        if (bnkHits.getInt("trkID", i) >0) {  
            return hit;
        }
        return null;
    }
    
    private static List<FittedCluster> recomposeClusters(List<FittedHit> fhits) {
        ClusterFitter cf = new ClusterFitter();
        Map<Integer, List<FittedHit>> grpHits = new HashMap<>();
        List<FittedCluster> clusters = new ArrayList<>();

        for (FittedHit hit : fhits) {
            int clusterId = hit.get_AssociatedClusterID();
            int trackId   = hit.get_AssociatedHBTrackID();

            // Skip if no valid cluster or track
            if (clusterId == -1 || trackId == -1) continue;

            int index = trackId * 10000 + clusterId;
            grpHits.computeIfAbsent(index, k -> new ArrayList<>()).add(hit);
        }

        for (Integer entry : grpHits.keySet()) {
            List<FittedHit> hits = grpHits.get(entry);
            if (hits.size() <= 3) continue;

            FittedHit seed = hits.get(0);
            Cluster cluster = new Cluster(seed.get_Sector(), seed.get_Superlayer(), seed.get_AssociatedClusterID());
            FittedCluster fcluster = new FittedCluster(cluster);
            fcluster.addAll(hits);

            fcluster.forEach(h -> h.set_TrkgStatus(1));

            cf.SetFitArray(fcluster, "TSC");
            cf.Fit(fcluster, true);
            cf.SetResidualDerivedParams(fcluster, true, false, WireIneffAnalViewer.dcDetector);
            cf.Fit(fcluster, false);
            cf.SetSegmentLineParameters(fcluster.get(0).get_Z(), fcluster);

            clusters.add(fcluster);
        }
        return clusters;
    }
    
    public static List<Segment> getSegments(int run, List<FittedCluster> allClusters, DataEvent event, DCGeant4Factory DcDetector) {
        List<Segment> segList = new ArrayList<>();

        for (FittedCluster fClus : allClusters) {
            if (fClus.size() > 12 || fClus.get_TrkgStatus() == -1) continue;

            Segment seg = new Segment(fClus);
            seg.set_fitPlane(DcDetector);

            boolean hasTDCBank = event.hasBank("DC::tdc"); 
            Map<Integer, List<int[]>> map = getWireMap(event, hasTDCBank);

            getLayerEfficiencies(run, seg, map, DcDetector);

            if (seg.get_Status() != -1) {
                double sumRes = 0, sumTime = 0;
                for (FittedHit h : seg) {
                    sumRes += h.get_TimeResidual();
                    sumTime += h.get_Time();
                }
                seg.set_ResiSum(sumRes);
                seg.set_TimeSum(sumTime);
                segList.add(seg);
            }
        }
        return segList;
    }
    public static double[] cellSizes=new double[6];
    private static WireTrajectoryResult getWireOnTrajectory(Segment seg, int layer, double trkX, DCGeant4Factory DcDetector) {
        int sector = seg.get_Sector();
        int superlayer = seg.get_Superlayer();
        //double cellSize = 2.0*Constants.getInstance().wpdist[superlayer-1];
        
        double x1 = DcDetector.getWireMidpoint(sector - 1, superlayer - 1, layer - 1, 1).x;
        double x0 = DcDetector.getWireMidpoint(sector - 1, superlayer - 1, layer - 1, 0).x;
        double deltaX = Math.abs(x1 - x0);
        double cellSize = (deltaX / 2);
        if(cellSizes[superlayer-1]==0)
            cellSizes[superlayer-1]=cellSize;
        
        double xFirstCell = x0;
        int nearestWire = (int) Math.ceil((trkX - xFirstCell + deltaX / 2.0) / deltaX);
        nearestWire      = adjustedWire(nearestWire); //make sure it's between 1 and 112
        
        int nearestWireM = adjustedWire(nearestWire - 1);
        int nearestWireP = adjustedWire(nearestWire + 1);

        double xn  = DcDetector.getWireMidpoint(sector - 1, superlayer - 1, layer - 1, nearestWire - 1).x;
        double ndoca  = getNormalizedDoca(trkX, xn, cellSize);

        double xnM = DcDetector.getWireMidpoint(sector - 1, superlayer - 1, layer - 1, nearestWireM - 1).x;
        double ndocaM  = getNormalizedDoca(trkX, xnM, cellSize);

        double xnP = DcDetector.getWireMidpoint(sector - 1, superlayer - 1, layer - 1, nearestWireP - 1).x;
        double ndocaP  = getNormalizedDoca(trkX, xnP, cellSize);

        int[] wires;
        double[] docas;
        if (Math.abs(trkX - xnP) < Math.abs(trkX - xnM)) {
            wires = new int[]{nearestWire, nearestWireP};
            docas = new double[]{ndoca, ndocaP};
        } else {
            wires = new int[]{nearestWireM, nearestWire};
            docas = new double[]{ndocaM, ndoca};
        }
        return new WireTrajectoryResult(wires, docas);
    }

    private static int adjustedWire(int wire) {
        return Math.max(1, Math.min(wire, 112)); // clamps value into [1,112]
    }

    private static Map<Integer, List<int[]>> getWireMap(DataEvent event, boolean hasTDCBank) {
        Map<Integer, List<int[]>> map = new HashMap<>();
        DataBank bank;
        String wireField;
        boolean layerIs1to36;

        if (hasTDCBank && event.hasBank("DC::tdc")) {
            bank = event.getBank("DC::tdc");
            wireField = "component";
            layerIs1to36 = true;
        } else if (!hasTDCBank && event.hasBank("TimeBasedTrkg::TBHits")) {
            bank = event.getBank("TimeBasedTrkg::TBHits");
            wireField = "wire";
            layerIs1to36 = false;
        } else {
            return map; // No relevant banks
        }

        for (int i = 0; i < bank.rows(); i++) {
            int sector = bank.getByte("sector", i);
            int layer  = bank.getByte("layer", i);
            int wire   = bank.getShort(wireField, i);

            int superlayer = layerIs1to36 ? (layer - 1) / 6 + 1 : bank.getByte("superlayer", i);
            layer = layerIs1to36 ? (layer - 1) % 6 + 1 : bank.getByte("layer", i);

            int ssl = getSectorSuperlayer(sector, superlayer);
            map.computeIfAbsent(ssl, k -> new ArrayList<>()).add(new int[]{layer, wire, i + 1});
        }
        return map;
    }

    private static int getSectorSuperlayer(int sector, int superlayer) {
        return (sector - 1) * 6 + superlayer;
    }
    
    private static void getLayerEfficiencies(int run, Segment seg, Map<Integer, List<int[]>> map, DCGeant4Factory DcDetector) {
        if (seg == null) return;
        float[] result = new float[3];
        CalSegmentTrajectory trj = new CalSegmentTrajectory();
        trj.set_SegmentId(seg.get_Id());
        trj.set_Superlayer(seg.get_Superlayer());
        trj.set_Sector(seg.get_Sector());

        double[] trkDocas = new double[6];
        double[] trkBfield = new double[6];
        int[] matchHits   = new int[6];

        int[][] matchedHits   = new int[3][6];
        double[][] matchedDocas = new double[3][6];
        Arrays.stream(matchedHits).forEach(row -> Arrays.fill(row, -1));

        int ssl = getSectorSuperlayer(seg.get_Sector(), seg.get_Superlayer());

        for (int l = 0; l < 6; l++) {
            double z = DcDetector.getWireMidpoint(seg.get_Sector() - 1, seg.get_Superlayer() - 1, l, 0).z;
            double trkX = seg.get_fittedCluster().get_clusterLineFitSlope() * z
                        + seg.get_fittedCluster().get_clusterLineFitIntercept();
            
            
            if (trkX == 0) continue;
            WireTrajectoryResult wtr = getWireOnTrajectory(seg, l + 1, trkX, DcDetector);
            int[] trjWire = wtr.intValues;
            double[] trjDocas = wtr.doubleValues;

            for (int[] j : map.getOrDefault(ssl, Collections.emptyList())) {
                if (j[0] == l + 1) {
                    for (int wo = 0; wo < trjWire.length; wo++) {
                        if (Math.abs(trjWire[wo] - j[1]) == 0) {
                            matchedHits[wo][l] = j[2];
                            matchedDocas[wo][l] = trjDocas[wo];
                        }
                    }
                }
            }

            matchHits[l] = -1;
            for (int wo = 0; wo < trjWire.length; wo++) {
                if (matchedHits[wo][l] > 0) {
                    matchHits[l] = matchedHits[wo][l];
                    trkDocas[l] = matchedDocas[wo][l];
                }
            }
            if (matchHits[l] == -1) {
                int trjW=-1;
                if(trjDocas[0]<trjDocas[1]) {
                    trkDocas[l] = trjDocas[0]; 
                    trjW=trjWire[0];
                } else {
                    trkDocas[l] = trjDocas[1];
                    trjW=trjWire[1];
                }
                boolean statusOK = EventProcessor.getWireStatus(run, seg.get_Sector(), seg.get_Superlayer(), l+1, trjW);
                if(!statusOK) return;
            }
            double B = 0;
            if(seg.get_Region()==2) {
                swim.Bfield(seg.get_Sector(), trkX, 0, z, result);
                B = Math.sqrt(result[0]*result[0]+result[1]*result[1]+result[2]*result[2]);
            }
            trkBfield[l]=B;
        }

        trj.setTrkDoca(trkDocas);
        trj.setTrkBfield(trkBfield);
        trj.setMatchedHitId(matchHits);
        seg.set_Trajectory(trj);
    }

    
    
    public static List<SegmentTrajectoryRecord> getLayIneffBank(DataEvent event, int count, int runNumber) {
        DataBank bank = event.getBank("RUN::config");
        int run = bank.getInt("run", 0);
        List<FittedHit> hits = new ArrayList<>();
        List<FittedCluster> clusters = new ArrayList<>();
        List<Segment> segments = new ArrayList<>();
        // get segment property
        DataBank bnkHits = event.getBank("TimeBasedTrkg::TBHits");
        FittedHit hit = null;
        for (int i = 0; i < bnkHits.rows(); i++) {
            hit = getHit(bnkHits, i);
            if(hit!=null)
                hits.add(hit);
        }
        clusters = recomposeClusters(hits);
        segments =  getSegments(run, clusters, event, WireIneffAnalViewer.dcDetector);
        List<SegmentTrajectoryRecord> trajectoryRecords = new ArrayList<>();

        for (Segment aSeglist : segments) {
            if (aSeglist.get_Id() == -1) continue;
            CalSegmentTrajectory trj = (CalSegmentTrajectory) aSeglist.get_Trajectory();

            for (int l = 0; l < 6; l++) {
                trajectoryRecords.add(new SegmentTrajectoryRecord(
                    trj.get_SegmentId(),
                    trj.get_Sector(),
                    trj.get_Superlayer(),
                    (l+1),
                    trj.getMatchedHitId()[l],
                    trj.getTrkDoca()[l],
                    trj.getTrkBfield()[l]
                ));
            }
        }

        return trajectoryRecords;
    }

    private static double getNormalizedDoca(double trkX, double xn, double cellSize) {
        return Math.abs(trkX - xn) / cellSize;
    }

    
    public static class WireTrajectoryResult {
        public final int[] intValues;
        public final double[] doubleValues;

        public WireTrajectoryResult(int[] intValues, double[] doubleValues) {
            this.intValues = intValues;
            this.doubleValues = doubleValues;
        }
    }
}

