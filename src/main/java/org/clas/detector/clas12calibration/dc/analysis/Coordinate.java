package org.clas.detector.clas12calibration.dc.analysis;


/**
 * 
 * @author ziegler
*/
import java.util.Arrays;

public class Coordinate {
    private final int[] size;  

    // Constructor
    public Coordinate(int... size) {
        this.size = size.clone(); 
    }

    @Override
    public int hashCode() {
        return Arrays.hashCode(size);  
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;  // Check for the same reference
        if (obj == null || getClass() != obj.getClass()) return false;  // Type check

        Coordinate other = (Coordinate) obj;
        return Arrays.equals(size, other.size);  
    }
}
