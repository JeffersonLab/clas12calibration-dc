#!/bin/sh -f

OS=$(uname)

# set up environment
classPath="$PWD/target/*:$PWD/target/classes/*"

# run 
#java -Xmx32g -Xms24g  -Dsun.java2d.pmoffscreen=false -Djava.util.logging.config.file=$CLAS12DIR/etc/logging/debug.properties -cp $classPath org.clas.detector.clas12calibration.viewer.Driver
java -Xmx8g -Xms4g  -Dsun.java2d.pmoffscreen=false -Djava.util.logging.config.file=$CLAS12DIR/etc/logging/debug.properties -cp $classPath org.clas.detector.clas12calibration.viewer.Driver

