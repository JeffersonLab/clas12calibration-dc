#!/bin/sh -f

OS=$(uname)

# set up environment
classPath="$PWD/target/*:$PWD/target/classes/*:lib/clas12detector-dc-1.0-SNAPSHOT.jar:lib/coat-libs-11.0.2-SNAPSHOT.jar"

# run 
java -Xmx4g -Xms1024m  -Dsun.java2d.pmoffscreen=false  -cp $classPath org.clas.detector.clas12calibration.viewer.Driver
