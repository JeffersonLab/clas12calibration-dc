#!/bin/sh -f

OS=$(uname)

# set up environment
classPath="$PWD/lib/*"

# install local jar
mvn install:install-file -Dfile=lib/clas12detector-dc-1.0-SNAPSHOT.jar -DgroupId=org.jlab.clas12.detector -DartifactId=clas12detector-dc -Dversion=1.0-SNAPSHOT -Dpackaging=jar
mvn install:install-file -Dfile=lib/coat-libs-11.0.2-SNAPSHOT.jar -DgroupId=org.jlab.coat -DartifactId=coat-libs -Dversion=11.0.2-SNAPSHOT -Dpackaging=jar
mvn install
