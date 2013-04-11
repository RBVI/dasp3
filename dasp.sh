#! /bin/bash
export JAVA_HOME=/usr/java/default
# java -Xmx16G -Xss25m -jar /usr/local/lib/dasp.jar $@
java -Xmx16G -Xss25m -jar jar/dasp.jar $@
