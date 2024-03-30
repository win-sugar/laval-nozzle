
rm -rf 0 postProcessing
rm -rf 0\.[0-9]* log.*

cp system/probes_org system/probes
cp system/sample_org system/sample

cp -r 0_org 0

blockMesh > log.blockMesh

setFields > log.setFields

rhoCentralFoam > /dev/null

