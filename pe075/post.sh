
sed '/enabled/ s/false/true/' system/sample_org > system/sample

sed '/enabled/ s/true/false/' system/probes_org > system/probes

rhoCentralFoam  -postProcess -latestTime

