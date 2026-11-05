#!/bin/bash
#SBATCH --job-name=alos2-pairtest
#SBATCH --partition=long
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=10:00:00
#SBATCH --exclude=xcng[0-1],xgpj0
#SBATCH --output=/home/t/tommy/pilot/logs/pairtest.log
set -uo pipefail
for i in $(seq 1 18); do [ -f ~/pilot/READY_FOR_ALOSSTACK ] && break; echo "waiting for archive ($i/18)"; sleep 300; done
[ -f ~/pilot/READY_FOR_ALOSSTACK ] || { echo "archive never became ready"; exit 1; }
REF=2018-10-12; SEC=2018-11-09
OUT=~/pilot/pairs/${REF}_${SEC}; mkdir -p $OUT
WORK=/temp/${SLURM_JOB_ID}; mkdir -p $WORK 2>/dev/null || { WORK=/tmp/pt_$$; mkdir -p $WORK; }
trap "cp $WORK/isce.log $OUT/ 2>/dev/null; cp $WORK/alos2Proc.xml $OUT/ 2>/dev/null; rm -rf $WORK" EXIT
echo "workdir $WORK on $(hostname)"
cp -r ~/pilot/ceos/$REF ~/pilot/ceos/$SEC $WORK/
cp ~/pilot/dem/elevation.dem.wgs84* ~/pilot/dem/elevation.dem.hdr $WORK/ 2>/dev/null
cd $WORK
cat > alos2App.xml <<XMLEOF
<?xml version="1.0" encoding="UTF-8"?>
<alos2App>
  <component name="alos2insar">
    <property name="reference directory">$WORK/$REF</property>
    <property name="secondary directory">$WORK/$SEC</property>
    <property name="reference polarization">HH</property>
    <property name="secondary polarization">HH</property>
    <property name="starting swath">2</property>
    <property name="ending swath">4</property>
    <property name="dem for coregistration">$WORK/elevation.dem.wgs84</property>
    <property name="dem for geocoding">$WORK/elevation.dem.wgs84</property>
  </component>
</alos2App>
XMLEOF
export PATH=$HOME/micromamba/envs/insar/bin:$HOME/micromamba/envs/insar/lib/python3.11/site-packages/isce/applications:$PATH
export OMP_NUM_THREADS=16
set +e
python3 $HOME/micromamba/envs/insar/lib/python3.11/site-packages/isce/applications/alos2App.py alos2App.xml
rc=$?
set -e
echo "alos2App exit: $rc"
find $WORK -maxdepth 2 \( -name "*.cor*" -o -name "filt_*" -o -name "*.unw*" -o -name "*.ion*" -o -name "*.los*" -o -name "*.xml" -o -name "*.vrt" \) -size -600M -exec cp --parents {} $OUT/ \; 2>/dev/null || true
du -sh $OUT
if [ $rc -eq 0 ]; then
  touch $OUT/DONE; echo PAIRTEST_OK
else
  echo PAIRTEST_FAILED
  exit $rc
fi
