#!/bin/tcsh
#
# Example: ./cal.sh orbit number
#
echo $1
#
if (! -e /pub/c/oura/HRSC-sim/s_$1_sim.bin) then
cat ./SHARAD-bin/s_$1_rgram.bin | ./bin-bscan > obs.txt
cat ./SHARAD-bin/s_$1_rgram.bin | ./bin-orb   > obs-orb.txt
rm -rf data
mkdir data
cat ./SHARAD-bin/s_$1_rgram.bin | ./bin-gm_HRSC 0 10000
foreach i (./data/????.gm)
echo $i:r
cat $i:r.gm | ./gm-de > $i:r.de
end
cat ./SHARAD-bin/s_$1_rgram.bin | ./de-bin > /pub/c/oura/HRSC-sim/s_$1_sim.bin
cat /pub/c/oura/HRSC-sim/s_$1_sim.bin | ./bin-bscan > sim.txt
#
endif
#
