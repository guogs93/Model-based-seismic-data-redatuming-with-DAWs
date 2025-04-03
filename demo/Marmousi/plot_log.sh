
dt=0.002
nr=851
nt=4801
f2=-100
d2=20

g=00001

ibs=`bc << OUT
$nt*4
OUT`

for skip in 50 150 250 350 450 550
do

rm obs_${skip}.bin
dd if=dataobs_new_${g}.bin of=obs_${skip}.bin skip=$skip count=1 ibs=$ibs

rm red_${skip}.bin
dd if=dataobs_redatum_${g}.bin of=red_${skip}.bin skip=$skip count=1 ibs=$ibs

rm error_${skip}.bin
$MISC_BIN/diffmod << OUT
red_${skip}.bin obs_${skip}.bin error_${skip}.bin
$nt 1 1
OUT

rm all_${skip}.bin
cat obs_${skip}.bin red_${skip}.bin error_${skip}.bin > all_${skip}.bin

$SU_BIN/psgraph < all_${skip}.bin > all_${skip}.ps n1=$nt f1=0 d1=$dt pairs=2 nplot=3 hbox=3.8 wbox=1.5 f1num=0 d1num=2 f2num=-1.0 d2num=0.5 label1='Time (s)' label2='Amplitude' \
labelfont='Times-Roman' titlefont='Times-Roman' style=seismic linecolor=black,blue,red titlesize=18 labelsize=18

done

rm fig_datafit_${g}.ps
sc=1
$SU_BIN/psmerge in=all_50.ps in=all_150.ps in=all_250.ps in=all_350.ps in=all_450.ps in=all_550.ps \
translate=0,0 translate=1.5,0 translate=3,0 \
translate=4.5,0 translate=6,0 translate=7.5,0 \
scale=$sc,$sc scale=$sc,$sc > fig_datafit_${g}.ps

convert fig_datafit_${g}.ps fig_datafit_${g}.pdf
evince fig_datafit_${g}.pdf
