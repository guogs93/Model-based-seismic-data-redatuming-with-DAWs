n1=176
n2=851
h=0.02
nt=4801
dt=0.002

g=00001

bclip=0.04 
wclip=-0.04

for skip in 200 400 600 800
do
ibs=`bc << OUT
$n1*$n2*4
OUT`

rm snapobs_${skip}.bin
dd if=wavefield_obs${g}.bin of=snapobs_${skip}.bin skip=$skip count=1 ibs=$ibs

rm snapobs_${g}_${skip}.ps
psimage < snapobs_${skip}.bin > snapobs_${g}_${skip}.ps n1=$n1 f1=0 d1=$h f2=0 d2=$h legend=0 d1num=1 d2num=2 f1num=0 f2num=0 lwidth=0.1 lheight=1  lstyle=vertright \
labelsize=10 titlesize=10 wbox=4.85 hbox=1 bclip=$bclip wclip=$wclip  grid1=dot grid2=dot \
bps=24 label1='Depth (km)' label2='Distance (km)' titlefont='Times-Roman' labelfont='Times-Roman' \
labelsize=14 titlesize=14 curve=ftopo.dat npair=2 curvecolor=yellow curvewidth=1

rm snapr_${skip}.bin
dd if=wavefield_ini${g}.bin of=snapr_${skip}.bin skip=$skip count=1 ibs=$ibs

rm snapr_${g}_${skip}.ps
psimage < snapr_${skip}.bin > snapr_${g}_${skip}.ps n1=$n1 f1=0 d1=$h f2=0 d2=$h legend=0 d1num=1 d2num=2 f1num=0 f2num=0 lwidth=0.1 lheight=1  lstyle=vertright \
labelsize=10 titlesize=10 wbox=4.85 hbox=1 bclip=$bclip wclip=$wclip grid1=dot grid2=dot \
bps=24  label1='Depth (km)' label2='Distance (km)' titlefont='Times-Roman' labelfont='Times-Roman' \
labelsize=14 titlesize=14 curve=ftopo.dat npair=2 curvecolor=yellow  curvewidth=1

rm snape_${skip}.bin
dd if=wavefield_exd${g}.bin of=snape_${skip}.bin skip=$skip count=1 ibs=$ibs

rm snape_${g}_${skip}.ps
psimage < snape_${skip}.bin > snape_${g}_${skip}.ps n1=$n1 f1=0 d1=$h f2=0 d2=$h legend=0 d1num=1 d2num=2 f1num=0 f2num=0 lwidth=0.1 lheight=1  lstyle=vertright \
labelsize=10 titlesize=10 wbox=4.85 hbox=1 bclip=$bclip wclip=$wclip  grid1=dot grid2=dot \
bps=24  label1='Depth (km)' label2='Distance (km)' titlefont='Times-Roman' labelfont='Times-Roman' \
labelsize=14 titlesize=14 curve=ftopo.dat npair=2 curvecolor=yellow  curvewidth=1

done

rm fig_wavefield_${g}.ps
sc=1
psmerge in=snapobs_${g}_200.ps in=snapr_${g}_200.ps in=snape_${g}_200.ps \
in=snapobs_${g}_400.ps in=snapr_${g}_400.ps in=snape_${g}_400.ps \
in=snapobs_${g}_600.ps in=snapr_${g}_600.ps in=snape_${g}_600.ps \
in=snapobs_${g}_800.ps in=snapr_${g}_800.ps in=snape_${g}_800.ps \
translate=0,18 translate=0,17 translate=0,16 \
translate=4.85,18 translate=4.85,17 translate=4.85,16 \
translate=0,15 translate=0,14 translate=0,13 \
translate=4.85,15 translate=4.85,14 translate=4.85,13 \
scale=$sc,$sc scale=$sc,$sc > fig_wavefield_${g}.ps

convert fig_wavefield_${g}.ps fig_wavefield_${g}.pdf
evince fig_wavefield_${g}.ps

