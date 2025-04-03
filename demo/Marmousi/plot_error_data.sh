# calculate the true constrast source

dt=0.006
nr=227
nt=2401
f2=-100
d2=75
vr=3500
nto=1250

f2km=`bc -l << OUT
$f2*0.001
OUT`

d2km=`bc -l << OUT
$d2*0.001
OUT`

g=00001

name1=dataobs_org_$g.bin
name2=datae_${g}.bin
nameo=dataeerr_$g.bin

$MISC_BIN/diffmod << OUT
$name1 $name2 $nameo
$nt $nr 1
OUT

bclip=19.7775 
wclip=-19.626

rm dataobs_org_$g.binr
$MISC_BIN/reducetime1 << OUT
dataobs_org_$g.bin dataobs_org_$g.binr
$nt $nr $dt
$f2 $d2
0 $vr
$nto
OUT

rm dataobs_org_$g.binrg
$MISC_BIN/gain << OUT
dataobs_org_$g.binr dataobs_org_$g.binrg
$nto $nr $dt
1
0.5
0 1.
0
$f2 $d2
OUT

rm gather_dataobs_org_$g.ps
psimage < dataobs_org_$g.binrg > gather_dataobs_org_$g.ps n1=$nto f1=0 d1=$dt f2=$f2km d2=$d2km legend=0 d1num=2 f2num=-1000 d1num=2 d2num=5 \
label1='Time-off/3.5(s)' label2='Offset(km)' labelsize=10 titlesize=10 wbox=2 hbox=2 perc=98 \
titlefont='Times-Roman' labelfont='Times-Roman' bps=24 grid1=dot grid2=dot n1tic=2 n2tic=2

rm datar_$g.binr
$MISC_BIN/reducetime1 << OUT
datar_$g.bin datar_$g.binr
$nt $nr $dt
$f2 $d2
0 $vr
$nto
OUT

rm datar_$g.binrg
$MISC_BIN/gain << OUT
datar_$g.binr datar_$g.binrg
$nto $nr $dt
1
0.5
0 1.
0
$f2 $d2
OUT

rm gather_datar_$g.ps
psimage < datar_$g.binrg > gather_datar_$g.ps n1=$nto f1=0 d1=$dt f2=$f2km d2=$d2km legend=0 f1num=-100 d1num=2 f2num=-1000 d2num=5 \
label2='Offset(km)' labelsize=10 titlesize=10 wbox=2 hbox=2 bclip=$bclip wclip=$wclip \
titlefont='Times-Roman' labelfont='Times-Roman' bps=24 grid1=dot grid2=dot n1tic=2 n2tic=2

rm datae_$g.binr
$MISC_BIN/reducetime1 << OUT
datae_${g}.bin datae_$g.binr
$nt $nr $dt
$f2 $d2
0 $vr
$nto
OUT

rm datae_$g.binrg
$MISC_BIN/gain << OUT
datae_$g.binr datae_$g.binrg
$nto $nr $dt
1
0.5
0 1.
0
$f2 $d2
OUT

rm gather_datae_$g.ps
psimage < datae_$g.binrg > gather_datae_$g.ps n1=$nto f1=0 d1=$dt f2=$f2km d2=$d2km legend=0 f1num=-100 d1num=2  f2num=-1000 d2num=5 \
label2='Offset(km)' labelsize=10 titlesize=10  wbox=2 hbox=2 bclip=$bclip wclip=$wclip \
titlefont='Times-Roman' labelfont='Times-Roman' bps=24 grid1=dot grid2=dot n1tic=2 n2tic=2

rm dataeerr_$g.binr
$MISC_BIN/reducetime1 << OUT
dataeerr_$g.bin dataeerr_$g.binr
$nt $nr $dt
$f2 $d2
0 $vr
$nto
OUT

rm dataeerr_$g.binrg
$MISC_BIN/gain << OUT
dataeerr_$g.binr dataeerr_$g.binrg
$nto $nr $dt
1
0.5
0 1.
0
$f2 $d2
OUT

rm gather_dataeerr_$g.ps
psimage < dataeerr_$g.binrg > gather_dataeerr_$g.ps n1=$nto f1=0 d1=$dt f2=$f2km d2=$d2km legend=0 f1num=-100 d1num=2  f2num=-5 d2num=5 \
label2='Offset(km)' labelsize=10 titlesize=10  wbox=2 hbox=2 bclip=$bclip wclip=$wclip \
titlefont='Times-Roman' labelfont='Times-Roman' bps=24 grid1=dot grid2=dot n1tic=2 n2tic=2

sc=1
psmerge in=gather_dataeerr_$g.ps in=gather_datae_$g.ps in=gather_datar_$g.ps in=gather_dataobs_$g.ps \
translate=8,6 translate=6,6 translate=4,6 translate=2,6 translate=0,6 \
scale=$sc,$sc,$sc scale=$sc,$sc,$sc > fig_error_assimilation_$g.ps

convert fig_error_assimilation_$g.ps fig_error_assimilation_$g.pdf 
evince fig_error_assimilation_$g.ps

