
dt=0.006
nr=341
nt=2401
f2=-100
d2=50
vr=3500
nto=1250

f2km=`bc -l << OUT
$f2*0.001
OUT`

d2km=`bc -l << OUT
$d2*0.001
OUT`

g=00001

bclip=39.8034 
wclip=-37.2624

rm dataobs_new_$g.binr
$MISC_BIN/reducetime1 << OUT
dataobs_new_$g.bin dataobs_new_$g.binr
$nt $nr $dt
$f2 $d2
0 $vr
$nto
OUT

rm dataobs_new_$g.binrg
$MISC_BIN/gain << OUT
dataobs_new_$g.binr dataobs_new_$g.binrg
$nto $nr $dt
1
0.5
0 1.
0
$f2 $d2
OUT

rm gather_dataobs_new_$g.ps
psimage < dataobs_new_$g.binrg > gather_dataobs_new_$g.ps n1=$nto f1=0 d1=$dt f2=$f2km d2=$d2km legend=0 d1num=2 f2num=-1000 d1num=2 d2num=5 \
label1='Time-off/3.5(s)' label2='Offset(km)' labelsize=10 titlesize=10 wbox=2 hbox=2 perc=98 \
titlefont='Times-Roman' labelfont='Times-Roman' bps=24 grid1=dot grid2=dot n1tic=2 n2tic=2

rm dataobs_redatum_$g.binr
$MISC_BIN/reducetime1 << OUT
dataobs_redatum_$g.bin dataobs_redatum_$g.binr
$nt $nr $dt
$f2 $d2
0 $vr
$nto
OUT

rm dataobs_redatum_$g.binrg
$MISC_BIN/gain << OUT
dataobs_redatum_$g.binr dataobs_redatum_$g.binrg
$nto $nr $dt
1
0.5
0 1.
0
$f2 $d2
OUT

rm gather_dataobs_redatum_$g.ps
psimage < dataobs_redatum_$g.binrg > gather_dataobs_redatum_$g.ps n1=$nto f1=0 d1=$dt f2=$f2km d2=$d2km legend=0 f1num=-100 d1num=2 f2num=-1000 d2num=5 \
label2='Offset(km)' labelsize=10 titlesize=10 wbox=2 hbox=2 bclip=$bclip wclip=$wclip \
titlefont='Times-Roman' labelfont='Times-Roman' bps=24 grid1=dot grid2=dot n1tic=2 n2tic=2

rm diff.binrg
$MISC_BIN/diffmod << OUT
dataobs_new_$g.binrg dataobs_redatum_$g.binrg diff.binrg
$nto $nr 1
OUT

rm gather_diff_$g.ps
psimage < diff.binrg > gather_diff_$g.ps n1=$nto f1=0 d1=$dt f2=$f2km d2=$d2km legend=0 f1num=-100 d1num=2 f2num=-1000 d2num=5 \
label2='Offset(km)' labelsize=10 titlesize=10 wbox=2 hbox=2 bclip=$bclip wclip=$wclip \
titlefont='Times-Roman' labelfont='Times-Roman' bps=24 grid1=dot grid2=dot n1tic=2 n2tic=2

rm fig_error_redatum_$g.ps
sc=1
psmerge in=gather_diff_$g.ps in=gather_dataobs_redatum_$g.ps in=gather_dataobs_new_$g.ps \
translate=4,6 translate=2,6 translate=0,6 \
scale=$sc,$sc,$sc scale=$sc,$sc,$sc > fig_error_redatum_$g.ps

convert fig_error_redatum_$g.ps fig_error_redatum_$g.pdf 
evince fig_error_redatum_$g.ps

