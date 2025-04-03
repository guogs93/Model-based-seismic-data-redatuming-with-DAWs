###########################################################################################################################
# Output is a highpass filtered Ricker wavelet (ricker_bl.bin) of nto sample and of sampling rate dt (s)
#Bandpass filter is controlled by fstoplo and fpasslo
###########################################################################################################################

nt=2048
dt=0.006
freqrick=4
nto=2401

###########################################################################################################################

ndt=`bc << OUT
$nto-$nt
OUT`

nt2=`bc << OUT
2*$nt
OUT`

nth=`bc << OUT
$nt2/2
OUT`

freqmax=`bc -l << OUT
0.5/$dt
OUT`

dfreq=`bc -l << OUT
$freqmax/($nt-1)
OUT`

nt1=`bc << OUT
$nt+1
OUT`

echo DTA = $dta

rm ricker.bin
$MISC_BIN/ricker << OUT
$nt $dt $freqrick
OUT

rm rickera.bin
$MISC_BIN/augmentmodel << OUT
ricker.bin rickera.bin
$nt 1
0 0 0 $nt
0
OUT

mv rickera.bin ricker.bin

rm ricker_f.bin
$PROBS_BIN/filmartin_v1 << OUT > out
ricker.bin ricker_f.bin
$nt2 1
$dt
1.7 5 0.5 1.5
4
1
OUT

xgraph < ricker_f.bin n1=$nt2 f1=0 d1=$dt pairs=2 grid1=dot grid2=dot

rm fspec
$PROBS_BIN/spectrum << OUT
ricker.bin fspec
0
$nt2
$dt
OUT

$MISC_BIN/dfreq << OUT
$nt2 $dt
OUT

rm fspec1
$PROBS_BIN/spectrum << OUT
ricker_f.bin fspec1
0
$nt2
$dt
OUT

#$MISC_BIN/dfreq << OUT
#$nt2 $dt
#OUT

rm fspec2
cat fspec fspec1 > fspec2

xgraph < fspec2 n1=$nth f1=0 d1=$dfreq pairs=2 style=normal nplot=2 linecolor=0,2  linewidth=4,2 grid1=dot grid2=dot

rm ricker_fw.bin
$MISC_BIN/window << OUT
ricker_f.bin ricker_fw.bin
$nt2 1 1
$nt1 $nt2 1
1 1 1
1 1 1
OUT

rm zero.bin
$MISC_BIN/genmod3d << OUT
zero.bin
$ndt 1 1
0.
0.
1 1 1
1 1 1
OUT

rm ricker_bl${freqrick}.bin
cat ricker_fw.bin zero.bin > ricker_bl${freqrick}.bin

xgraph < ricker_bl${freqrick}.bin n1=$nto f1=0 d1=$dt pairs=2  grid1=dot grid2=dot
