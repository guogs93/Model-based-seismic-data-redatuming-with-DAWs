n1=71
n2=341
h=50

ns=1
xs0=100
ds=8000
zs=450

nr=341
xr0=0
dr=50
zr=50

rm ftopo.ascii
rm ftopo.bin

$MISC_BIN/detectbathy << OUT > detect.out
vtrue.bin
$n1 $n2 1 $h
1500
$ns 1 $zs $xs0 0. $ds 0.
$nr 1 $zr $xr0 0. $dr 0.
0
OUT

rm acqui_org_datum.txt
$MISC_BIN/buildacqui2d_ffwi2d << OUT
$n1 $n2 $h
$ns $zs $xs0 0 $ds
$nr $zr $xr0 0 $dr
OUT

zr=450

mv acqui.txt acqui_org_datum.txt

rm acqui_new_datum.txt
$MISC_BIN/buildacqui2d_ffwi2d << OUT
$n1 $n2 $h
$ns $zs $xs0 0 $ds
$nr $zr $xr0 0 $dr
OUT

mv acqui.txt acqui_new_datum.txt


rm vgrad1.bin
$MISC_BIN/genmodgrad1 << OUT
$n1 $n2 $h 1500 0 1131 0.8196 450
OUT

mv vgrad1.bin vinit.bin

echo 
echo ------------------------------------------------------
echo BUILD ACQUISITION FILE
echo ------------------------------------------------------

echo
echo ===================================================================
echo ===================================================================
echo XL XC NS NR OFFSET SPREADL = $xl $xc $ns $nr $offset $spreadl
echo ===================================================================
echo ===================================================================
echo
