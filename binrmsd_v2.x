#! /bin/csh -f

# get an RMSD initial value
set firstrmsd = `awk 'NR>1{exit};{print $2}' rmsd.dat`

cat<<eof>b.awk
BEGIN {
val=$firstrmsd;
line=1;
}
{
if( \$1 < val )
    {I
#      print \$1 " " \$2
#        if( \$1 < val )
#        {
          val=\$1;
# store the line with the lowest energy score
          line=\$0;
#        }
     } 
}
END {
print line
}

eof
#echo "here"
egrep -v "nan" rmsd.dat | awk '{print $2 " " $1}' > tmp.dat
set loweststruct = `awk -f b.awk tmp.dat`
#echo "$loweststruct"
set lowpdb = `echo "$loweststruct" | awk '{print $2}'`
set reuval = `grep "$lowpdb" score.dat | awk '{print $1}'`
echo "$loweststruct REU $reuval"
#exit
#echo $lowpdb

# get rid of the old one
rm -f *.lowrmsd.pdb
cp ./struct/$lowpdb.pdb $lowpdb.lowrmsd.pdb

rm -f b.awk

if (-e qscore.dat) then
  awk '{printf "%5.4f\n", $2}' rmsd.dat | sort -n | sed "s/0\./ /g" | sed "s/\.//g" > tmp.dat
else
  awk '{printf "%4.3f\n", $2}' rmsd.dat | sort -n | sed "s/\.//g" | sed "s/ -000/ -/g" | sed "s/ -00/ -/g" | sed "s/ -0/ -/g" > tmp.dat
endif
egrep -v "scoring" tmp.dat | egrep -v "score" | sed "s/ 000/ /g" | sed "s/ 00/ 0/g" | sed "s/ 0/ /g" > sorted.dat

rm -f tmp.dat

# find highest and lowest values
set val_high = 0
set val_low = 0
#set step_size = 0

# number of bins to calculate for
set number_of_bins = 100
set setps = 0
set start = 1
set myval = 0

foreach val ("`cat sorted.dat`")
  @ myval=`expr $val`
#   echo "$myval"
  if ($start == 1) then
#    echo "hilo being $myval $val_high $val_low"
    @ val_low = $myval
    @ val_high = $myval
#   echo "hilo assigned  $val_high $val_low"
   @ start = 0
#   exit
  else
     if ($val_high < $myval) then
        @ val_high = $myval 
     else if ($val_low > $myval) then
        @ val_low = $myval
     endif
#     exit
  endif
end 

echo "hilo assigned $val_high $val_low"
@ steps = $val_high - $val_low
#echo "$steps steps"
set step_size = 0
@ step_size = `expr $steps / $number_of_bins`

#set an upper limit on the step size
if ($step_size > 1000) then # RSMD 1 Angstrom
   @ step_size = 1000
   # expand the number of bins
   @ number_of_bins = `expr $steps / 1000`
endif

# set a maximum energy of interest 200 Angstroms --unlikely
set max_REU = 200000 #miliREU's
#echo " $step_size "
set num_cnt = 0
set bin_cnt = 0
set bin_num = 0
set eval = 0
set tmpv = 0

@ bin_cnt += $step_size

foreach val ("`cat sorted.dat`")
   @ myval = `expr $val` 
   @ eval =  $val_low + $bin_cnt
   if ($myval <=  $eval) then
      @ num_cnt++;
#   echo "$eval > $myval"
  # exit
# set threshold
   else if ($myval < $max_REU) then
#   else
       # dump data
#      echo "$eval  $num_cnt"
      echo "$eval  $num_cnt" >> tmp.dat
       # find next step
        @ next_val = $eval 
      while ($myval > $next_val) 
        @ tmpv = $next_val + $step_size   
        @ next_val =  $tmpv
        @ tmpv = $bin_cnt + $step_size
        @ bin_cnt = $tmpv
        # populate empty bins
        if ($myval > $next_val) then
#          echo "$next_val  0"
          echo "$next_val  0" >> tmp.dat
        endif
#        echo "inside"
#        exit
      end
#      echo "bin cnt $bin_cnt for $myval"
      @ num_cnt = 1 # reset count
      @ eval = $val_low + $bin_cnt
#     exit
   endif
#exit
end

#echo "$eval  $num_cnt"
echo "$eval  $num_cnt" >> tmp.dat
rm -f sorted.dat

#awk '{print $1/1000 " " $2}' tmp.dat > bin_rmsd.dat
#rm -f tmp.dat

if (-e qscore.dat) then
   awk '{print $1/10000 " " $2}' tmp.dat > bin_rmsd.dat
   set xaxis_label = "TM-score %"
else
   awk '{print $1/1000 " " $2}' tmp.dat > bin_rmsd.dat
   set xaxis_label = "RMSD (Angstroms)"
endif

rm -f tmp.dat

cat<<eof>bfile_bin_rmsd.xmg
READ NXY "bin_rmsd.dat"
title "Frequency of Structures by RMSD to Lowest"
subtitle "$number_of_bins Bins, $loweststruct, REU $reuval"
s0 legend "structs"
xaxis label "$xaxis_label"
yaxis label "Frequency"
VIEW XMIN 0.25
VIEW XMAX 0.65
VIEW YMIN 0.45
VIEW YMAX 0.85
LEGEND 0.15, 0.35
AUTOSCALE
PRINT TO "bin_rmsd.png"
HARDCOPY DEVICE "PNG"
PAGE SIZE 512, 512
DEVICE "PNG" DPI 300
DEVICE "PNG" OP "transparent:on"
DEVICE "PNG" OP "compression:9"
eof

awk 'NR>100{exit};{print $0}' bin_rmsd.dat > bin_rmsd_100.dat

sed "s/bin_rmsd/bin_rmsd_100/g" bfile_bin_rmsd.xmg > bfile_bin_rmsd_100.xmg

chmod 777 *
