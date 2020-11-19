#! /bin/bash -f

#
# Data input file should just be a single column of numbers 
# Usage:
# bash calc_stats.bash <datafile>
#
# $ bash calc_stats.bash tmpf32.dat
# 82                   -- number of lines/entries
# 46020.020000         -- total sum
# 25940520.272400      -- sum of squares
# 561.219756,37.150956 -- average, std deviation
#
#
#
#
if [ $# -gt 1 ]
   then
      echo "greater than zero $#"
fi




if [ $# -gt 1 ]
    then
        echo "hey $#"
        argv1=$1
        argv2=$2
        n2=$(wc -l < $argv2)

        echo "$argv1, $argv2"
        if [ -e $argv2 ]
            then
               echo "okay, data is in second column of file $argv2"
               awk '{print $2}' $argv2 > $argv1
            else
               echo "$argv2 does not exist"
               exit
        fi
     else
     argv1=$1
     nint=$(wc -l < $argv1)
     
fi

echo "$nint"

avg=$(awk '{num += $0 } END { printf("%f",num) }' $argv1)
# set avg = `cat tmp1.dat`
echo "$avg"
sqrs=$(awk '{num += $0^2} END {printf("%f",num) }' $argv1)
# set sqrs = `cat tmp1.dat`
echo "$sqrs"
echo "$nint $avg $sqrs" > tmp1.dat
awk '{num = (($3/$1) - ($2/$1)^2)^0.5; num1 = $2/$1 } END {printf("%f,%f",num1,num) }' tmp1.dat

# remove temp file
rm -f tmp1.dat
