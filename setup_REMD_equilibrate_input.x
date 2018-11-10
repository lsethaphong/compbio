#!/bin/bash -f

if[ -f groupfile ];
then
  rm groupfile
fi

nrep='wc temperatures.dat | awk '{print $1}'
echo $nrep
count=0
for TEMP in 'cat temperatures.dat'
do
  let COUNT+=1
  REP='printf "%03d" $COUNT'
  echo "TEMPERATURE: $TEMP K ==> FILE: prod.mdin.$REP"
  sed "s/XXXXX/$TEMP/g" REMD_RNA.mdin > temp
  sed "s/RANDOM_NUMBER/$RANDOM/g" temp > prod.mdin.$REP
  echo "-O -rem 0 -i prod.mdin.$REP -o prod.mdout.$REP -c heat2.rst.$REP -ref heat2.rst.$REP -r prod.rst.$REP -x prod.mdcrd.$REP -inf mdinfo.$REP -p generic.prmtop" >> production.groupfile

  rm -f temp

  echo "-O -rem 0 -i min1.in -o min1.out.$REP -c generic.inpcrd -ref generic.inpcrd -r min1.rst.$REP -x min1.mdcrd.$REP -inf min1.mdinfo.$REP -p generic.prmtop" >> min1.groupfile

  sed "s/XXXXX/$TEMP/g" REMD_RNA.heat > temp
  sed "s/RANDOM_NUMBER/$RANDOM/g" temp > heat.in.$REP
  echo "-O -rem 0 -i heat.in.$REP -o heat.out.$REP -c min1.rst.$REP -ref min1.rst.$REP -r heat.rst.$REP -x heat.mdcrd.$REP -inf heat.mdinfo.$REP -p generic.prmtop" >> heat.groupfile

  rm -f temp

  sed "s/XXXXX/$TEMP/g" REMD_RNA.mdpme > temp
  sed "s/RANDOM_NUMBER/$RANDOM/g" temp > mdpme.in.$REP
  echo "-O -rem 0 -i mdpme.in.$REP -o mdpme.mdout.$REP -c heat.rst.$REP -ref heat.rst.$REP -r mdpe.rst.$REP -x mdpme.mdcrd.$REP -inf mdpme.mdinfo.$REP -p generic.prmtop" >> mdpme.groupfile

  rm -f temp

  echo "-O -rem 0 -i min2.in -o min2.out.$REP -c mdpme.rst.$REP -ref mdpme.rst.$REP -r min2.rst.$REP -x min2.mdcrd.$REP -inf min2.mdinfo.$REP -p generic.prmtop" >> min2.groupfile
  
  sed "s/XXXXX/$TEMP/g" REMD_RNA.mdpme2 > temp
  sed "s/RANDOM_NUMBER/$RANDOM/g" temp > mdpme2.in.$REP
  echo "-O -rem 0 -i mdpme2.in.$REP -o mdpme2.mdout.$REP -c min2.rst.$REP -ref min2.rst.$REP -r mdpme2.rst.$REP -x mdpme2.mdcrd.$REP -inf mdpme2.mdinfo.$REP -p generic.prmtop" >> mdpme2.groupfile

  rm -f temp

  echo "-O -rem 0 -i min3.in -o min3.out.$REP -c mdpme2.rst.$REP -ref mdpme2.rst.$REP -r min3.rst.$REP -x min3.mdcrd.$REP -inf min3.mdinfo.$REP -p generic.prmtop" >> min3.groupfile

  echo "-O -rem 0 -i min4.in -o min4.out.$REP -c min3.rst.$REP -ref min3.rst.$REP -r min4.rst.$REP -x min4.mdcrd.$REP -inf min4.mdinfo.$REP -p generic.prmtop" >> min4.groupfile

  echo "-O -rem 0 -i min5.in -o min5.out.$REP -c min4.rst.$REP -ref min4.rst.$REP -r min5.rst.$REP -x min5.mdcrd.$REP -inf min5.mdinfo.$REP -p generic.prmtop" >> min5.groupfile

  echo "-O -rem 0 -i min6.in -o min6.out.$REP -c min5.rst.$REP -ref min5.rst.$REP -r min6.rst.$REP -x min6.mdcrd.$REP -inf min6.mdinfo.$REP -p generic.prmtop" >> min6.groupfile

  echo "-O -rem 0 -i min7.in -o min7.out.$REP -c min6.rst.$REP -ref min6.rst.$REP -r min7.rst.$REP -x min7.mdcrd.$REP -inf min7.mdinfo.$REP -p generic.prmtop" >> min7.groupfile



  sed "s/XXXXX/$TEMP/g" REMD_RNA.heat2 > temp
  sed "s/RANDOM_NUMBER/$RANDOM/g" temp > heat2in.$REP
  echo "-O -rem 0 -i heat2.in.$REP -o heat2.out.$REP -c min7.rst.$REP -ref min7.rst.$REP -r heat2.rst.$REP -x heat2.mdcrd.$REP -inf heat2.mdinfo.$REP -p generic.prmtop" >> heat2.groupfile

  rm -f temp

done

echo "#" >> groupfile

echo "N REPLICANTS = $nrep"
echo " Done."
