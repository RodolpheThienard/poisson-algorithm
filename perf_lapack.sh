#!/bin/sh


# Removing older files
rm -f "dgbtrf&dgbtrs.dat"
rm -f "dgbtrftridiag&dgbtrs.dat"
rm -f "perf.dat"
rm -f "perf_iter.dat"

# 8 -> 8^7
size=10
i=1
h=10000
val=0

Direct method
while [ $i -le $size ]
do
  val=`echo $i \* $h | bc `
  ./bin/tpPoisson1D_direct $val >> perf.dat
  ((i++))
done
i=1
# Iter method
while [ $i -le $size ]
do
  val=`echo $i \* 0.1 | bc `
  ./bin/tpPoisson1D_iter $h $val >> perf_iter2.dat
  ((i++))
done


# all performances measured file
file="perf.dat"

# file with the addition of dgbtrf and dgbtrs
while IFS= read -r line
do
  if [[ "$line" = *"dgbtrf :"* ]]
  then 
    num=`echo "$line" | awk -F ':' '{print $1}'`
    dgbtrf=`echo "$line" | awk -F ':' '{print $3}'`
  fi
  if [[ "$line" = *"dgbtrs"* ]]
  then
    dgbtrs=`echo "$line" | awk -F ':' '{print $3}'`
    val=`echo "$dgbtrf + $dgbtrs" | bc`
    echo "$num : $val">> "dgbtrf&dgbtrs.dat"
  fi
 if [[ "$line" = *"dgbtrftridiag"* ]]
  then 
    dgbtrftridiag=`echo "$line" | awk -F ':' '{print $3}'`
    val=`echo "$dgbtrftridiag + $dgbtrs" | bc`
    echo "$num : $val">> "dgbtrftridiag&dgbtrs.dat"
  fi
  
done<$file

