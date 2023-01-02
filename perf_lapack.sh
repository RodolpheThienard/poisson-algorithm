#!/bin/sh

# 8 -> 8^7
size=20
i=1
val=2

while [ $i -le $size ]
do
  val=`echo 2 \* $val | bc `
  ./bin/tpPoisson1D_direct $val >> perf.dat
  ((i++))
done


# all performances measured file
file="perf.dat"

# Removing older files
rm -f "dgbtrf&dgbtrs.dat"
rm -f "dgbtrftridiag&dgbtrs.dat"

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



# file with the addition of dgbtrftridiag and dgbtrs
# while IFS= read -r line
# do
#   if [[ "$line" = *"dgbtrftridiag :"* ]]
#   then 
#     num=`echo "$line" | awk -F ':' '{print $1}'`
#     dgbtrftridiag=`echo "$line" | awk -F ':' '{print $3}'`
#   fi
#   if [[ "$line" = *"dgbtrs"* ]]
#   then
#     dgbtrs=`echo "$line" | awk -F ':' '{print $3}'`
#     val=`echo "$dgbtrftritiag + $dgbtrs" | bc`
#     echo "$num : $val">> "dgbtrftridiag&dgbtrs.dat"
#   fi
  
# done<$file
