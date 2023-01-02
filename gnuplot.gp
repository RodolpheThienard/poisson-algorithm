set terminal png size 700,500
set xtics rotate by 45 right
#set style fill solid 1 border 0
#set style data histograms
set datafile separator ':'
#set yrange [0:*]
set ylabel "Time in nanosec"
set xlabel "Size N"
set key left top
#set boxwidth 0.9
set logscale y 2

#set output 'dgbtrf_dgbtrs_dgbsv.png' 
#set multiplot layout 3,1 rowsfirst
#set xtics norotate center 
#plot "perf.dat" every 3::0 using 3:xtic(1) ti "dgbtrf" with lp
#plot "perf.dat" every 3::1 using 3:xtic(1) ti "dgbtrs" with lp
#plot "perf.dat" every 3::2 using 3:xtic(1) ti "dgbsv" with lp
#unset multiplot

set output 'dgbtrf_dgbtrftridiag.png'
set title "comparaison du temps d'exécution entre dgbtrf et dgbtrftridiag en log2"
plot "perf.dat" every 4::0 using 3:xtic(1) ti "dgbtrf" with lp,\
"perf.dat" every 4::3 using 3:xtic(1) ti "dgbtrftridiag" with lp 

set output 'solution_perf.png'
set title "comparaison du temps d'exécution entre dgbsv dgbtrs + dgbtrftridiag en log2"
plot "perf.dat" every 4::2 using 3:xtic(1) ti "dgbsv" with lp,\
"dgbtrftridiag&dgbtrs.dat" using 2:xtic(1) ti "dgbtrftridiag + dgbtrs" with lp


set output 'dgbsv_dgbtrf.png'
set title "comparaison du temps d'exécution entre dgbsv, dgbtrs + dgbtrf en log2"
plot "perf.dat" every 4::2 using 3:xtic(1) ti "dgbsv" with lp,\
"dgbtrf&dgbtrs.dat" using 2:xtic(1) ti "dgbtrf + dgbtrf" with lp,\
