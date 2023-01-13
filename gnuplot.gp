set terminal svg size 700,500
set xtics rotate by 45 right
#set style fill solid 1 border 0
#set style data histograms
set datafile separator ':'
#set yrange [0:*]
set ylabel "Erreur"
set xlabel "Taille N*N"
set key left top
#set boxwidth 0.9

set output 'forward_error.svg'
set title "evoulution de l'erreur avant en fonction de la taille de la matrice"
plot "dgbtrf&dgbtrs_err.dat" using 2:xtic(1) ti "dgbtrf + dgbtrs" with lp,\
"dgbtrftridiag&dgbtrs_err.dat" using 2 ti "dgbtrftridiag + dgbtrs" with lp,\
"dgbsv_err.dat" using 2 ti "dgbsv" with lp

set ylabel "Temps en seconde"

set output 'dgbtrf_dgbtrftridiag.svg'
set title "comparaison du temps d'exécution entre dgbtrf et dgbtrftridiag"
plot "perf.dat" every 5::0 using 3:xtic(1) ti "dgbtrf" with lp,\
"perf.dat" every 5::3 using 3:xtic(1) ti "dgbtrftridiag" with lp

set output 'solution_perf.svg'
set title "comparaison du temps d'exécution entre dgbsv dgbtrs + dgbtrftridiag"
plot "perf.dat" every 5::2 using 3:xtic(1) ti "dgbsv" with lp,\
"dgbtrftridiag&dgbtrs.dat" using 2:xtic(1) ti "dgbtrftridiag + dgbtrs" with lp,\
"dgbtrf&dgbtrs.dat" using 2:xtic(1) ti "dgbtrf + dgbtrs" with lp


set output 'dgbsv_dgbtrf.svg'
set title "comparaison du temps d'exécution entre dgbsv, dgbtrs + dgbtrf"
plot "perf.dat" every 5::2 using 3:xtic(1) ti "dgbsv" with lp,\
"dgbtrf&dgbtrs.dat" using 2:xtic(1) ti "dgbtrf + dgbtrs" with lp

set key right top

set logscale x 2
set logscale y 2
set output "convergence.svg"
set title "Analyse de la convergence en Itération"
set xlabel "Itération"
set ylabel "Résidu"
plot "RESVEC_alpha.dat" using 1 ti "Richardson Alpha" with lp,\
"RESVEC_Jacobi.dat" using 1 ti "Richardson Jacobi" with lp,\
"RESVEC_Gauss.dat" using 1 ti "Richardson Gauss-Seidel" with lp

unset logscale
set output "convergence_temps.svg"
set title "Analyse de la convergence en temps"
set xlabel "Taille N*N"
set ylabel "Temps en seconde"
plot "perf_iter.dat" every 3::0 using 3:xtic(1) ti "Richardson Alpha" with lp,\
"perf_iter.dat" every 3::1 using 3 ti "Richardson Jacobi" with lp,\
"perf_iter.dat" every 3::2 using 3 ti "Richardson Gauss-Seidel" with lp

set datafile separator " "
set output "err.svg"
set title "Delta entre Analytique sol et sol MB"
set ylabel "Delta Erreur"
set xlabel "Taille N*N"
plot "SOL_JACOBI.dat" using ($2-$1) ti "Richardson Jacobi" with lp,\
"SOL_GAUSS.dat" using ($2-$1) ti "Richardson Gauss" with lp,\
"SOL_ALPHA.dat" using ($2-$1) ti "Richardson Alpha" with lp

set datafile separator ":"
set output "temps_scilab.svg"
set title "Temps d'éxécution Richardson alpha en fonction de N"
set ylabel "Temps en seconde"
set xlabel "Taille N*N"
plot "scilab.dat" using 2:xtic(1) ti "Scilab" with lp,\
"dgbtrftridiag&dgbtrs.dat" using 2:xtic(1) ti "C" with lp
