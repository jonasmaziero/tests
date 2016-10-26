reset
set terminal postscript enhanced color "Times-Roman" 18
set output "werner.eps"

set xrange [0:1]
set yrange [-0.01:0.7]
set key left top
#set xlabel 'w'

plot 'werner2.dat' u 1:2 w lp ps 1 t 'E_{n}', 'werner2.dat' u 1:3 w lp ps 1 t 'E_{hs}',\
'werner3.dat' u 1:2 w lp ps 1 t 'E_{n}', 'werner3.dat' u 1:3 w lp ps 1 t 'E_{hs}'