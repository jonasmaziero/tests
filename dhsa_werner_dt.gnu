reset
set terminal postscript enhanced color "Times-Roman" 18
set output "dhsa_werner_dt.eps"

set xrange [2:100]
set yrange [0:930]
#set key right bottom
#set xlabel '1/h'
#set ylabel 't(s)'

plot 'dt_amp.dat' u 1:4 w lp ps 2 pt 13 notitle