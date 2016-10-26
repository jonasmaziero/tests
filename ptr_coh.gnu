reset
set terminal postscript enhanced color "Times-Roman" 18
set output "ptr_coh.eps"

set xrange [1.5:5]
set yrange [0:]
#set key right bottom
#set xlabel '1/h'
#set ylabel 't(s)'

plot 'Cnl3.dat' w lp ps 1 t 'n=3', 'Cnl4.dat' w lp ps 1 t 'n=4', 'Cnl5.dat' w lp ps 1 t 'n=5',\
'Cnl6.dat' w lp ps 1 t 'n=6', 'Cnl7.dat' w lp ps 1 t 'n=7',\
'Cnl8.dat' w lp ps 1 t 'n=8', 'Cnl9.dat' w lp ps 1 t 'n=9', 'Cnl10.dat' w lp ps 1 t 'n=10'