reset
set terminal postscript enhanced color "Times-Roman" 18
set output "ptr_time.eps"

#set xrange [1.5:5]
#set yrange [0:]
#set key right bottom
#set xlabel '1/h'
#set ylabel 't(s)'

#set logscale y

plot 'ptr_time.dat' u 1:4 w lp notitle