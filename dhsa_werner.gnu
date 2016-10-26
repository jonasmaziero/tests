reset
set terminal postscript enhanced color "Times-Roman" 18
set output "dhsa_werner.eps"

set xrange [-1:1]
set yrange [0:1]
#set key right bottom
#set xlabel '1/h'
#set ylabel 't(s)'

plot 'dhsa_werner_d4.dat' u 1:2 w l notitle, 'dhsa_werner_d4.dat' u 1:3 w p ps 1 notitle,\
'dhsa_werner_d9.dat' u 1:2 w l notitle, 'dhsa_werner_d9.dat' u 1:3 w p ps 1 notitle,\
'dhsa_werner_d16.dat' u 1:2 w l notitle, 'dhsa_werner_d16.dat' u 1:3 w p ps 1 notitle,\
'dhsa_werner_d25.dat' u 1:2 w l notitle, 'dhsa_werner_d25.dat' u 1:3 w p ps 1 notitle,\
'dhsa_werner_d36.dat' u 1:2 w l notitle, 'dhsa_werner_d36.dat' u 1:3 w p ps 1 notitle,\
'dhsa_werner_d49.dat' u 1:2 w l notitle, 'dhsa_werner_d49.dat' u 1:3 w p ps 1 notitle,\
'dhsa_werner_d64.dat' u 1:2 w l notitle, 'dhsa_werner_d64.dat' u 1:3 w p ps 1 notitle,\
'dhsa_werner_d81.dat' u 1:2 w l notitle, 'dhsa_werner_d81.dat' u 1:3 w p ps 1 notitle,\
'dhsa_werner_d100.dat' u 1:2 w l notitle, 'dhsa_werner_d100.dat' u 1:3 w p ps 1 notitle