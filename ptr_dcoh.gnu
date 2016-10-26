reset
set terminal postscript enhanced color "Times-Roman" 18
set output "ptr_dcoh.eps"

set xrange [1.5:4.5]
set yrange [-0.36:-0.08]
set key right bottom
#set xlabel '1/h'
#set ylabel 't(s)'

x0 = NaN
y0 = NaN

plot 'Cnl3.dat' u (dx=$1-x0,x0=$1,$1-dx/2):(dy=$2-y0,y0=$2,dy/dx)  w p ps 1 t 'n=3',\
'Cnl4.dat' u (dx=$1-x0,x0=$1,$1-dx/2):(dy=$2-y0,y0=$2,dy/dx) w p ps 1 t 'n=4',\
'Cnl5.dat' u (dx=$1-x0,x0=$1,$1-dx/2):(dy=$2-y0,y0=$2,dy/dx) w p ps 1 t 'n=5',\
'Cnl6.dat' u (dx=$1-x0,x0=$1,$1-dx/2):(dy=$2-y0,y0=$2,dy/dx) w p ps 1 t 'n=6',\
'Cnl7.dat' u (dx=$1-x0,x0=$1,$1-dx/2):(dy=$2-y0,y0=$2,dy/dx) w p ps 1 t 'n=7',\
'Cnl8.dat' u (dx=$1-x0,x0=$1,$1-dx/2):(dy=$2-y0,y0=$2,dy/dx) w p ps 1 t 'n=8',\
'Cnl9.dat' u (dx=$1-x0,x0=$1,$1-dx/2):(dy=$2-y0,y0=$2,dy/dx) w p ps 1 t 'n=9',\
'Cnl10.dat' u (dx=$1-x0,x0=$1,$1-dx/2):(dy=$2-y0,y0=$2,dy/dx) w p ps 1 t 'n=10'