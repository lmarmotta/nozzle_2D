#!/usr/bin/gnuplot

set term pdf monochrome size 15.0cm,9.0cm

set encoding iso_8859_1

set datafile separator ","

set output "cp.pdf"

set origin 0.0,0.0

set grid

set xtics font "Times-Roman, 15"
set ytics font "Times-Roman, 15"

set xlabel "X" center font "Times-Roman, 16" 
set ylabel "P/Pt" center font "Times-Roman, 16" 

set yrange [0.2:1.20]
set xrange [-2.5:1.50]

set border lw 2

set pointsize 0.3

set key right box font "Times-Roman,14"

plot "pressure_exp.csv" using 1:2:(0.010) lt 1.5 lw 1.0 with circles title " Exp",\
     "pressure_com.dat" using 1:3 w l title " Numerical"
