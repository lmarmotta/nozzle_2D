#!/usr/bin/gnuplot 

set term pdf monochrome size 15.0cm,9.0cm

set output "comp_time.pdf"

set encoding iso_8859_1

set grid

set xtics font "Times-Roman, 13"
set ytics font "Times-Roman, 13"

set xlabel "tempo [s]" center
set ylabel "Quantidade de pontos na direção x" center
set title "Desempenho computacional do método"

set offset graph 0.10,0.10,0.10,0.10

plot 'times.dat' w lp lw 1.5

