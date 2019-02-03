#!/usr/bin/gnuplot 

set term pdf monochrome size 15.0cm,9.0cm

set encoding iso_8859_1

set output "residues.pdf"

set grid

set xtics font "Times-Roman, 13"
set ytics font "Times-Roman, 13"

set xlabel "iterations" center
set ylabel "Norm(L_inf)" center
set title "Resíduos obtidos com cada malha"

#set offset graph 0.10,0.10,0.10,0.10
set xrange[0:10000]

set border lw 2

set pointsize 0.3

set key box font ",12"

plot "residue_5050.dat" u 1:2 with points lt -1 lw 0.1 pt 4 title "malha 50x50",\
     "residue_8080.dat" u 1:2 with lines lw 1.5 title "malha 80x80",\
     "residue_110110.dat" u 1:2 with linespoints lw 1.5 title "malha 110x110",\
     "residue_200200.dat" u 1:2 lw 1.9 title "malha 200x200"

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

set term pdf monochrome size 15.0cm,9.0cm

set output "pressures.pdf"

set encoding iso_8859_1

set grid

set xtics font "Times-Roman, 13"
set ytics font "Times-Roman, 13"

set xlabel "x" center
set ylabel "Pressão" center
set title "Comparação de pressões em y=0.5"

set offset graph 0.10,0.10,0.10,0.10

set xrange[0:3.5]
set border lw 2

set pointsize 0.3

set key left box font ",12"

plot "exact.dat" u 1:2 with linespoints lt -1 lw 0.1 pt 4 title "Exata",\
     "mid_pressure_8080.dat" u 1:2 with lines lw 1.5 title "malha 80x80",\
     "mid_pressure_8080.dat" u 1:2 with lines lw 1.6 title "malha 80x80",\
     "mid_pressure_110110.dat" u 1:2 with linespoints lw 1.7 title "malha 110x110",\
     "mid_pressure_200200.dat" u 1:2 lw 1.8 title "malha 200x200"


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

plot 'times.dat' w l

