#!/usr/bin/gnuplot -persist

# Script gnuplot to iteratively plot residuals of nozzle_2d.

set termopt enhanced
set xtics font "Times-Roman, 10"
set ytics font "Times-Roman, 10"

set xlabel "Iterations [number]"
set ylabel "Log_{10} RHS_{Max}"
set title "Nozzle 2D Residue"

set style line 1 lt rgb "black" 
set style line 2 lt rgb "red" 

plot "residue.dat" u 1:2 with lines title "Continuity Residue"

pause 10; refresh; reread
