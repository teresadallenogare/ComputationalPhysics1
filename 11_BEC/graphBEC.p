set title 'Wave functions'
set title font "Helvetica,16"
set xlabel 'r'
set xrange[0:6]
set xlabel font "Helvetica,14"
set ylabel 'wave function'
set ylabel font "Helvetica,14"
set grid 
plot 'wf_gs.txt' u 1:2 w l t 'xi = 1', 'wf_gs.txt' u 1:3 w l t 'xi = 10', 'wf_gs.txt' u 1:4 w l t 'xi = 100'

pause -1