plot "dataT.txt" using 1:2 w l title "|T|^2 with xi = 1", "dataT.txt" using 1:4 w l title "|T|^2 with xi = 9", \
"dataT.txt" using 1:6 w l title "|T|^2 with xi = 49", 
set xlabel "E^*"
set ylabel "|T|^2"
set title "TRANSMISSION COEFFICIENTS ASYMMETRIC BARRIER 2"
pause -1