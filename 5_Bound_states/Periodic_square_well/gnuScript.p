plot "dataE_K_numerov.txt" u 1:2 w l t "n = 0", "dataE_K_numerov.txt" u 1:3 w l t "n = 1", "dataE_K_numerov.txt" u 1:4 w l t "n = 2"
set xlabel "k"
set ylabel "E"

#plot "dataE_K.txt" u 1:2 w l, "dataE_K.txt" u 1:3 w l, "dataE_K.txt" u 1:4 w l
pause -1