plot "dataTrans.txt" using 1:2 w l title "Rectangular potential", "dataTrans.txt" using 1:3 w l title "Gaussian potential",\
"dataTrans.txt" using 1:4 w l title "Asymm potential 1",  "dataTrans.txt" using 1:5 w l title "Asymm potential 2";
set xlabel "E/V0"
set title "Transmission coefficient"
set ylabel "T"
pause -1