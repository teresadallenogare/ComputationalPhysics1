plot "dataX.txt" using 1:2 w l title "Body 1", "dataX.txt" using 3:4 w l title "Body 2", "dataX.txt" using 5:6 w l title "Body 3"
set xlabel "x"
set ylabel "y"
set title "Trajectories three bodies with alpha1 =  0.306893 and alpha2 = 0.125507"
set style line 2 lc rgb '#0060ad' pt 7 
set object circle at first 0,0 size scr 0.01 \
    fillcolor rgb 'black’ fillstyle solid
set object circle at first -1,0 size scr 0.01 \
    fillcolor rgb 'black’ fillstyle solid
set object circle at first 1,0 size scr 0.01 \
    fillcolor rgb 'black’ fillstyle solid
pause -1