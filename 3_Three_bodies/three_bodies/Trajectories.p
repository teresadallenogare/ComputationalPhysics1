set xrange[-2:2]
set yrange [-2:2]
set pointsize 10
set style line 2 lc rgb '#0060ad' pt 7 
set object circle at first 0,0 size scr 0.01 \
    fillcolor rgb 'black’ fillstyle solid
set object circle at first -1,0 size scr 0.01 \
    fillcolor rgb 'black’ fillstyle solid
set object circle at first 1,0 size scr 0.01 \
    fillcolor rgb 'black’ fillstyle solid
 ntail = 5000 #number of points to draw in the tail
 ninc = 10  #number of increments between framee
do for [ii = 1: 10000 :ninc] {
 im  = ((ii - ntail) < 0 ? 1:ii-ntail)
 plot 'dataX.txt' using 1:2  every ::im::ii with lines title "Body 1", \
 'dataX.txt' using 3:4  every ::im::ii with lines title "Body 2", \
 'dataX.txt' using 5:6  every ::im::ii with lines title "Body 3"
}

pause -1