plot "dataU.txt" using 1:2 w l title "potential energy system", "dataU.txt" using 1:3 w l title "kinetic energy system",\
"dataU.txt" using 1:4  w l title "total energy system"
set arrow 1 from 0, 81.6 to 10, 81.6 nohead ls 2 linewidth 1 linecolor rgb "#bd2129"
set arrow 2 from 0, 234 to 10, 234 nohead ls 2 linewidth 1 linecolor rgb "#0000ff"
set xlabel "t"
set ylabel "Energies(t)"
set title "ENERGY OF THE SYSTEM DURING EQUILIBRATION (rho = 0.7, T = 0.5)"
pause -1