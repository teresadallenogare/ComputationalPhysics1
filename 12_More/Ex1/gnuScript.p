plot "dataT.txt" using 1:2 w l title "|T|^2  for csi = 0.025", "dataT.txt" using 1:4 w l title " |T|^2 for csi = 0.005","dataT.txt" using 1:6 w l title "|T|^2 for csi = 0.0025"
set title "TRANSMISSION COEFFICIENT FOR DIFFERENT CSI"
set xlabel "E/|V_0|"
set ylabel "|T|^2"

pause -1