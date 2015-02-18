reset


set xlabel "Millions of Steps"
set ylabel "Magnetization, m"

set key right top
set key outside

plot "Results1.0" using 1:2 with line title "T=1.0"
replot "Results2.0" using 1:2 with line title "T=2.0"
replot "Results2.1" using 1:2 with line title "T=2.1"
replot "Results2.2" using 1:2 with line title "T=2.2"
replot "Results2.3" using 1:2 with line title "T=2.3"
replot "Results2.4" using 1:2 with line title "T=2.4"
replot "Results2.5" using 1:2 with line title "T=2.5"
replot "Results2.6" using 1:2 with line title "T=2.6"
replot "Results2.7" using 1:2 with line title "T=2.7"
replot "Results2.8" using 1:2 with line title "T=2.8"
replot "Results2.9" using 1:2 with line title "T=2.9"
replot "Results3.0" using 1:2 with line title "T=3.0"
replot "Results4.0" using 1:2 with line title "T=4.0"
set terminal png size 400,300 enhanced font ",20"
set output "m_vs_Steps_plot.png"
replot
