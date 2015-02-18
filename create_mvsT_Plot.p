
reset

set xlabel "Temperature, T"
set ylabel "Magnetization, m"

unset key

plot "avg_m_values" lt rgb "black" with linespoints,  "avg_m_values"  using 1:2:3 with yerrorbars
set terminal png size 400,300 enhanced font ",20"

set term png
set output "m_vs_T_plot.png"
replot


