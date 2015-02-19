To execute Ising model follow below instructions (download all the codes from github)
goto https://github.com/tridip66/ising or you found this 'text' here
1) Download ising.f from src folder
2) compile it in HPCC development node by:
	f95 ising.f -o ising.exe
3) Execute it by ./ising.exe

4) Download makeplots.R
Load R module to plot different parameters with temperature
5) module load R
6) R < makeplots.R [--save --no-save --vanilla]

To visualize spin 
7) Download spinviz.m and plotspin.m
8) module load Matlab
9) matlab -nodesktop < spinviz.m
10) you can download 'report' folder 
11) Execute report by: pdflatex report.tex

