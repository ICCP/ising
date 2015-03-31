rm data/*
rm images/*
rm spins.gif
make clean -f ising_makefile
make ising2.exe -f ising_makefile
./ising2.exe
python3 plot_spins.py
convert -delay 10 -loop 0 images/SpinData*.png spins.gif
