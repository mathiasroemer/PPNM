libs = ../matlib

CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(CODE)
TIME = time --output=$@ --append --format "%e %U"

default : Out.txt Out.wave.svg Out.convergence.svg

Out.txt Out.wave.data: main.exe
	mono $< 1>Out.txt 2>Out.wave.data

Out.conv.rmin.data : mainConv.exe
	for rmin in $$(seq 0.025 0.025 0.9); do\
		mono $< >> $@ -rmin:$$rmin -rmax:8 -acc:0.01 -eps:0.01; \
	done

Out.conv.rmax.data : mainConv.exe
	for rmax in $$(seq 3 0.5 10); do\
		mono $< >> $@ -rmin:0.1 -rmax:$$rmax -acc:0.01 -eps:0.01; \
	done


Out.conv.acc.data : mainConv.exe
	for acc in $$(seq 0.0025 0.0025 0.075); do\
		mono $< >> $@ -rmin:0.1 -rmax:8 -acc:$$acc -eps:0.01; \
	done


Out.conv.eps.data : mainConv.exe
	for eps in $$(seq 0.005 0.005 0.075); do\
		mono $< >> $@ -rmin:0.1 -rmax:8 -acc:0.01 -eps:$$eps; \
	done

Out.convergence.svg : Out.conv.rmin.data Out.conv.rmax.data Out.conv.acc.data Out.conv.eps.data
	echo '\
	set terminal svg size 1200,700 background "white" ;\
	set out "$@" ;\
	unset key ;\
	set key right top ;\
	set xlabel "rmin [Bohr radius]" ;\
	set ylabel "Energy E_0 [Hartree]" ;\
	set tics out ;\
	set grid ;\
	set title "Energy convergence with varying rmin, rmax = 8, acc = 0.01 and eps = 0.01" ;\
	set multiplot layout 2,2 rowsfirst ;\
	plot \
	"Out.conv.rmin.data" using 1:5 with points pointtype 6 notitle ;\
	set xlabel "rmax [Bohr radius]" ;\
	set title "Energy convergence with rmin = 0.1, varying rmax, acc = 0.01 and eps = 0.01" ;\
	plot \
	"Out.conv.rmax.data" using 2:5 with points pointtype 6 notitle ;\
	set xlabel "acc" ;\
	set title "Energy convergence with rmin = 0.1, rmax = 8, varying acc and eps = 0.01" ;\
	plot \
	"Out.conv.acc.data" using 3:5 with points pointtype 6 notitle ;\
	set xlabel "eps" ;\
	set title "Energy convergence with rmin = 0.1, rmax = 8, acc = 0.01 and varying eps" ;\
	plot \
	"Out.conv.eps.data" using 4:5 with points pointtype 6 notitle ;\
	' | tee log.gpi | gnuplot

Out.wave.svg : Out.wave.data 
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	set key right ;\
	set xlabel "r" ;\
	set ylabel "f" ;\
	set tics out ;\
	set grid ;\
	set title "Reduced radial wave-function for hydrogen" ;\
	f(x) = x*exp(-x) ;\
	plot \
	"$<" index 0 using 1:2 with points pointtype 6 title "Wavefunction at lowest found root energy" \
	,f(x) with lines title "Exact result x*Exp(x)";\
	' | tee log.gpi | gnuplot	


extras.dll : ../linear_equations/main.cs $(libs)/matrix/matrix.cs $(libs)/matrix/vector.cs $(libs)/genlist.cs ODE.cs roots.cs ; $(MKLIB)

main.exe : main.cs extras.dll ; $(MKEXE)

mainConv.exe : mainConv.cs extras.dll ; $(MKEXE)

clean:
	$(RM) *.exe *.dll *.log [Oo]ut* [Ll]og* [Ee]rr
