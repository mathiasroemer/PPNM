libs = ../matlib/matrix

CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(CODE)
TIME = time --output=$@ --append --format "%e %U"

default : Out.txt Out.convergence.svg Out.egfnc.svg

Out.egfnc.svg : Out.egfnc.data
	echo '\
	set terminal svg background "white" ;\
	set out "$@" ;\
	unset key ;\
	set key right top ;\
	set xlabel "r [Bohr radius]" ;\
	set ylabel "Probability" ;\
	set tics out ;\
	set grid ;\
	set title "Two lowest eigen-functions at rmax = 15 and dr = 0.1 vs analytic function" ;\
	plot \
	"$<" index 0 using ($$1):($$2)**2 with points pointtype 6 title "Lowest eigenvalue function" \
	,"$<" index 1 using ($$1):($$2)**2 with points pointtype 6 title "Second lowest eigenvalue function" \
	, (x*(2*exp(-x)))**2 with lines title "n = 1 reduced radial wave function" \
	, (x*((1/sqrt(2))*(1-x/2)*exp(-x/2)))**2 with lines title "n = 2 reduced radial wave function" ;\
	' | tee log.gpi | gnuplot	

Out.convergence.svg : Out.rmax.data Out.dr.data
	echo '\
	set terminal svg size 1200,700 background "white" ;\
	set out "$@" ;\
	unset key ;\
	set key right top ;\
	set xlabel "rmax [Bohr radius]" ;\
	set ylabel "Lowest energy [Hartree]" ;\
	set tics out ;\
	set grid ;\
	set title "Energy convergence with varying rmax and dr = 0.3" ;\
	set multiplot layout 1,2 rowsfirst ;\
	plot \
	"Out.rmax.data" using 1:3 with points pointtype 6 notitle ;\
	set xlabel "dr [Bohr radius]" ;\
	set title "Energy convergence with varying dr and rmax = 10" ;\
	plot \
	"Out.dr.data" using 2:3 with points pointtype 6 notitle ;\
	' | tee log.gpi | gnuplot	



Out.txt: main.exe
	mono $< > $@ -rmax:10 -dr:0.3

Out.egfnc.data: hydrogen.exe	
	mono $< > $@ -rmax:15 -dr:0.1 -egfnc:y

Out.dr.data : hydrogen.exe
	for dr in $$(seq 0.1 0.1 2); do\
		mono $< >> $@ -rmax:10 -dr:$$dr; \
	done

Out.rmax.data : hydrogen.exe
	for rmax in $$(seq 0.6 0.5 20); do\
		mono $< >> $@ -rmax:$$rmax -dr:0.3; \
	done

main.exe : main.cs matrix.dll ; $(MKEXE)
hydrogen.exe : hydrogen.cs matrix.dll ; $(MKEXE)

matrix.dll : $(libs)/matrix.cs $(libs)/vector.cs jacobi.cs; $(MKLIB)

clean:
	$(RM) *.exe *.dll [Oo]ut* [Ll]og*
