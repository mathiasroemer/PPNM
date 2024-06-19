libs = ../matlib

CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(CODE)
TIME = time --output=$@ --append --format "%e %U"

default : Out.txt Out.AB.svg Out.C.svg

Out.txt Out.data.txt : main.exe
	mono $< 1> Out.txt 2> Out.data.txt

Out.AB.svg : Out.data.txt
	echo '\
	set terminal svg size 1200,700 background "white" ;\
	set out "$@" ;\
	unset key ;\
	set key right top ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set tics out ;\
	set grid ;\
	set title "Network vs training data" ;\
	set multiplot layout 2,2 rowsfirst ;\
	plot [-1:1]\
	"$<" index 0 using 1:2 with lines lt 8 lw 2 dashtype 2 title "Linear interpolation of network values" \
	,"$<" index 0 using 1:3 with lines lt 7 lw 1 dashtype 3 title "Training function (Cos(5x-1)*Exp(-x*x))" \
	,"$<" index 1 using 1:2 with points pointtype 6 lc rgb "red" title "Training points" ;\
	set title "Linear interpolation of network vs analytical first derivative" ;\
	plot [-1:1]\
	"$<" index 2 using 1:2 with lines lt 8 lw 2 dashtype 2 title "Network" \
	,"$<" index 2 using 1:3 with lines lt 7 lw 1 dashtype 3 title "Analytical" ;\
	set title "Linear interpolation of network vs analytical second derivative" ;\
	plot [-1:1]\
	"$<" index 3 using 1:2 with lines lt 8 lw 2 dashtype 2 title "Network" \
	,"$<" index 3 using 1:3 with lines lt 7 lw 1 dashtype 3 title "Analytical" ;\
	set title "Linear interpolation of network vs numerical anti derivative" ;\
	plot [-1:1]\
	"$<" index 4 using 1:2 with lines lt 8 lw 2 dashtype 2 title "Network" \
	,"$<" index 4 using 1:3 with lines lt 7 lw 1 dashtype 3 title "Numerical" ;\
	' | tee log.gpi | gnuplot	

Out.C.svg : Out.data.txt
	echo '\
	set terminal svg size 1200,700 background "white" ;\
	set out "$@" ;\
	unset key ;\
	set key right top ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set tics out ;\
	set grid ;\
	set title "Network ODE y´´ = -5y where y(0)=0, y´(0)=1" ;\
	set multiplot layout 2,2 rowsfirst ;\
	plot \
	"$<" index 5 using 1:2 with lines lt 8 lw 2 dashtype 2 title "Network ODE prediction" \
	,"$<" index 5 using 1:3 with lines lt 7 lw 1 dashtype 3 title "Analytical solution" ;\
	' | tee log.gpi | gnuplot	


extra.dll : $(libs)/matrix/matrix.cs $(libs)/matrix/vector.cs minimization.cs QRGS.cs integration.cs; $(MKLIB)

main.exe : main.cs extra.dll; $(MKEXE)

clean:
	$(RM) *.exe *.dll *.log [Oo]ut* [Ll]og* [Ee]rr
