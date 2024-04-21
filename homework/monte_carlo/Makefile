libs = ../matlib

CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(CODE)
TIME = time --output=$@ --append --format "%e %U"

default : Out.txt Out.aoc.svg Out.error_scaling.svg

Out.txt Out.data.txt : main.exe
	mono $< 1>Out.txt 2>Out.data.txt

Out.aoc.svg : Out.data.txt 
	echo '\
	set terminal svg size 1200,500 background "white" ;\
	set out "$@" ;\
	set key left ;\
	set xlabel "N" ;\
	set ylabel "Area" ;\
	set tics out ;\
	set grid ;\
	set title "Area of unit circle" ;\
	set multiplot layout 1,2 rowsfirst;\
	f(x)=1/sqrt(x);\
	plot \
	"$<" index 0 using 1:2 with points pointtype 6 title "Area with N points" ;\
	set ylabel "Error" ;\
	set title "Error of area of unit circle" ;\
	plot \
	f(x) title "1/Sqrt(x)" \
	,"$<" index 0 using 1:3 with points pointtype 6 title "Estimated error" \
	,"$<" index 0 using 1:4 with points pointtype 2 title "Exact error" ;\
	' | tee log.gpi | gnuplot	

Out.error_scaling.svg : Out.data.txt
	echo '\
	set terminal svg size 1200,500 background "white" ;\
	set out "$@" ;\
	set key left ;\
	set xlabel "N" ;\
	set ylabel "Error" ;\
	set tics out ;\
	set grid ;\
	set title "Error scaling of area of unit circle" ;\
	set multiplot layout 1,2 rowsfirst;\
	plot \
	"$<" index 0 using 1:3 with points pointtype 6 title "Estimated error of Plain mc" \
	,"$<" index 1 using 1:3 with points pointtype 6 title "Estimated error of Quasi mc" ;\
	plot \
	"$<" index 2 using 1:2 with points pointtype 2 title "Exact error of Plain mc" \
	,"$<" index 2 using 1:3 with points pointtype 2 title "Exact error of Quasi mc" ;\
	' | tee log.gpi | gnuplot	



vector.dll : $(libs)/matrix/vector.cs ; $(MKLIB)

main.exe : main.cs vector.dll ; $(MKEXE)

clean:
	$(RM) *.exe *.dll *.log [Oo]ut* [Ll]og* [Ee]rr
