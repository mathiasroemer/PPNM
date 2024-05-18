libs = ../matlib/matrix

CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(CODE)
TIME = time --output=$@ --append --format "%e %U"

default : Out.LinSpline.svg Out.QSpline.svg Out.CSpline.svg Out.CSplineAndBuiltIn.svg


Out.LinSpline.svg : OutLinSpline.txt
	echo '\
		set terminal svg background "white" ;\
		set out "$@" ;\
		set xlabel "x";\
		set ylabel "y";\
		set grid ;\
		set title "Linear interpolation of Cos(x)";\
		plot "$<" index 0 using 1:2 pointtype 6 title "Data" \
		,"$<" index 1 using 1:2 with lines title "Function [Cos(x)]" \
		,"$<" index 2 using 1:2 pointtype 1 title "Linear interpolation" \
		,"$<" index 3 using 1:2 pointtype 6 ps 0.5 title "Integral [Sin(x) ?]" ;\
	'|tee log.gpi |gnuplot

Out.QSpline.svg : OutQSpline.txt
	echo '\
		set terminal svg background "white" ;\
		set out "$@" ;\
		set xlabel "x";\
		set ylabel "y";\
		set grid ;\
		set title "Quadratic interpolation of Sin(x)";\
		plot "$<" index 0 using 1:2 pointtype 6 title "Data" \
		,"$<" index 1 using 1:2 with lines title "Function [Sin(x)]" \
		,"$<" index 2 using 1:2 pointtype 1 title "Quadratic interpolation" \
		,"$<" index 3 using 1:2 pointtype 6 ps 0.5 lc rgb "red" title "Derivative [Cos(x) ?]" \
		,"$<" index 4 using 1:2 pointtype 6 ps 0.5 lc rgb "black" title "Integral [-Cos(x) ?]" ;\
	'|tee log.gpi |gnuplot


Out.CSpline.svg : OutCSpline.txt
	echo '\
		set terminal svg background "white" ;\
		set out "$@" ;\
		set xlabel "x";\
		set ylabel "y";\
		set grid ;\
		set title "Cubic interpolation of Sin(x)";\
		plot "$<" index 0 using 1:2 pointtype 6 title "Data" \
		,"$<" index 1 using 1:2 with lines title "Function [Sin(x)]" \
		,"$<" index 2 using 1:2 pointtype 1 title "Cubic interpolation" \
		,"$<" index 3 using 1:2 pointtype 6 ps 0.5 lc rgb "red" title "Derivative [Cos(x) ?]" \
		,"$<" index 4 using 1:2 pointtype 6 ps 0.5 lc rgb "black" title "Integral [-Cos(x) ?]" ;\
	'|tee log.gpi |gnuplot


Out.CSplineAndBuiltIn.svg : OutCSpline.txt
	echo '\
		set terminal svg background "white" ;\
		set out "$@" ;\
		set xlabel "x";\
		set ylabel "y";\
		set grid ;\
		set title "Gnuplot built-in cubic splines vs our implentation";\
		plot "$<" index 5 using 1:2 pointtype 6 title "Data" \
		,"$<" index 6 using 1:2 with lines title "Function [Sin(2*x)*Exp(-x/3.0)+Cos(x)]" \
		,"$<" index 7 using 1:2 pointtype 1 title "Cubic interpolation" \
		,"$<" index 5 using 1:2 smooth csplines title "Gnuplot Cubic splines" ;\
	'|tee log.gpi |gnuplot


OutLinSpline.txt : linspline.exe
	mono $< > $@


OutQSpline.txt : qspline.exe
	mono $< > $@

OutCSpline.txt : cspline.exe
	mono $< > $@

vector.dll : $(libs)/vector.cs ; $(MKLIB)

linspline.exe : linspline.cs vector.dll ; $(MKEXE)

qspline.exe : qspline.cs vector.dll ; $(MKEXE)

cspline.exe : cspline.cs vector.dll ; $(MKEXE)

clean:
	$(RM) *.exe *.dll *.log [Oo]ut* [Ll]og* [Ee]rr
