CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(CODE)

default : Out.SignSpline.svg Out.CosDerivAndInteg.svg Out.SecondDerivative.svg

Out.SecondDerivative.svg : Out.Splines.data
	echo '\
		set terminal svg size 1100,600 background "white" ;\
		set out "$@" ;\
		set key left ;\
		set key box width -8 ;\
		set xlabel "x";\
		set ylabel "y";\
		set grid ;\
		set title "Subspline of Sign(x) using order 3";\
		set multiplot layout 1,2 rowsfirst ;\
		plot \
		"$<" index 1 using 1:2 with lines lc rgb "black" title "Function [Sign(x)]" \
		,"$<" index 2 using 1:2 pointtype 1 title "Subspline interpolation" \
		,"$<" index 9 using 1:2 pointtype 6 ps 0.5 title "First derivative" \
		,"$<" index 10 using 1:2 pointtype 6 ps 0.5 title "Second derivative" \
		,"$<" index 0 using 1:2 pointtype 7 ps 0.7 lc rgb "black" title "Data points" ;\
		set title "Subspline of Sign(x) using order 4";\
		plot \
		"$<" index 1 using 1:2 with lines lc rgb "black" title "Function [Sign(x)]" \
		,"$<" index 11 using 1:2 pointtype 1 title "Subspline interpolation" \
		,"$<" index 12 using 1:2 pointtype 6 ps 0.5 title "First derivative" \
		,"$<" index 0 using 1:2 pointtype 7 ps 0.7 lc rgb "black" title "Data points" ;\
	'|tee log.gpi |gnuplot


Out.SignSpline.svg : Out.Splines.data
	echo '\
		set terminal svg background "white" ;\
		set out "$@" ;\
		set key left ;\
		set key box width -9 ;\
		set xlabel "x";\
		set ylabel "y";\
		set grid ;\
		set title "Subspline interpolation of Sign(x)";\
		plot \
		"$<" index 1 using 1:2 with lines lc rgb "black" title "Function [Sign(x)]" \
		,"$<" index 2 using 1:2 pointtype 1 title "Subspline interpolation" \
		,"$<" index 3 using 1:2 pointtype 6 ps 0.5 title "Cubic spline interpolation" \
		,"$<" index 0 using 1:2 pointtype 7 ps 0.7 lc rgb "black" title "Data points" ;\
	'|tee log.gpi |gnuplot

Out.CosDerivAndInteg.svg : Out.Splines.data
	echo '\
		set terminal svg background "white" ;\
		set out "$@" ;\
		set key right ;\
		set key bottom ;\
		set key box width -12 ;\
		set xlabel "x";\
		set ylabel "y";\
		set grid ;\
		set title "Subspline interpolation of Sign(x)";\
		plot \
		"$<" index 5 using 1:2 with lines lc rgb "black" title "Function [Cos(x)]" \
		,"$<" index 6 using 1:2 pointtype 1 title "Subspline interpolation" \
		,"$<" index 7 using 1:2 pointtype 6 ps 0.5 title "Subspline first derivative [-Sin(x) ?]" \
		,"$<" index 8 using 1:2 pointtype 6 ps 0.5 title "Subspline integral [Sin(x) ?]" \
		,"$<" index 4 using 1:2 pointtype 7 ps 0.7 lc rgb "black" title "Data points" ;\
	'|tee log.gpi |gnuplot


Out.Splines.data : main.exe
	mono $< > $@

extras.dll : bettersubspline.cs subspline.cs cspline.cs vector.cs ; $(MKLIB)

main.exe : main.cs extras.dll ; $(MKEXE)

clean:
	$(RM) *.exe *.dll *.log [Oo]ut* [Ll]og* [Ee]rr
