CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(CODE)
TIME = time --output=$@ --append --format "%e %U"

default : Out.SignSpline.svg Out.CosDerivAndInteg.svg

Out.SignSpline.svg : Out.Splines.data
	echo '\
		set terminal svg background "white" ;\
		set out "$@" ;\
		set key left ;\
		set key box width -8 ;\
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
		set key left ;\
		set key box width -12 ;\
		set xlabel "x";\
		set ylabel "y";\
		set grid ;\
		set title "Subspline interpolation of Sign(x)";\
		plot \
		"$<" index 5 using 1:2 with lines lc rgb "black" title "Function [Cos(x)]" \
		,"$<" index 6 using 1:2 pointtype 1 title "Subspline interpolation" \
		,"$<" index 7 using 1:2 pointtype 6 ps 0.5 title "Subspline derivative [-Sin(x) ?]" \
		,"$<" index 8 using 1:2 pointtype 6 ps 0.5 title "Subspline integral [Sin(x) ?]" \
		,"$<" index 4 using 1:2 pointtype 7 ps 0.7 lc rgb "black" title "Data points" ;\
	'|tee log.gpi |gnuplot


Out.Splines.data : main.exe
	mono $< > $@

extras.dll : subspline.cs cspline.cs vector.cs ; $(MKLIB)

main.exe : main.cs extras.dll ; $(MKEXE)

clean:
	$(RM) *.exe *.dll *.log [Oo]ut* [Ll]og* [Ee]rr
