libs = ../matlib/matrix

CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:, $(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(CODE)
TIME = time --output=$@ --append --format "%e %U"

default : Out.txt 


Out.txt : main.exe
	mono $< > $@

main.exe : main.cs matrix.dll ; $(MKEXE)

matrix.dll : $(libs)/matrix.cs $(libs)/vector.cs ; $(MKLIB)

clean:
	$(RM) *.exe *.dll [Oo]ut* [Ll]og*
