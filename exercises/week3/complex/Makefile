Out.txt: main.exe
	mono > $@ $<

main.exe : main.cs cmath.dll
	mcs -reference:cmath.dll -target:exe -out:main.exe main.cs

cmath.dll : cmath.cs complex.cs
	mcs -target:library -out:./cmath.dll $^

clean:
	rm -f main.exe Out.txt
