Out.txt: main.exe input.txt Makefile
	mono main.exe < input.txt > $@

main.exe: main.cs genlist.dll
	mcs -target:exe -out:$@ -reference:genlist.dll main.cs

genlist.dll: genlist.cs
	mcs -target:library -out:$@ genlist.cs

clean:
	rm -f main.exe Out.txt genlist.dll

