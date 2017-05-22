XML = fann.xml

all: html html-single dvi pdf ps rtf tex txt

html: fann.xml Makefile
	jw -b html -o html $(XML)
	(cd html && perl -p -i -e 's/<META/<link href="..\/style.css" rel="stylesheet" type="text\/css"><META/' *.html)

html-single:
	jw -u -b html $(XML)
	perl -p -i -e 's/<META/<link href="style.css" rel="stylesheet" type="text\/css"><META/' fann.html

dvi: fann.xml Makefile
	jw -u -b dvi $(XML)

pdf: fann.xml Makefile
	jw -u -b pdf $(XML)

ps: fann.xml Makefile
	jw -u -b ps $(XML)

rtf: fann.xml Makefile
	jw -u -b rtf $(XML)

tex: fann.xml Makefile
	jw -u -b tex $(XML)

txt: fann.xml Makefile
	jw -u -b txt $(XML)

sf: pdf html html-single
	perl -p -i -e 's/><\/BODY/> <br><div><A href="http:\/\/sourceforge.net"> <IMG src="http:\/\/sourceforge.net\/sflogo.php?group_id=93562&amp;type=5" width="210" height="62" border="0" alt="SourceForge.net Logo" \/><\/A><\/div><br> <\/body/' fann.html
	(cd html && perl -p -i -e 's/><\/BODY/> <br><div><A href="http:\/\/sourceforge.net"> <IMG src="http:\/\/sourceforge.net\/sflogo.php?group_id=93562&amp;type=5" width="210" height="62" border="0" alt="SourceForge.net Logo" \/><\/A><\/div><br> <\/body/' *.html)
	tar czf sf.tgz *.html search.php style.css fann.pdf html

spell:
	ispell -h -d american -p ./personal.dict $(XML)

clean:
	rm -rf \
		html html-single dvi pdf ps rtf tex txt \
		fann.html \
		fann.dvi \
		fann.pdf \
		fann.ps \
		fann.rtf \
		fann.tex \
		fann.txt *~

distclean: clean
