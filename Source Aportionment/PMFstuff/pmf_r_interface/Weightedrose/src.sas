data src2;
	infile "F:/EPA/Weightedrose/src2.txt" firstobs=2;
	input date: IS8601DA. lead fireworks copper sumsec zinc soil steel winsec mobile;
run;

data winddir;
	infile "F:/EPA/Weightedrose/winddir2.csv" dlm=',';
	input date: MMDDYY10. winddir;
run;

proc sort data=src2 out=src2;
	by date;
run;

proc sort data=winddir out=winddir;
	by date;
run;


data src3;
	merge src2 winddir;
	by date;
	if lead='.' then delete;
	format date IS8601DA.;
	file 'F:\EPA\Weightedrose\src3.csv' dlm=',';
	put date: IS8601DA. winddir lead fireworks copper sumsec zinc soil steel winsec mobile;
run;

proc print;
run;
