

ods html close;
ods html;

/* &folder contains the 4ti2 input file (dimension d=3,4,5,6) and output files (dimension d=3,4,5) */
%let folder=C:\Users\roberto\Politecnico di Torino Staff Dropbox\Roberto Fontana\DOE\algstat26\software\4ti2;

/* the lib "b" connects to SAS datasets which contains counting vectors of 
strength 2 non-isomorphic OAs generated using the 
oa package available at http://www.pietereendebak.nl/oapackage/series.html */ 
libname b "C:\Users\roberto\Politecnico di Torino Staff Dropbox\Roberto Fontana\DOE\algstat26\software\oas";

/* &caso must must be equal to x12 (1st and 2nd order moments) */
%let caso=x12;

/* dimension dd can be chosen between 3 and 6 */
*%let dd=3;
%let dd=4;
*%let dd=5;
%let dd=6;

/* strength must must be equal to 2 */
%let strength=2;

/* some string of variable names are generated */
%let nfattori=%eval(&dd);
%put &nfattori;
%global nomivars nomivars1 y time nomivars2 xxcol exxcol;
%macro nomivar(nf);
%put "NFATTORI" &nf.;
%let p=;
%let s=;
%do i=1 %to &nf.;
%let p=&p p&i.;
%let s=&s s&i.;
%end;
%let y=;
%let m=;
%let time=;
%let nd=;
%let nff=%eval(2**&nf.);
%let xxcol=;
%let exxcol=permfirst &p &s;
%do i=1 %to &nff.;
%let y=&y y&i.;
%let m=&m m&i.;
%let time=&time t&i.;
%let nd=&nd nd&i.;
%let xxcol=&xxcol col&i.;
%let exxcol=&exxcol col&i.;
%end;
%let nomivars=noa permfirst idpc ids &p &s &y &m &time;
%let nomivars1=noa permfirst idpc ids ndx &m;
%let nomivars2=&nomivars &nd ndx;
%mend;

%nomivar(&nfattori);
%put &nomivars;
%put &nomivars1;
%put &y;
%put &time;
%put &xxcol;
%put &exxcol;



proc iml;

/* START OF MODULE DEFINITIONS */
start ffmix(lev);
nfac=ncol(lev);          /*number of factors*/
npt=1;
do i=1 to nfac;
npt=npt*lev[i];          /* cardinality of the design*/ 
end;

ff=j(npt,nfac,0);
do i=0 to npt-1;
tmp=i;
	do j=1 to nfac;
	ff[i+1,nfac-j+1]= mod(tmp,lev[nfac-j+1]);   /*conversion base algorithm*/
    tmp=int(tmp/lev[nfac-j+1]);	
	end;
end;
return(ff);
finish;


/* END OF MODULE DEFINITIONS */

/********************* START OF THE MAIN **********/
d=&dd.;
strength=&strength.;
print "DIMENSION" d "STRENGTH" strength;

/* the moments of order 1,...,d associated to the pmf p are computed as M*p 
where M is the 2^d * 2^d moment matrix. For d=3
M= 
1 1 1 1 1 1 1 1 
0 1 0 1 0 1 0 1 
0 0 1 1 0 0 1 1 
0 0 0 1 0 0 0 1 
0 0 0 0 1 1 1 1 
0 0 0 0 0 1 0 1 
0 0 0 0 0 0 1 1 
0 0 0 0 0 0 0 1 

M=can be cmputed as M=M0@M0@M0 where is the kronecker product and M1 is
M1=
1 1
0 1

the matrix M (M1) is denoted by mom (mom0) 
*/
mom0={1 1, 0 1};
mom=mom0;
do i=2 to d;
mom=mom@mom0;
end;
print mom;


/* AFTER computing the extreme rays with 4ti2 */
/* For d <= 5, 4ti2 can compute the complete set of extreme rays */
/* The resulting 4ti2 output file is then read */

if d<=5 then do;
	nomefile="H_x12_d_&dd._t_&strength..mat.ray";
	print nomefile;


	submit nomefile;

			data n_ray;
			infile  "&folder\&nomefile"
			firstobs=1 obs=1;
			input nr nc;
			RUN;

	endsubmit;
	use n_ray;
	read all var{nr} into nr;
	close n_ray;
	/* number of rays found */
	print "NUMBER OF RAYS" nr;


	nc=2**d;
	print nc nomefile;
	submit nc nomefile;

					data h_ray;
						infile  "&folder\&nomefile"
						firstobs=2;
						input x1-x&nc;
					RUN;

	endsubmit;
	use h_ray;
	read all into h_ray;
	close h_ray;

	if d<=4 then print "THE FULL SET OF RAYS" h_ray;

end;

/* for d=6 4ti2 can not compute the full set of rays. The OAs of size 8,12,16,20,24 are used */
if d=6 then do;
	use b.cnt_2_6_t2_8;
	read all into noniso1;
	close b.cnt_2_6_t2_8;

	use b.cnt_2_6_t2_12;
	read all into noniso2;
	close b.cnt_2_6_t2_12;
	
	use b.cnt_2_6_t2_16;
	read all into noniso3;
	close b.cnt_2_6_t2_16;

	use b.cnt_2_6_t2_20;
	read all into noniso4;
	close b.cnt_2_6_t2_20;

	use b.cnt_2_6_t2_24;
	read all into noniso5;
	close b.cnt_2_6_t2_24;

	noniso=noniso1//noniso2//noniso3//noniso4//noniso5;
	noniso=t(noniso);
	n_noniso=ncol(noniso);

end;

/* the matrix er contains, as columns, 
- for d=3,4,5  the full set of extremal rays :
- for d=6 the counting vectors of OAs of size 8,12,16,20,24 (d=6) */
if d<=5 then er=t(h_ray);
if d=6 then er=noniso;

print "THERE ARE " (ncol(er)) "AVAILABLE PMFS ";
/* the rays are normalized */	
do i=1 to ncol(er);
		er[,i]=er[,i]/er[+,i];
end;

/* for d=6 it can happen that different OAs give the same pmf. Such duplicated pmfs 
are eliminated */
if d=6 then do;
	ter=t(er);
	create ter from ter;
	append from ter;
	close ter;

	submit;
		proc sort data=ter nodupkey;
		by &xxcol;
		run;
	endsubmit;
	use ter;
	read all into ter;
	close ter;
	er=t(ter);
end;
print "THEre ARE " (ncol(er)) "AVAILABLE (AFTER ELIMINATION OF REPEATED PMFS)";	

/* only "small" cases are printed */
if d<=4 then print er;


/* computation of the full set of moments corresponding to the extreme pmfs */
moments=mom*er;
if d<=4 then do;
	print er, moments;
end;
print "END OF MOMENTS COMPUTATIONS FOR THE AVAILABLE PMFS"; 

/* ff contains the full factorial designs {0,1}^d */
/* it is used to represents all the moments
e.g. 000 --> E[1], 001 --> E[X3], 101 --> E[X1*X3]
*/

ff=ffmix(j(1,d,2));
if d<=4 then print ff;

/* for all the moments (E[X1], E[X2], E[X1*X2], ecc.) the minimum and the maximum 
of the moments over the set of all the available pmfs are computed */ 
minmom=j(nrow(moments),1,.);
maxmom=j(nrow(moments),1,.);
ord=j(nrow(moments),1,.);
free labels;
do i=1 to nrow(moments);
	ord[i]=ff[i,+];
	labels=labels//(ff[i,]);
	xmin=min(moments[i,]);
	xmax=max(moments[i,]);
	minmom[i]=xmin;
	maxmom[i]=xmax;
end;
/* labels are the moments, ord their order, minmom and maxmom the minimum (maximum) of
the moments for all the pmfs in the matrix er */
print labels ord, minmom maxmom;

/* For moments of all orders, the minimum and maximum values are computed */
free best;
do i=1 to d;
	righe=loc(ord=i);
	xord=ord[righe];
	xminmom=minmom[righe];
	xmaxmom=maxmom[righe];
	print i xord xminmom xmaxmom;

	best=best//(i||min(xminmom)||max(xmaxmom));
end;
print best;
print "END OF THE COMPUTATION OF THE BOUNDS OF THE MOMENTS";

/* effect of switching on the moments of the pmfs in er */
/* the switch matrix contains all the possible switch, (\sigma vectors in the paper)*/
switch=ffmix(j(1,d,2));

/* for all the the new moments are computed. This are stored in the matrix srtmom */
free bests;
do k=1 to nrow(switch);
	supp=mod(ff+switch[k,],2);
	supper=supp||er;
	call sortndx(ndx, supper, 1:d);
	sorted = supper[ndx,];
	newer=sorted[,d+1:ncol(sorted)];
	srtmom=mom*newer;

	minmom=j(nrow(srtmom)-1,1,.);
	maxmom=j(nrow(srtmom)-1,1,.);
	ord=j(nrow(srtmom)-1,1,.);
	do i=2 to nrow(srtmom);
		ord[i-1]=ff[i,+];
		xmin=min(srtmom[i,]);
		xmax=max(srtmom[i,]);

		minmom[i-1]=xmin;
		maxmom[i-1]=xmax;
	end;

	*print ord minmom maxmom;

	free best;
	do i=1 to d;
		righe=loc(ord=i);
		xord=ord[righe];
		xminmom=minmom[righe];
		xmaxmom=maxmom[righe];
		*print i xord xminmom xmaxmom;

		best=best//(i||min(xminmom)||max(xmaxmom));
	end;
	*print best;
	bests=bests//(j(nrow(best),1,k)||best);
end;
/* matrix bests
in the first column the id of the switch is reported 
		e.g. 1 --> \sigma=(0,...,0,0) i.e. no switching
			 2 --> \sigma=(0,...,0,1) i.e only variable Xd is switched
			...
			2^d --> \sigma=(1,\ldots,1,1) i.e. all the variables are switched
in the second column the order of the moments (from 1 to d)
in the third (fourth) column the minimum (maximum) of the moment over the pmfs that
are obtainded as the switch of the pmfs of er */
print bests;

/* The matrix allbest contains the minimum (maximum) values of the moments over all switches */
free allbest;
do i=1 to d;
	righe=loc(bests[,2]=i);
	xminmom=bests[righe,3];
	xmaxmom=bests[righe,4];
	print i xminmom xmaxmom;
			allbest=allbest//(i||min(xminmom)||max(xmaxmom));
end;
print allbest;
print "A 1";

/* END OF COMPUTATION OF THE BOUNDS OF THE MOMENTS
/* For d = 3, 4, 5, all extreme rays are available, and switching does not provide tighter bounds */
/* This part of the code can be used for verification */

/* For d = 6, only the pmfs corresponding to the selected OAs are available */
/* In this case, this part of the code can be used to improve the bounds */


/* FURTHER COMPUTATIONS, mainly for d = 6 */
/* Generation of the H matrix used as input for 4ti2 */

tt=strength;
ff=ffmix(j(1,d,2));
xx=1-2*ff;
free h ords xcols;
do i=2 to nrow(ff);
	ord=ff[i,+];
	if ord<=tt then do;
		ords=ords//ord;
		cols=loc(ff[i,]);
		tmp=j(1,d,.);
		tmp[1:ncol(cols)]=cols;
		xcols=xcols//tmp;
		uno=j(nrow(ff),1,1);
		do j=1 to ncol(cols);
			uno=hdir(uno,xx[,cols[j]]);
		end;
		h=h||uno;
	end;
end;
hbase=t(h);
print ords, xcols hbase;
/* end of H matrix computation*/

/* Identification of the extreme pmfs among those stored in the er matrix */
free ranks nonisoray;
/* The matrix er is copied into the matrix noniso */
/* Note: the name noniso is appropriate for d = 6; for d = 3, 4, 5, it is not,
   since er contains the complete set of extreme pmfs */
noniso=er;
israyg=j(1,ncol(noniso),.);
do i=1 to ncol(noniso);
	exray=noniso[,i];
	/* for each p=exray, the AA matrix is the submatrix of H_d^p containing the columns of H_d
   	corresponding to the positive entries of p */
	AA=hbase[,loc(exray)];
	rank_AA = round(trace(ginv(AA) * AA));
	isray=0;
	/* p = exray is an extreme ray if rank(AA) = ncol(AA) - 1 */
	if ncol(AA)-rank_AA=1 then isray=1;
	*print rank_AA (ncol(AA));
	ranks=ranks//(rank_AA||ncol(AA)||isray);
	/* the matrix nonisoray contains the extreme rays */
	if isray=1 then nonisoray=nonisoray||exray;
	israyg[i]=isray;
end;
n_nonisoray=ranks[+,3];
/* n_nonisoray: number of extreme pmfs */
print n_nonisoray (israyg[+]);
if d<=4 then print ranks nonisoray;
chk=hbase*nonisoray;
if d<=4 then print chk;

/* computation of the moments over the extreme pmfs */
moments=mom*nonisoray;

if d<=4 then print moments;

ff=ffmix(j(1,d,2));
if d<=4 then print ff;

minmom=j(nrow(moments),1,.);
maxmom=j(nrow(moments),1,.);
ord=j(nrow(moments),1,.);
do i=1 to nrow(moments);
	ord[i]=ff[i,+];
	xmin=min(moments[i,]);
	xmax=max(moments[i,]);

	minmom[i]=xmin;
	maxmom[i]=xmax;
end;

print ord minmom maxmom;

do i=1 to d;
righe=loc(ord=i);
xord=ord[righe];
xminmom=minmom[righe];
xmaxmom=maxmom[righe];
print i xord xminmom xmaxmom;
end;
print "END B";

free best;
do i=1 to d;
righe=loc(ord=i);
xord=ord[righe];
xminmom=minmom[righe];
xmaxmom=maxmom[righe];
best=best//(i||min(xminmom)||max(xmaxmom));
end;
print best;
print "END B.1";

/* all the 2^d switches are applied to the available extreme pmfs which 
are stored in nonisoray*/
free bests;
do k=1 to nrow(switch);
	supp=mod(ff+switch[k,],2);
	supper=supp||nonisoray;
	call sortndx(ndx, supper, 1:d);
	sorted = supper[ndx,];
	newer=sorted[,d+1:ncol(sorted)];
	srtmom=mom*newer;

	minmom=j(nrow(srtmom)-1,1,.);
	maxmom=j(nrow(srtmom)-1,1,.);
	ord=j(nrow(srtmom)-1,1,.);
	do i=2 to nrow(srtmom);
		ord[i-1]=ff[i,+];
		xmin=min(srtmom[i,]);
		xmax=max(srtmom[i,]);

		minmom[i-1]=xmin;
		maxmom[i-1]=xmax;
	end;

	*print ord minmom maxmom;

	free best;
	do i=1 to d;
		righe=loc(ord=i);
		xord=ord[righe];
		xminmom=minmom[righe];
		xmaxmom=maxmom[righe];
		*print i xord xminmom xmaxmom;

		best=best//(i||min(xminmom)||max(xmaxmom));
	end;
	*print best;
	bests=bests//(j(nrow(best),1,k)||best);
end;
print bests;

free allbest;
do i=1 to d;
	righe=loc(bests[,2]=i);
	xminmom=bests[righe,3];
	xmaxmom=bests[righe,4];
	print i xminmom xmaxmom;
			allbest=allbest//(i||min(xminmom)||max(xmaxmom));
end;
print allbest;
/* END OF FURTHER COMPUTATIONS*/


