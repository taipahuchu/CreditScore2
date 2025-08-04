/* ******************************************************************* */
/*        This is a script that explores LGD calculations              */
/*        Done by Taipa                                                */
/*        December 2024                                                */
/*        Methods 1) Linear regression,                                */
/*                2) Fractional response using many link functions,    */
/*                3) Inflated Beta regression                          */
/*                4) Mixed Effects Method                              */
/*                5) Logistic regression                              */
/* *********************************************************************/
Options validvarname=V7;

/*                1.DATA  COLLECTION                   */
/* assigning libraries to the folder which has data */
libname lgd '/home/u49308301/Research/LGD';
%let path = %sysfunc(pathname(lgd));
%put &path;

/* Assigning file names to the files with the loan data  */
FILENAME train "/home/u49308301/Research/LGD/data.csv";

/* Importing the files from csv  to a sas data set */
PROC IMPORT DATAFILE=train dbms=csv replace out=lgd;
	getnames=Yes;
RUN;

/* Marking the response variable between 0 and 1 and making it look like LGD data */
 
data LGD;
SET LGD;
if target <= 6.5 then target = 6.5;
if target >=  8.5 THEN TARGET = 8.5 ;
target1 = (target- 6.5)/(8.5 - 6.5);
RUN;

PROC UNIVARIATE DATA =lgd noprint ; 
HISTOGRAM target;
RUN;

PROC UNIVARIATE DATA =lgd noprint ; 
HISTOGRAM target1;
RUN;

title "Splitting the data into training and testing datasets";
PROC SURVEYSELECT noprint data=work.LGD out=LGD 
		samprate=.7 seed=1234 outall;
run;

/* verify sampling */
proc freq data=LGD;
	table selected/ nocum nopercent;
run;

/* dropping unwanted columns produced by the surveyselect procedure and the identity column and target*/
data train (drop=selected SelectionProb SamplingWeight) test (drop=selected 
		SelectionProb SamplingWeight);
	set LGD;

	if selected=1 then
		output train;
	else
		output test;
		drop id target;
		RENAME target1=LGD;
run;



/*            Program 1.      LINEAR REGRESSION               */

proc reg data=train outest=train_model;
	model LGD =cont1--cont14;
quit;

proc score data=test score=train_model out=reg_output (keep=LGD
		model1 ) predict type=parms;
	var cont1--cont14;
run;

/* Plot the predicted data */
PROC UNIVARIATE DATA=reg_output noprint;
	histogram LGD / odstitle=title;
run;

title;
title 'This is a plot of the predicted values removing the outlies';

PROC UNIVARIATE DATA=reg_output NOPRINT;
	WHERE MODEL1 between 0 and 1;
	HISTOGRAM MODEL1/ odstitle=title;
RUN;

title;

data reg_output;
set reg_output ;
rename model1=Pred;
RUN;


/*            Program 2.      FRACTIONAL RESPONSE REGRESSION               */


proc nlmixed data=train tech=newrap maxiter=3000 maxfunc=3000 qtol=0.0001;
	parms b0-b14=0.0001;
	cov_mu=b0+b1*Cont1+b2*Cont2+…+b14*Cont14;
	mu=logistic(cov_mu);
	loglikefun=LGD*log(mu)+(1-LGD)*log(1-mu);
	model LGD~general(loglikefun);
	predict mu out=frac_resp_output (keep=instrument_id LGD pred);
run;


/*            Program 3.       Inflated beta regression                   */

proc nlmixed data=train tech=quanew maxiter=3000 maxfunc=3000 qtol=0.0001;
	parms b0-b14=0.0001 pie=0.2 kesai=0.3 phi=2;
	cov_mu=b0+b1*Cont1+b2*Cont2+…+b14*Cont14;
	mu=logistic(cov_mu);

	if LGD=0 then
		loglikefun=log(pie)+log(1-kesai);

	if LGD>=1 then
		loglikefun=log(pie)+log(kesai);

	if 0<LGD<1 then
		loglikefun=log(1-pie)+lgamma(phi)-lgamma(mu*phi)-lgamma((1-mu)*phi) 
			+(mu*phi-1)*log(LGD)+((1-mu)*phi-1)*log(1-LGD);
	predict pie*kesai+(1-pie)*mu out=Inf_beta_output (keep=instrument_id LGD pred);
	model LGD~general(loglikefun);
run;


/* Putting all the output table names into a table   */
DATA tables;
input col $32.;
datalines;
frac_resp_output
Inf_beta_output
reg_output
;
run;

/* Putting all names into a macro */

data _null_;
set TABLES end = eof;
CALL SYMPUTX(CAT("Col",_n_),COL);
if eof then CALL SYMPUTX('myobs',_n_);
RUN;


%put &myobs;



%MACRO myrun;
%do i= 1 %to &myobs;
/* Performance metrics  */
%let output_data=&&Col&i.;
title "Performance metrics for &&Col&i";
proc iml;
	use &output_data;
	read all var {LGD, pred} into data;
	close &output_data;
	m=nrow(data);
	LGD=data[, 1];
	LGD_estimate=data[, 2];
	start var(x);
	mean=x[:, ];
	countn=j(1, ncol(x));

	do i=1 to ncol(x);
		countn[i]=sum(x[, i]^=.);
	end;
	var=(x-mean)[##, ]/(countn-1);
	return (var);
	finish;
	resid=LGD_estimate-LGD;
	SS_error=resid`*resid;
	SS_total=var(LGD)#(m-1);
	R_square=1-SS_error/SS_total;
	RMSE=sqrt(resid`*resid/m);
	MAE=sum(abs(resid))/m;
	print R_square RMSE MAE;
	create output_measure var {R_square RMSE MAE};
	append;
	quit;

%end;
title;
%MEND myrun;

%myrun



















