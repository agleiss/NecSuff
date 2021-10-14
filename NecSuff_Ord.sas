*******************************************************************;
*
* NECSUFF_ORD
* ===========
*
* SAS-macro for computation of degrees of necessity and sufficiency
* for an ordinal outcome variable (assuming proportional odds)
*
* Gleiss, A., Henderson, R., Schemper, M., Degrees of necessity 
* and of sufficiency: further results and extensions, with an 
* application to covid-19 mortality in Austria
* accepted by Statistics in Medicine
*
* Author:  Andreas Gleiss
* Version: 1.1 (freq option)
* Date:    06 Oct 2020
*
* Macro parameters:
* =================
* 
* data		name of SAS data set
* y			name of ordinal outcome variable with values from ylist 
* ylist		outcome codes sorted from "worst" to "best" outcome
*			(lowest code corresponds to "best" outcome)
* x			names of independent variables (or linear predictor)
* xclass	list of independent class variable names
* freq		name of optional frequency variable (default: none)
* odssel	control output from proc logistic 
*			(default: ParameterEstimates OddsRatios)
* print 	=1 to print results (default), =0 (e.g., for simulations)
*
*******************************************************************;

%let mypath=E:\Andreas\Makros & Vorlagen\SAS-Macros; * change to path containing NecSuff.SAS;
%include "&mypath\NecSuff.SAS";

%macro necsuff_ord(data=, y=, ylist=2 1 0, x=, xclass=, freq=,
					odssel=ParameterEstimates OddsRatios, print=1);
%let ny=0;
%do %while(%scan(&ylist,&ny+1)~=);
	%let ny=%eval(&ny+1);
	%let y&ny=%scan(&ylist,&ny);
	%end;

proc format;
	value cat -999="W.sum";
	run;

proc freq data=&data noprint;
	tables &y / out=_tab outcum;
	%if &freq~= %then weight &freq;;
	run;
proc datasets noprint;
	delete _all;
	run;
data _works;
	set &data;
	_zeil=_n_;
	%if &freq= %then freq=1;
			   %else freq=&freq;;
	format &y f3.0; * quasi unformat;
	run;
ods select &odssel;
proc logistic data=_works;
	class &y &xclass / order=internal;
	model &y(descending)=&x; 
	freq freq;
	output out=_estpo pred=p_i_dach xbeta=xbeta;
	run;

%do iy=1 %to %eval(&ny-1);
	data _po&iy;
		set _estpo;
		where _level_=&&y&iy;
		if not missing(&y) then _yd=(&y >= &&y&iy); * since lowest code = best outcome;
		run;
	data _po_in&iy;
		set _po&iy(where=(not missing(_yd) and not missing(xbeta)));
		run;
	%necsuff(data=_po&iy, y=_yd, x=xbeta, freq=&freq, refcat=first, 
			 inpred=_po_in&iy, inpredvar=p_i_dach, 
			 odssel=none, print=0);

	data _erg;
		set _necsuff_;
		cat=&&y&iy;
		keep cat dn: ds: ev: or p_bar alpha;
		run;
	proc append data=_erg base=_all force;
		run;
	%end;
proc means data=_tab(where=(&y>&&y&ny)) noprint;
	var count;
	output out=_tabsum sum=cum_freq;
	run;
data _tab_;
	set _tab(rename=(&y=cat) where=(cat>&&y&ny) drop=cum:);
	if _n_=1 then set _tabsum(keep=cum_freq); output;
	run; 
proc sort data=_all;
	by cat;
	run;
data _all_;
	merge _all
		  _tab_;
	by cat;
	weight=count/cum_freq;
	dn1_term=dn_1*weight;
	ds1_term=ds_1*weight;
	dn2_term=dn_2*weight;
	ds2_term=ds_2*weight;
	ev_term=ev_indir*weight;
	run;
proc means data=_all_ sum noprint;
	var ev_term dn1_term ds1_term dn2_term ds2_term;
	output out=_wsums sum=;
	run; * auch EV als (gewichtetes) Mittel aus Einzel-EVs;
data _all2;
	set _all_
		_wsums(in=wsum rename=(dn1_term=dn_1 ds1_term=ds_1 
							   dn2_term=dn_2 ds2_term=ds_2 ev_term=ev_indir));
	if wsum then cat=-999;
	format cat cat.;
	run;
%if &print=1 %then %do;
	title "Degrees of necessity and sufficiency";
	title2 "Ordinal outcome = &y";
	proc print data=_all2 noobs label;
		var cat p_bar alpha or ev_indir  
			dn_1 ds_1 dn_2 ds_2 weight;
		format or p_bar alpha ev_indir dn: ds: f5.3;
		label dn_1="DN1" ds_1="DS1" dn_2="DN2" ds_2="DS2";
		run;
		%end;
title; title2;
%mend;

*%necsuff_ord(data=goeg.goeg, y=ord, ylist=2 1 0, x=altergr_, xclass=altergr_,
						odssel=none, print=1);
*%necsuff_ord(data=goeg.goeg, y=ord, ylist=2 1 0, x=alter65,
						odssel=none, print=1);
