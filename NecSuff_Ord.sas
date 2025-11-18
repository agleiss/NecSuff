*******************************************************************;
*
* NECSUFF_ORD
* ===========
*
* SAS-macro for computation of degrees of necessity and sufficiency
* for an ordinal outcome variable (assuming proportional odds or
* based on input predictions)
*
* Gleiss, A., Henderson, R., Schemper, M., Degrees of necessity 
* and of sufficiency: further results and extensions, with an 
* application to covid-19 mortality in Austria
* Statistics in Medicine 2021, 40:3352–3366;
*
* Author:  Andreas Gleiss
* Version: 1.2 (direct input of predictions via inpredo dataset,
*				modopt and by arguments)
* Date:    16 Oct 2025
*
* Macro parameters:
* =================
* 
* data		name of SAS data set
* y			name of ordinal outcome variable with values from ylist 
* ylist		outcome codes sorted from worst to best outcome
*			(lowest code corresponds to best outcome)
* x			names of independent variables (or linear predictor)
* xclass	list of independent class variable names
* modopt	optional options for model statement in proc logistic
* inpredo	dataset containing predictions
*			(x can be omitted, needs _level_ variable)
* inpredvaro	name of variable in inpredo dataset containing 
*			predictions (cumulative predicted probabilities per _level_)
*			(caution: y and inpredvar must be compatible)
* odssel	control output from proc logistic 
*			(default: ParameterEstimates OddsRatios)
* print 	=1 to print results (default), =0 (e.g., for simulations)
* by		optional variable name for by processing
*
*******************************************************************;

%let mypath=E:\Andreas\Makros & Vorlagen\SAS-Macros; * change to path containing NecSuff.SAS;
%include "&mypath\NecSuff.SAS";

%macro necsuff_ord(data=, y=, ylist=2 1 0, x=, xclass=, modopt=%str(link = logit),
					inpredo=, inpredovar=,
					odssel=ResponseProfile ParameterEstimates OddsRatios, print=1,
					by=);

%let ny=0;
%do %while(%scan(&ylist,&ny+1)~=);
	%let ny=%eval(&ny+1);
	%let y&ny=%scan(&ylist,&ny);
	%end;

proc format;
	value cat 999="W.sum";
	run;
%if &by= %then %do;
	proc freq data=&data noprint;
		tables &y / out=_tab;* outcum;
		run;
	data _tab;
		set _tab;
		_all =1;
		run;
	%end;
%else %do;
	proc freq data=&data noprint;
		tables &by * &y / out=_tab;*(rename = (pct_row=cumpct)) outpct;
		run;
	%end;

%if &by= %then %let by = _all;;

proc datasets noprint;
	delete _all;
	run;
data _works;
	set &data;
	_zeil=_n_;
	_all=1;
	format &y f3.0; * quasi unformat;
	run;
proc sort data=_works;
	by &by;
	run;

%if &inpredo= %then %do;
	ods select &odssel;
	proc logistic data=_works outest=_parest;
		class &y &xclass / order=internal;
		model &y(descending)=&x / &modopt; 
		output out=_estpo pred=p_i_dach /* cumulative */ xbeta=xbeta;
		by &by;
		run;
	%end;
%else %do;
	data _estpo;
		set &inpredo;
		rename &inpredovar = p_i_dach; * cumulative;
		_all=1;
		run;
	proc sort data=_estpo;
		by &by;
		run;
	%end;

%do iy=1 %to %eval(&ny-1);
	data _po&iy;
		set _estpo;
		where _level_=&&y&iy;
		if not missing(&y) then _yd=(&y >= &&y&iy); * since lowest code = best outcome;
		run;
	data _po_in&iy;
		set _po&iy(where=(not missing(_yd) /*and not missing(xbeta)*/));
		run;
	%necsuff(data=_po&iy, y=_yd, /*x=xbeta,*/ freq=, refcat=first, 
			 inpred=_po_in&iy, inpredvar=p_i_dach, 
			 odssel=none, print=0, by=&by);

	data _erg;
		set _necsuff_;
		cat=&&y&iy;
		keep &by cat dn: ds: ev: or p_bar alpha;
		run;
	proc append data=_erg base=_all force;
		run;
	%end;
proc means data=_tab(where=(&y>&&y&ny)) noprint nway;
	class &by;
	var count;
	output out=_tabsum sum=cum_freq;
	run;
data _tab_;
	merge _tab(rename=(&y=cat) where=(cat>&&y&ny)) /*drop=cum:)*/
		  _tabsum(keep=cum_freq &by); 
	by &by;
	run; 
proc sort data=_all;
	by &by cat;
	run;
data _all_;
	merge _all
		  _tab_;
	by &by cat;
	weight=count/cum_freq;
	dn1_term=dn_1*weight;
	ds1_term=ds_1*weight;
	dn2_term=dn_2*weight;
	ds2_term=ds_2*weight;
	ev_term=ev_indir*weight;
	run;
proc means data=_all_ sum noprint nway;
	class &by;
	var ev_term dn1_term ds1_term dn2_term ds2_term;
	output out=_wsums sum=;
	run; * auch EV als (gewichtetes) Mittel aus Einzel-EVs;
data _all2;
	set _all_
		_wsums(in=wsum rename=(dn1_term=dn_1 ds1_term=ds_1 
							   dn2_term=dn_2 ds2_term=ds_2 ev_term=ev_indir));
	if wsum then cat=999;
	format cat cat.;
	run;
proc sort data = _all2;
	by &by cat;
	run;
%if &print=1 %then %do;
	title "Degrees of necessity and sufficiency";
	title2 "Ordinal outcome = &y";
	proc print data=_all2 noobs label;
		var &by cat p_bar alpha or ev_indir  
			dn_1 ds_1 dn_2 ds_2 weight;
		format or p_bar alpha ev_indir dn: ds: f5.3;
		label dn_1="DN1" ds_1="DS1" dn_2="DN2" ds_2="DS2";
		run;
		%end;
title; title2;
%mend;

*%necsuff_ord(data=goeg.goeg, y=ord, ylist=2 1 0, x=altergr_, xclass=altergr_,
						odssel=none, print=1);
