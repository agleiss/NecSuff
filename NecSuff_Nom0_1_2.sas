*******************************************************************;
*
* NECSUFF_NOM0
* ============
*
* SAS-macro for computation of degrees of necessity and sufficiency
* for a nominal outcome variable with reference category
*
* Gleiss, A., Henderson, R., Schemper, M., Degrees of necessity 
* and of sufficiency: further results and extensions, with an 
* application to covid-19 mortality in Austria
* accepted by Statistics in Medicine
*
* Author:  Andreas Gleiss
* Version: 1.1 (bug fixed)
* Date:    15 Oct 2020
*
* Macro parameters:
* =================
* 
* data		name of SAS data set
* y			name of ordinal outcome variable with values from ylist 
* ylist		outcome codes
* refcat	reference category (default = 0) 
* x			names of independent variables (or linear predictor)
* odssel	control output from proc logistic 
*			(default: ParameterEstimates OddsRatios)
* print 	=1 to print results (default), =0 (e.g., for simulations)
*
*******************************************************************;

%let mypath=E:\Andreas\Makros & Vorlagen\SAS-Macros; * change to path containing NecSuff.SAS;
%include "&mypath\NecSuff.SAS";

%macro necsuff_nom0(data=, y=, ylist=, refcat=0, x=, 
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
	run;
proc datasets noprint;
	delete _all;
	run;
data _works;
	set &data;
	_zeil=_n_;
	freq=1;
	format &y f3.0; * quasi unformat;
	run;
ods select &odssel;
proc logistic data=_works;
	model &y(descending ref="&refcat")=&x / link=glogit; 
	output out=_estpo pred=p_i_dach xbeta=xbeta;
	run;

%do iy=1 %to &ny;
	%if &&y&iy~=&refcat %then %do;
		data _po&iy;
			set _estpo;
			/*%if &refcat~= %then %do;
				where &y in (&refcat,&&y&iy);
				%end;*/
			where _level_=&&y&iy;
			if not missing(&y) then _yd=(&y = &&y&iy);
			run;
		data _po_in&iy;
			set _po&iy(where=(not missing(_yd) and not missing(xbeta)));
			run;
		proc means data=_estpo noprint nway;
			class _zeil;
			where _level_ in (&refcat,&&y&iy);
			var p_i_dach;
			output out=_sums sum=_sum;
			run;
		data _po_in&iy._;
			merge _po_in&iy
				 _sums;
			by _zeil;
			p_i_dach_=p_i_dach/_sum;
			run;
		%necsuff(data=_po&iy, y=_yd, x=xbeta, freq=, refcat=first, 
				 inpred=_po_in&iy._, inpredvar=p_i_dach_, 
				 odssel=none, print=0);

		data _erg;
			set _necsuff_;
			cat=&&y&iy;
			keep cat dn: ds: ev: or p_bar alpha;
			run;
		proc append data=_erg base=_all force;
			run;
		%end;
	%end;
proc means data=_tab(
		%if &refcat= %then %do;
			where=(not missing(&y)) 
			%end;
		%else %do;
			where=(not missing(&y) and &y~=&refcat)
			%end;
		) noprint;
	var count;
	output out=_tabsum sum=cum_freq;
	run;
data _tab_;
	set _tab(rename=(&y=cat) 
		%if &refcat= %then %do;
			where=(not missing(cat)) 
			%end;
		%else %do;
			where=(not missing(cat) and cat~=&refcat)
			%end;
		drop=cum:);
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
	title2 "Nominal outcome = &y";
	proc print data=_all2 noobs label;
		var cat p_bar alpha or ev_indir  
			dn_1 ds_1 dn_2 ds_2 weight;
		format or p_bar alpha ev_indir dn: ds: f5.3;
		label dn_1="DN1" ds_1="DS1" dn_2="DN2" ds_2="DS2";
		run;
		%end;
title; title2;
%mend;

*%necsuff_nom0(data=goeg.goeg, y=ord, ylist=0 1 2, refcat=0, x=alter65, 
				odssel=all, print=1);

