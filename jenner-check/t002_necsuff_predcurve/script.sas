*******************************************************************;
*
* NECSUFF
* =======
*
* SAS-macro for computation of degrees of necessity and sufficiency
* for a dichotomous outcome variable
*
* Gleiss, A. & Schemper, M., Quantifying degrees of necessity and 
* sufficiency in cause-effect relationships with dichotomous and 
* survival outcome
* Statistics in Medicine 2019, 38:47334748
*
* Author:  Andreas Gleiss
* Version: 1.5 (macro parameter by added)
* Date:    27 Aug 2025
*
* Macro parameters:
* =================
* 
* data		name of SAS data set
* y			name of dichotomous outcome variable 
* x			name of independent variable (or linear predictor)
* freq		name of optional frequency variable (default: none)
* refcat	reference category of y (first=default or last)
* inpred	dataset containing predictions
*			(data and x can be omitted)
* inpredvar	name of variable in inpred dataset containing 
*			predictions 
*			(caution: y and inpredvar must be compatible)
* odssel	control output from proc logistic 
*			(default: ResponseProfile)
* print 	=1 to print results (default), =0 (e.g., for simulations)
* by		optional variable name for by processing
*
*******************************************************************;

%macro necsuff(data, y, x, freq=, refcat=first, inpred=, inpredvar=, 
				odssel=ResponseProfile, print=1, by=);

%if &by= %then %let by = _all;;

%if &inpred= %then %do;
	data _work;
		set &data;
		if missing(&x) or missing(&y) then delete;
		%if &freq= %then freq=1;
				   %else freq=&freq;;
		_all=1;
		run;
	proc sort data=_work;
		by &by;
		run;
	ods select &odssel;
	proc logistic data=_work outest=_betas;
		class &y / param=ref ref=&refcat;
		model &y=&x;
		output out=_est1 pred=p_i_dach;
		freq freq;
		by &by;
		run;
	data _betas;
		set _betas(rename=(intercept=beta0 &x=beta1));
		_dummy=1;
		keep _dummy beta0 beta1;
		run;
	proc means data=_work noprint nway;
		class &by;
		var &y;
		output out=_my1 mean=p_bar max=_dummy;
		freq freq;
		run;
	data _my1;
		set _my1;
		%if &refcat=last %then p_bar=1-p_bar;;
		run;
	%end;
%else %do;
	data _est1;
		set &inpred(where=(not missing(&inpredvar) and not missing(&y)));
		p_i_dach=&inpredvar;
		freq=1;
		_all=1;
		run;
	data _betas;
		_dummy=1; beta0=.; beta1=.; output;
		run;
	proc means data=_est1 noprint nway;
		class &by;
		var &y; *&inpredvar;
		output out=_my1 mean=p_bar;
		run;
	data _my1;
		set _my1;
		%if &refcat=last %then p_bar=1-p_bar;;
		_dummy=1;
		run;
	%end;
proc sort data=_est1;
	by &by descending p_i_dach;
	run;
data _est1;
	set _est1;
	_dummy=1;
	_row=_n_;
	run;
data _est1b;
	merge _est1 
		  _betas;
	by _dummy;
	run;
data _est1_;
	merge _est1b 
		  _my1(keep=&by p_bar _dummy _freq_ rename=(_freq_=n));
	by &by;
	retain dn1_sum ds1_sum dn2_sum ds2_sum n_smaller n_larger;

	if first.&by then do;
		dn1_sum=0; ds1_sum=0; dn2_sum=0; ds2_sum=0; n_smaller=0; n_larger=0;
		end;

	smaller=(p_i_dach<p_bar);
	larger=(p_i_dach>p_bar);
	n_smaller=n_smaller+smaller*freq;
	n_larger=n_larger+larger*freq;

	if p_bar~=0 then do;
		if smaller then dn1_sum=dn1_sum + freq*((p_bar - p_i_dach)/p_bar)**2;
		if smaller then dn2_sum=dn2_sum + freq*((p_bar - p_i_dach)/p_bar);
		end;
	if p_bar~=1 then do;
		if larger then ds1_sum=ds1_sum + freq*((p_i_dach - p_bar)/(1-p_bar))**2;
		if larger then ds2_sum=ds2_sum + freq*((p_i_dach - p_bar)/(1-p_bar));
		end;
	run; 

data _necsuff_;
	set _est1_;
	by &by _row;
	if not last.&by then delete;
	if n_smaller~=0 then do;
		dn_1_=dn1_sum/n_smaller; 
		dn_2=dn2_sum/n_smaller; 
		end;
	else do;
		dn_1_=0;
		dn_2=0;
		end;
	if n_larger~=0 then do;
		ds_1_=ds1_sum/n_larger;
		ds_2=ds2_sum/n_larger;
		end;
	else do;
		ds_1_=0;
		ds_2=0;
		end;
	weight_dn=p_bar/(1-p_bar) * n_smaller/n;
	weight_ds=(1-p_bar)/p_bar * n_larger/n;
	ev_indir=weight_dn*dn_1_ + weight_ds*ds_1_;
	dn_1=sqrt(dn_1_);
	ds_1=sqrt(ds_1_);
	alpha=n_larger/n;
	or=exp(beta1);
	progn_fact="&x";
	keep &by progn_fact p_bar alpha or ev_indir dn_1 ds_1 dn_2 ds_2;
	run;
ods select all;
%if &print=1 %then %do;
	title "Degrees of necessity and sufficiency";
	title2 "Outcome = &y";
	proc print data=_necsuff_ noobs label;
		var &by progn_fact p_bar alpha or ev_indir dn_1 ds_1 dn_2 ds_2;
		format or p_bar alpha ev_indir dn: ds: f5.3;
		label 	progn_fact="Prognostic factor"
				or="OR" ev_indir="EV" p_bar="est.P(D)"
				dn_1="DN1" ds_1="DS1" dn_2="DN2" ds_2="DS2";
		run;
	title; title2;
	%end;
%mend;

*******************************************************************;
*
* NECSUFF_PREDCURVE
* =================
*
* SAS-macro for plotting the predictiveness curve and showing
* areas that correspond to DN and DS;
* To be applied directly after call of %necsuff();
* Uses datasets _est1, _est1_ and _necsuff_ produced by %necsuff()
*
* Gleiss, A. Visualizing a markers degrees of necessity and of 
* sufficiency in the predictiveness curve, 
* BMC Medical Research Methodology (2025) 25:107.
*
* Author:  Andreas Gleiss
* Version: 1.0
* Date:    12 Sept 2024
*
* Macro parameters:
* =================
* 
* show_areas		select type of plot:
*					0 = show predictiveness curve only
*					1 = show areas A_N and A_S
*					2 = show numerators and denominators of DN, DS
*						(default)
* where				optional restriction (e.g. to (non-)events only)
* y					name of outcome variable in original dataset
*					(default = event)
* event				label for event (e.g. death) (default = event)
* linetype			step (for categorical predictor) or
*					series (for continuous predictor without areas)
* addhline			value in [0,1] for adding horizontal reference 
*					line
* xadd				additional description for x axis label
* showmeas			1 or 0 to show or not show total gain (TG)
*					and its standardizations STG and STG2
*
*******************************************************************;

%macro necsuff_predcurve(show_areas = 2, where = %str(1 = 1),
						 y = event, event = event, 
						 linetype = step, addhline =, xadd =,
					 	 showmeas = 1);

data _data;
	set _est1; /* produced by %necsuff() */
	xbeta = exp(p_i_dach) / (1 + exp(p_i_dach));
	where &where;
	run;

proc rank data = _data
		  out = _rest percent ties = low;
	var xbeta;
	ranks r_xbeta;
	run;

proc means data = _data noprint;
	var p_i_dach;
	output out = _p_bar mean = p_bar;
	run;
data _null_;
	set _p_bar;
	call symput('p_bar', p_bar);
	run;

data _rest_;
	set _rest;
	p_bar = &p_bar;
	q_xbeta = r_xbeta / 100; * quantile;
	low = (p_i_dach <= p_bar);
	high = (p_i_dach > p_bar);
	if low then q_xbeta_low = q_xbeta;
	if high then q_xbeta_high = q_xbeta;
	run;

proc freq data = _rest_ noprint;
	tables low / out = _ltab;
	run;
data _null_;
	set _ltab(where = (low = 1));
	q0 = percent / 100;
	call symput('q0', q0);
	run;
data _rest_;
	set _rest_;
	q0 = &q0;
	run;

proc sort data = _rest_ out = _rest_sorted;
	by xbeta;
	run;

data _rest_sorted;
	set _rest_sorted;
	_dummy = 1;
	_row = _n_;
	run;

* add corner points for propor plotting;
data _add1;
	set _rest_sorted(where = (missing(q_xbeta_high)));
	by _dummy;
	q_xbeta_low = 0;
	q_xbeta = 0;
	if not first._dummy then delete;
	run;
data _add2;
	set _rest_sorted(where = (missing(q_xbeta_high)));
	by _dummy;
	q_xbeta_low = q0;
	q_xbeta = q0;
	if not last._dummy then delete;
	run;
data _add3;
	set _rest_sorted(where = (missing(q_xbeta_low)));
	by _dummy;
	q_xbeta_high = q0;
	q_xbeta = q0;
	if not first._dummy then delete;
	run;
data _add4;
	set _rest_sorted(where = (missing(q_xbeta_low)));
	by _dummy;
	q_xbeta_high = 1;
	q_xbeta = 1;
	if not last._dummy then delete;
	run;

data _plot;
	set _rest_sorted
		_add1 _add2 _add3 _add4;
	run;

proc sort data = _plot;
	by xbeta;
	run;

ods graphics on / height = 450px width = 450px border = off;

proc sgplot data = _plot noautolegend;
	xaxis min = 0 max = 1 label = "Risk quantile q &xadd"
		  offsetmin = 0 offsetmax = 0;
	yaxis min = 0 max = 1 label = "Risk of &event"
		  offsetmin = 0 offsetmax = 0;

	%if &addhline ~= %then %do;
		refline &addhline / lineattrs = (pattern = 3);
		%end;

	%if &show_areas = 2 %then %do;
		refline &q0 / axis = x label = "q(*ESC*){unicode '2080'x}" labelattrs = GraphUnicodeText;
		band x = q_xbeta_low lower = 0 upper = p_bar / fillattrs = (color = grey);* transparency = 0.75);
		band x = q_xbeta_high lower = p_bar upper = 1 / fillattrs = (color = lightgrey);* transparency = 0.75);
		%end;
	%if &show_areas > 0 %then %do;
		band x = q_xbeta_low lower = p_i_dach upper = p_bar / fillattrs = (color = grey) /* transparency = 0.75)*/
															  fillpattern fillpatternattrs = (color = black pattern = L1)
																 	type = step;
		band x = q_xbeta_high lower = p_bar upper = p_i_dach / fillattrs = (color = lightgrey) /* transparency = 0.75)*/
															   fillpattern fillpatternattrs = (color = black pattern = L1)
																 	type = step;
		%end;

	refline &p_bar / lineattrs = (color = black) label = "P(D)";

	&linetype x = q_xbeta y = p_i_dach / lineattrs = (color = black);

	run;

%if &showmeas %then %do;

	data _allout;
		set _necsuff_; * produced in %necsuff();
		A_N = DN_2 * p_bar * (1 - alpha); * alpha = 1 - p_0;
		A_S = DS_2 * (1 - p_bar) * alpha;
		TG = A_N + A_S; * total gain;
		STG = TG / (2 * p_bar * (1 - p_bar));
		STG2 = (A_N + A_S) / (p_bar * (1 - alpha) + (1 - p_bar) * alpha);
		run;

	proc print data = _allout noobs label;
		var A_N A_S TG STG STG2; 
		run;

	%end;
%mend;



/* A fixed continuous-marker cohort: marker x and dichotomous outcome y, so the
   predictiveness curve has spread. %necsuff fits the logistic model and produces
   the datasets _est1/_est1_/_necsuff_; %necsuff_predcurve then draws the
   predictiveness curve shading the areas that correspond to DN and DS. */
data cohort;
	input x y;
	datalines;
-2.0 0
-1.7 0
-1.5 0
-1.2 0
-1.0 0
-0.8 0
-0.6 0
-0.4 0
-0.2 0
0.0 0
0.2 1
0.4 0
0.5 1
0.7 1
0.9 0
1.0 1
1.2 1
1.4 1
1.5 1
1.7 1
1.9 1
2.1 1
;
run;
%necsuff(data=cohort, y=y, x=x);
%necsuff_predcurve(show_areas=2, y=y, event=event, linetype=series, showmeas=0);
