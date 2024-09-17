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
* Gleiss, A. Visualizing a marker’s degrees of necessity and of 
* sufficiency in the predictiveness curve, submitted (2024)
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
