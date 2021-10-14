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
* Statistics in Medicine 2019, 38:4733–4748
*
* Author:  Andreas Gleiss
* Version: 1.4 (bug fixed for refcat=last)
* Date:    5 Oct 2021
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
*
*******************************************************************;

%macro necsuff(data, y, x, freq=, refcat=first, inpred=, inpredvar=, 
				odssel=ResponseProfile, print=1);
%if &inpred= %then %do;
	data _work;
		set &data;
		if missing(&x) or missing(&y) then delete;
		%if &freq= %then freq=1;
				   %else freq=&freq;;
		run;
	ods select &odssel;
	proc logistic data=_work outest=_betas;
		class &y / param=ref ref=&refcat;
		model &y=&x;
		output out=_est1 pred=p_i_dach;
		freq freq;
		run;
	data _betas;
		set _betas(rename=(intercept=beta0 &x=beta1));
		_dummy=1;
		keep _dummy beta0 beta1;
		run;
	proc means data=_work noprint;
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
		run;
	data _betas;
		_dummy=1; beta0=.; beta1=.; output;
		run;
	proc means data=_est1 noprint;
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
	by descending p_i_dach;
	run;
data _est1;
	set _est1;
	_dummy=1;
	_row=_n_;
	run;
data _est1_;
	merge _est1 
		  _my1(keep=p_bar _dummy _freq_ rename=(_freq_=n))
		  _betas;
	by _dummy;
	retain dn1_sum ds1_sum dn2_sum ds2_sum n_smaller n_larger (0 0 0 0 0 0);

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
	by _dummy _row;
	*where _row=n;
	if not last._dummy then delete;
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
	keep progn_fact p_bar alpha or ev_indir dn_1 ds_1 dn_2 ds_2;
	run;
ods select all;
%if &print=1 %then %do;
	title "Degrees of necessity and sufficiency";
	title2 "Outcome = &y";
	proc print data=_necsuff_ noobs label;
		var progn_fact p_bar alpha or ev_indir dn_1 ds_1 dn_2 ds_2;
		format or p_bar alpha ev_indir dn: ds: f5.3;
		label 	progn_fact="Prognostic factor"
				or="OR" ev_indir="EV" p_bar="est.P(D)"
				dn_1="DN1" ds_1="DS1" dn_2="DN2" ds_2="DS2";
		run;
	title; title2;
	%end;
%mend;

*%necsuff(data=sh_tab3.hl_prost, y=capsule, x=psa);
*%necsuff(data=sh_tab3.hl_prost, y=capsule, x=gleason, odssel=none);
/*
data lungca; * Swedish cohort study (Nilsson et al, 2001);
	do i=1 to 8120;
		x=0; y=0; output;
		end;
	do i=1 to 36;
		x=0; y=1; output;
		end;
	do i=1 to 4331;
		x=1; y=0; output;
		end;
	do i=1 to 177;
		x=1; y=1; output;
		end;
	run;
%necsuff(data=lungca, y=y, x=x);
*/
