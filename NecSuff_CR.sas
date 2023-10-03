*******************************************************************;
*
* NECSUFF_CR
* ==========
*
* SAS-macro for computation of degrees of necessity and sufficiency
* (DN, DS) and explained variation (EV, direct estimate) for 
* competing risks survival data based on Fine & Gray model
*
* This macro has been adapted from the macro %surev by M.Schemper 
* and G. Heinze 
*
* Output: marginal EV^CR, DS^CR and DS^CR
*
* Gleiss, A., Gnant, M., Schemper, M., Explained variation and 
* degrees of necessity and of sufficiency for competing risks 
* survival data. Submitted to Biometrical Journal.
*
*
* Author:  Andreas Gleiss
* Version: 1.2 (removed q and corrected p_ij, 
*				added indirect EV estimate)
* Date:    03 Oct 2023
*
* Macro parameters:
* =================
* 
* data		name of SAS data set
* time		name of time variable of survival outcome 
* status	name of status variable of survival outcome 
* resp		value of status variable indicating event of
*			interest (default=1)
* cr		value of status variable indicating comepting
*			event (default=2)
* varlist	list of independent variables 
*			(leave empty if only strata are used)
* strat		name of stratum variable 
*			(consecutive numbers starting from 1)
* out		optional name of output dataset
* useranks	if =1 replaces survival times by their ranks (default)
* print		if =1 prints results (default)
*
*******************************************************************;


%MACRO NecSuff_CR(DATA=, TIME=TIME, STATUS=STATUS, RESP=1, cr=2,
					varlist=, strat=, 
					out=_out, useranks=1,print=1);


%if &strat =  %then %let strat = _all;
data _data;
	set &data;
	_all = 1; * dummy stratum for all observations if no strata are given;
	run;

%if &varlist =  %then %let varlist = _all;

* transform varlist into single accessible macro variables;
%let nvar = 0;
%do %while(%scan(&varlist,&nvar + 1) ~= );
	%let nvar = %eval(&nvar + 1);
	%let v&nvar = %scan(&varlist,&nvar);
	%end;

%let vz = &nvar; * number of prognostic factors;

proc means data = _data noprint;
	var &strat;
	output out = _maxstrat max = maxstrat;
	run;
data _null_;
	set _maxstrat;
	call symput('maxstrat',maxstrat);
	run;
%let maxstrat = %eval(&maxstrat); * number of strata as macro variable;

DATA RX; 
	SET _DATA; 
	Z = &RESP; * status code for response (event of interest);
	z_ = &cr;  * status code for competing event;

	* set new status variable stat to 1 for event of interest and 2 for competing event;
	IF &STATUS = Z THEN STAT = 1;  
				   ELSE if &status = z_ then stat = 2; 
									    else STAT = 0; * censored;

	IF &TIME = . or STAT = . THEN DELETE;

	IF &TIME<0 THEN DO;
	   ERROR 'NEGATIVE TIME ENCOUNTERED, PROGRAM STOPS';
	   ABORT; END;

	* delete observations with missing value for prognostic factor or stratum;
	if missing(&strat) then delete;

	%DO J = 1 %TO &VZ; 
		IF missing(&&V&J) THEN DELETE;   
		%END;    
		NN + 1;
	run;

%if &useranks = 1 %then %do; * rank transformation to survival times;
	proc rank data = rx out = rx;
		var &time;
		ranks &time._r;
		run;
	%let timebu = &time;
	%let time = &time._r;
	%end;

* SUBTRACT MEANS FROM COVARIATES AND SORT BY ASCENDING TIME;
PROC STANDARD DATA = RX OUT = RR MEAN = 0;
	VAR %DO J = 1 %TO &VZ; 
		&&V&J 
		%END;  
		;
	run;

PROC SORT DATA = RR;
	BY &TIME;
	run;

data rr;
	set rr;
	_line_ = _n_;
	run;

data _covs; * dataset for baseline statement of proc phreg below;
	%DO J = 1 %TO &VZ; 
		&&V&J = 0; * corresponds to mean due to standardization above; 
		%END; 
	run;

* stratified Fine & Gray model for event of interest (EoI);
PROC PHREG DATA = RR OUTEST = V noprint; * store coefficients;
	MODEL &TIME*STAT(0) =  %DO J = 1 %TO &VZ; &&V&J %END;  
		/ eventcode = 1; 
	BASELINE covariates = _covs OUT = S cif = cif0;
	strata &strat;
run;
	* restructure dataset of baseline CIFs for EoI to have strata beside each other;
	proc sort data = s out = _s;
		by &time;
		run;
	data ss_;
		merge
		%do is = 1 %to &maxstrat; 
			_s(keep = &time cif0 &strat where = (&strat = &is) rename = (cif0 = cif0&is)) 
			%end;;
		by &time;
		drop &strat ;
		run;

	data ss_; * prepare for strata without EoI;
		set ss_; 
		%do is = 1 %to &maxstrat;
			if &time = 0 and cif0&is = . then cif0&is = 0;
			%end;
		run;

	data s_; * fill missing CIF values with last non-missing value;
		set ss_;
		retain cif0lag1-cif0lag&maxstrat;
		%do is=1 %to &maxstrat;
			if cif0&is=. then cif0&is=cif0lag&is;
			%end;
		output;
		%do is=1 %to &maxstrat;
			cif0lag&is=cif0&is;
			%end;
		run;


* stratified Fine & Gray model for competing event (CE);
* (needed for estimation of p(.|x) in estimation of M(.|x));
PROC PHREG DATA = RR OUTEST = V_CE noprint; * store coefficients;
	MODEL &TIME*STAT(0) =  %DO J = 1 %TO &VZ; &&V&J %END;  
	/ eventcode = 2;
	BASELINE covariates = _covs OUT = S_CE cif = cif0;
	strata &strat;
run;
	* restructure dataset of baseline CIFs for CE to have strata beside each other;
	proc sort data = s_CE out = _s_CE;
		by &time;
		run;
	data ss_CE;
		merge
		%do is = 1 %to &maxstrat; 
			_s_CE(keep = &time cif0 &strat where = (&strat = &is) rename = (cif0 = cif0&is)) 
			%end;;
		by &time;
		drop &strat ;
		run;

	data ss_CE; * prepare for strata without CE;
		set ss_CE; 
		%do is = 1 %to &maxstrat;
			if &time = 0 and cif0&is = . then cif0&is = 0;
			%end;
		run;
	data s_CE_; * fill missing CIF values with last non-missing value;
		set ss_CE;
		retain cif0lag1-cif0lag&maxstrat;
		%do is=1 %to &maxstrat;
			if cif0&is=. then cif0&is=cif0lag&is;
			%end;
		output;
		%do is=1 %to &maxstrat;
			cif0lag&is=cif0&is;
			%end;
		run;

ods select none;

* unconditional CIF for event of interest;
PROC LIFETEST DATA = RR outcif = t noprint;
	TIME  &TIME*STAT(0) / eventcode = 1; 
	run;
data t;
	set t;
	survival = 1 - cif; 
	run;

* number of different status values (1: no competing events, 2: competing events);
proc means data = rr noprint;
	var stat;
	output out = stats max = max_stat;
	run;
data _null_;
	set stats;
	call symput('max_stat',max_stat);
	run;

%if &max_stat = 2 %then %do; * competing events present;
	* unconditional CIF for CE (needed for estimation of p within M(.));
	PROC LIFETEST DATA = RR outcif = t_CE noprint;
	 TIME  &TIME*STAT(0) / eventcode = 2;
	run;
	data t_ce;
		set t_ce;
		survival = 1-cif; 
		run;
	%end;

* reverse KM estimate G(t);
PROC LIFETEST DATA = RR OUTSURV = B NOPRINT ;  
	TIME  &TIME*STAT(1,2); * censor both event types;
	run;

ods select all;

* prepare data for proc iml;
DATA U3(KEEP =  &TIME  KM cif) ;
 SET T;
 KM = SURVIVAL;
 IF KM ne .;
run;

DATA U2(KEEP =  &TIME  KM cif) ;
 SET U3;
 BY &TIME;
 IF LAST.&TIME;
run;

DATA U1(KEEP =  &TIME BKM) ;
 SET B; 
 if survival ne . then BKM = SURVIVAL;
 RETAIN BKM;
run;

DATA U;
 MERGE U1 U2(IN = W);
 BY &TIME;
 IF W;
run;

DATA u5;
 MERGE S_(IN = W) U ;
 BY  &TIME;
 IF W;
run;

DATA u6;
 SET U5 ;
 BY  &TIME;
 IF FIRST.&TIME;
 IF &TIME > 0.0; 
run;

data u6;
 set u6;
 _zeitp_ = &time;
 %do is = 1 %to &maxstrat;
	_cif0&is._ = cif0&is; * baseline CIF;
	%end;
 _cif_ = cif;	* uncond. CIF;
 _bkm_ = bkm;   * reverse KM;
 _line_ = _n_;
 keep _line_ _zeitp_ _cif0: _cif_ _bkm_;
run;

data u7;
	merge rr u6;
	by _line_; 
	run;


* analogue for CE;

%if &max_stat = 2 %then %do;
	DATA U3_(KEEP =  &TIME cif ) ;
	 SET T_CE;
	 cif = cif;
	 IF cif ne .;
	run;

	DATA U2_(KEEP =  &TIME cif) ;
	 SET U3_;
	 BY &TIME;
	 IF LAST.&TIME;
	run;

	DATA U4_;
	 MERGE S_CE_(IN = W) U2_;
	 BY  &TIME;
	 IF W;
	run;

	DATA u6_;
	 SET U4_;
	 BY  &TIME;
	 IF FIRST.&TIME;
	 IF &TIME > 0.0; 
	run;

	data u6_;
	 set u6_;
	 _zeitp_ = &time;
	  %do is = 1 %to &maxstrat;
		_cif0&is._ = cif0&is; * baseline CIF;
		%end;
	 _cif_ = cif;	  * uncond. CIF;
	 _line_ = _n_;
	 keep _line_ _zeitp_ _cif0: _cif_;
	run;
	%end;
%else %do;
	* u6_ empty;
	%end;


data _time_; * store starting time;
	x = time();
	output;
	run;

		PROC IML;

				start explvar(PAR, par_, DATALINES, datalines_, IOPARMS);
					*** code translated from original fortran SUREV.FOR by G. H. 2012 Apr 12;
					* par			regression coefficients from Fine & Gray model for EoI;
					* par_			regression coefficients from Fine & Gray model for CE;
					* datalines		big table containing various data (see below);
					* datalines_	table containing data for CE (see below);
					* ioparms		vector of input parameters (see below);

					* prepare vectors and matrices according to dimensions given in ioparms;
					zeit = repeat(0, ioparms[1]);
					exbeta = repeat(0, ioparms[1]);
					exbeta_ = repeat(0, ioparms[1]);
					istat = repeat(0, ioparms[1]);
					inde = repeat(0, ioparms[1]);
					inde_ = repeat(0, ioparms[1]);
					zeitp = repeat(0, ioparms[3]);
					cif0 = j(ioparms[3], ioparms[4], 0);	* baseline CIF;
					cif = repeat(0, ioparms[3]);			* uncond. CIF;

					if ioparms[5]>0 then do;
						zeitp_ = repeat(0, ioparms[5]);		* distinct failure times CE;
						cif0_ = j(ioparms[5], ioparms[4], 0);	* baseline CIF for CE;
						cif_ = repeat(0, ioparms[5]);			* uncond. CIF for CE (only at distinct CE failure times);
						end;
					bkm = repeat(0, ioparms[3]);			* reverse KM;
					dup = repeat(0, ioparms[3]);
				    cov = j(ioparms[1], ioparms[2], 0);

				   n = ioparms[1];   * # individuals;
				   nk = ioparms[2];  * # covariables;
				   np = ioparms[3];  * # distinct failure times for event of interest;
				   ns = ioparms[4];  * # strata;
   				   np_ = ioparms[5]; * # distinct failure times for competing event;
				    												/* for each individual: */
				   zeit = datalines[, 1]; 							/* time variable */
				   istat = datalines[, 2];							/* status variable */
				   cov = datalines[, 3:(nk + 2)]; 					/* covariate values */
				   strt = datalines[, nk + 3]; 						/* stratum number */
				   zeitp = datalines[, nk + 4]; 					/* distinct EoI times */
				   cif0 = datalines[, (nk + 5):(nk + 5 + ns-1)];	/* baseline CIF for EoI */
				   cif = datalines[, (nk + 5 + ns)];				/* uncond. CIF for EoI */
				   bkm = datalines[, (nk + 5 + ns + 1)];			/* reverse KM estimate */

				   if np_>0 then do;
					    zeitp_ = datalines_[, 1]; 					/* distinct CE times */
				   		cif0_ = datalines_[, 2:(ns + 1)];			/* baseline CIF for CE */
   				   		cif_ = datalines_[, (ns + 2)];				/* unconditional CIF for CE */
						end;

				   * end of all input;

				   	  * construction of pointer to distinct EoI failure times (inde) 
					  * and of counter of tied failure times of interest (dup = d_j);
				      IZ = 0;                                                             
				      do j  =  1 to np; * loop over distinct times to event of interest t(j);
				      	  dup[j] = 0;                                                      
						  mark31:
					      if iz  =  n then goto mark32;                                            
						  iz = iz + 1;                                                           
					      if istat[iz] = 0 then do;                                        
						      inde[iz] = j;                                                        
						      goto mark31;                                                           
					      end;                                                              
					      if istat[iz] = 1 & zeit[iz] = zeitp[j] then do; 
						      dup[j] = dup[j] + 1;                                                   
						      inde[iz] = j;                                                        
						      goto mark31;                                                           
					      end;                                                       
					      if istat[iz] = 2 then do; 
						      inde[iz] = j;                                                        
						      goto mark31;                                                           
					      end;                                                              
					      if zeit[iz] > zeitp[j] then iz = iz-1; 
				          mark30:
				      end;  
					  mark32:
					  ;
				   	  * construction of pointer (inde_) for CE;
				      IZ = 0;                                                             
				      do j  =  1 to np_; * loop over distinct times to competing event t(j)_;
						  mark310:
					      if iz  =  n then goto mark320;                                            
						  iz = iz + 1;                                                           
					      if istat[iz] = 0 then do;                                        
						      inde_[iz] = j;                                                        
						      goto mark310;                                                           
					      end;                                                              
					      if istat[iz] = 2 & zeit[iz] = zeitp_[j] then do; 
						      inde_[iz] = j;                                                        
						      goto mark310;                                                           
					      end;                                                       
					      if istat[iz] = 1 then do;
						      inde_[iz] = j;                                                     
						      goto mark310;                                                           
					      end;                                                              
					      if zeit[iz] > zeitp_[j] then iz = iz-1; 
				          mark300:
				      end;  
					  mark320:
					  ;

					  d1u_s = 0; * for dir. est. of V;                                                      
				      d1c_s = 0;                                                       
				      wei_s = 0;                                                        
					  dn = 0; * for DN;
		  			  ds = 0; * for DS;
					  d1u = 0; * für indir. est. of V;                                                   
				      d1c = 0;                                                       

					  * prepare prognostic index (xbeta for EoI, xbeta_ for CE)
					  *	from covariate values cov and regression coefficients par and par_;
					  do i = 1 to n;                                                      
					      xbeta = 0; /* EoI */                                                    
					      do i1 = 1 to nk; 
							 mark19:
							 xbeta = xbeta + cov[i, I1] * par[i1];                                    
					      end;
					   	  mark18:
						  exbeta[i] = exp(xbeta);

						  if np_ > 0 then do; * CEs present;
					      	xbeta_ = 0; /* CE */                                                     
					      	do i1 = 1 to nk; 
							 	mark190:
							 	xbeta_ = xbeta_ + cov[i, I1] * par_[i1];                                    
					      	end;
					   	  	mark180:
						  	exbeta_[i] = exp(xbeta_);
						  end;
				      end; 
						                            
				      do j = 1 to np; * BIG LOOP OVER DISTINCT EoI FAILURE TIMES; 
					      resc = 0;		* for summing conditional terms for direct estimate;
					      resk = 0;     * for summing unconditional terms for direct estimate;
			  			  term_n = 0;	* for summing DN terms; 
			  			  term_s = 0; 	* for summing DS terms;
						  n_kl = 0; 	* counter of individuals in protective X set;
						  n_gr = 0;		* counter of individuals in harmful X set;
					      resc_ind = 0; * for summing conditional terms for indirect estimate;
					      resk_ind = 0; * for summing unconditional terms for indirect estimate;
 
						  * estimate M(.) and M(.|x);
					      do i = 1 to n; * loop over individuals; 
 
							* direct estimate of V;
						      fg = 1 - (1 - cif0[j, strt[i]]) ** exbeta[i]; * Fine & Gray CIF estimate C(t_(j));

						      if zeit[i] > zeitp[j] then do; * alive at t(j);                                   
							      resc = resc + fg;                                                 
							      resk = resk + cif[j];
							      goto mark12;                                                          
						      end;                                                             
						      if zeit[i]  <=  zeitp[j] & istat[i]  = 1 then do; * event of interest at t_(j);
							      resc = resc + (1 - fg);                                                  
							      resk = resk + (1 - cif[j]);      
							      goto mark12;                                                           
						      end;                                                             
						      if zeit[i]  <=  zeitp[j] & istat[i]  =  0 then do; * censored at t_(j);             
								  ife = inde[i]; 

							      fgc = 1 - (1 - cif0[ife, strt[i]]) ** exbeta[i]; * Fine & Gray CIF estimate C(t_i);

								  if np_>0 then do; * CEs present;
								      ife_ = inde_[i]; 
									  cif_ife_ = cif_[ife_];

								      fgc_ = 1 - (1 - cif0_[ife_, strt[i]]) ** exbeta_[i]; * Fine & Gray CIF estimate C_CE(t_i);

									  end;
								  else do; * no CEs present;
								  	  cif_ife_ = 0;
								      fgc_ = 0;
									  end;

							      if fgc + fgc_ < 0.999999 
									then RESC = RESC + (fg*(1 - fg - fgC_)/(1 - fgC - fgC_) 
														+  (1 - fg)*(1 - (1 - fg - fgC_)/(1 - fgC - fgC_))); 

							      if cif[ife] + cif_ife_ < 0.999999 
									then  RESK = RESK + (cif[J]*(1 - cif[J] - cif_ife_)/(1 - cif[IFE] - cif_ife_) 
													  +  (1 - cif[J])*(1 - (1 - cif[J] - cif_ife_)/(1 - cif[IFE] - cif_ife_))); 
							  end;                     
						      if zeit[i]  <=  zeitp[j] & istat[i]  =  2 then do; * competing event at t_(j);
							      ife = inde[i];                                                  
							      RESC = RESC  + fg; 
							      RESK = RESK  + cif[J]; 
						      end;  
					          mark12:

							* indirect estimate of V;
							  resc_ind = resc_ind + 2 * fg * (1 - fg);
							  resk_ind = resk_ind + 2 * cif[j] * (1 - cif[j]); * independent of i  = > n times;

							* estimate of DN,  DS;
							  if fg < cif[j] then do; /* >>> NECESSARY <<< */
							  	n_kl = n_kl + 1;
								if cif[j]<1 then 
									term_n = term_n  +  ((cif[j]-fg)/cif[j])**2;
								end;
							  if fg > cif[j] then do; /* >>> SUFFICIENT <<< */
							  	n_gr = n_gr + 1;
								if cif[j]>0 then 
									term_s = term_s  +  ((fg-cif[j])/(1-cif[j]))**2;
								end;

				          end; /* i-loop */   

		      			  D1C = D1C + (RESC_ind / n) * DUP[J] / BKM[J]; * numerator for indirect V;
						  D1U = D1U + (RESK_ind / n) * DUP[J] / BKM[J]; * denom. for indirect V; 
	
						  D1C_S = D1C_S + (RESC / n) * DUP[J] / BKM[J]; * numerator for dir. V;
					      D1U_S = D1U_S + (RESK / n) * DUP[J] / BKM[J]; * denom. for dir. V;
					      WEI_S = WEI_S + DUP[J] / BKM[J]; 				* weight;                                                                             

						  if n_kl > 0 then 
					      	dn = dn + sqrt(term_n / n_kl) * DUP[J] / BKM[J];                                
						  if n_gr>0 then
					      	ds = ds + sqrt(term_s / n_gr) * DUP[J] / BKM[J];                                

					      mark8:
				      end; /* j-loop */                                                      

					  ioparms[3] = 1 - D1C / D1U;   * indirect V (weight sum cancels); 

				      D = D1U_S / WEI_S;                                          
				      Dx = D1C_S / WEI_S;   
				      ioparms[4] = (D - Dx) / D; /* direct estimate of V */

				      ioparms[5] = dn / wei_s;                                          
				      ioparms[6] = ds / wei_s;

					  * set small negative values to zero;
					  if ioparms[3] < 0 then ioparms[3] = 0; 
					  if ioparms[4] < 0 then ioparms[4] = 0; 
					  if ioparms[5] < 0 then ioparms[5] = 0; 
					  if ioparms[6] < 0 then ioparms[6] = 0; 

					  mark546:

				finish; /* explvar */


		* prepare input data to "explvar";
		use u7;
		 read all var("&time" || "stat"  %do j = 1 %to &vz;  || "&&v&j" %end; || "&strat"	|| "_zeitp_" 
						%do is = 1 %to &maxstrat; || "_cif0&is._" %end;
						|| "_cif_"
						||"_bkm_") 
		  into datlin;
		 close u7;

		use u6;
		 read all var("_zeitp_") into u6;
		 close u6;

		use v;
		 read all var("&v1" %if &vz>1 %then %do; 
		   %do j = 2 %to &vz;  || "&&v&j" %end;
		  %end;) into par;
		 close v;

		ioparms = repeat(0,8,1);
		ioparms[1] = nrow(datlin);	* # individuals;
		ioparms[2] = &vz;			* # covariates;
		ioparms[3] = nrow(u6);		* # distinct failure times for event of interest;
 		ioparms[4] = &maxstrat;		* # strata;

		if ioparms[3]<ioparms[1] then do; * if necessary fill empty places with zeroes;
		  do i = ioparms[3] + 1 to ioparms[1];
		   do j = &vz + 4 to &vz + 4 + &maxstrat + 2;
		    datlin[i,j] = 0;
		   end;
		  end;
		 end;

		* analogue for CE;
		%if &max_stat = 2 %then %do;
			use u6_;
			 read all var("_zeitp_" 
							%do is = 1 %to &maxstrat; || "_cif0&is._" %end;
							|| "_cif_") 
			  into datlin_;
			 close u6_;
			use v_CE;
			 read all var("&v1" %if &vz>1 %then %do; 
			   %do j = 2 %to &vz;  || "&&v&j" %end;
			  %end;) into par_;
			 close v_CE;
			ioparms[5] = nrow(datlin_); /* CE */
			%end;
		%else %do;
			datlin_ = -999;
			par_ = -999;
			ioparms[5] = 0;
			%end;

		x = time(); * reset computation time;

		 call explvar(par,par_,datlin,datlin_,ioparms);

		y = time();
		V_ind = ioparms[3];
		V_dir = ioparms[4];
		DN = ioparms[5];
		DS = ioparms[6];
		outset = V_dir||V_ind||DN||DS;

		z = y - x; * computation time;
		create &out from outset [colname = {'V_dir', 'V_ind', 'DN', 'DS'}];
			append from outset;
			close &out;
		create _ttt_ from z [colname = {'time'}];
			append from z;
			close _ttt_;
		quit; /* proc iml */

data _ttt_;
	set _ttt_;
	time = floor(time*100 + 0.5)/100;
	put "NOTE: Execution of surev used " time "seconds.";
	run;

data &out;
	set &out;
	run;

%if &print %then %do;
	proc print noobs label;
	var V_dir V_ind DN DS;
	run;
	%end;%mend;

/*
%NecSuff_CR(DATA=abcsg.data, TIME=tt_dist, STATUS=stat_dist, RESP=1, cr=2,
			varlist=	Nstage1 Nstage2 Nstage3);
%NecSuff_CR(DATA=abcsg.data, TIME=tt_dist, STATUS=stat_dist, RESP=1, cr=2,
			varlist=, 
		  	strat=GradingC__);
%NecSuff_CR(DATA=abcsg.data, TIME=tt_dist, STATUS=stat_dist, RESP=1, cr=2,
			varlist=	log_TumorSizeN);
%NecSuff_CR(DATA=abcsg.data, TIME=tt_dist, STATUS=stat_dist, RESP=1, cr=2,
			varlist=	Nstage1 Nstage2 Nstage3
						log_TumorSizeN, 
		  	strat=GradingC__);
*/

