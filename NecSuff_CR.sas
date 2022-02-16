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
* survival data. Submitted to SMMR.
*
*
* Author:  Andreas Gleiss
* Version: 1.0
* Date:    16 Feb 2022
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


%if &strat= %then %let strat=_all;
data _data;
	set &data;
	_all=1;
	run;

%if &varlist= %then %let varlist=_all;

%let nvar=0;
%do %while(%scan(&varlist,&nvar+1)~=);
 %let nvar=%eval(&nvar+1);
 %let v&nvar=%scan(&varlist,&nvar);
%end;
%let vz=&nvar;

proc means data=_data noprint;
	var &strat;
	output out=_maxstrat max=maxstrat;
	run;
data _null_;
	set _maxstrat;
	call symput('maxstrat',maxstrat);
	run;
%let maxstrat=%eval(&maxstrat);

DATA RX; SET _DATA; Z=&RESP; z_=&cr;
IF &STATUS=Z THEN STAT=1;  
			 ELSE if &status=z_ then stat=2; else STAT=0; /*!!!!!!*/
IF &TIME=. or STAT=. THEN DELETE;
IF &TIME<0 THEN DO;
   ERROR 'NEGATIVE TIME ENCOUNTERED, PROGRAM STOPS';
   ABORT; END;
if missing(&strat) then delete;
%DO J=1 %TO &VZ; IF missing(&&V&J) THEN DELETE;   %END;    NN+1;
run;
%if &useranks=1 %then %do;
 proc rank data=rx out=rx;
 var &time;
 ranks &time._r;
 run;
 %let timebu=&time;
 %let time=&time._r;
%end;

* SUBTRACT MEANS FROM COVARIATES AND SORT BY ASCENDING TIME;
PROC STANDARD DATA=RX OUT=RR MEAN=0;
VAR %DO J=1 %TO &VZ; &&V&J %END;  ;
run;

PROC SORT DATA=RR;
 BY &TIME;
run;

data rr;
set rr;
_line_=_n_;
run;

	data _covs;
		%DO J=1 %TO &VZ; 
			&&V&J=0; /* also in presence of strata: standardize covariates in all strata to zero */ 
			%END; 
		run;

PROC PHREG DATA=RR OUTEST=V noprint;               ***** coefficients ;
MODEL &TIME*STAT(0)= %DO J=1 %TO &VZ; &&V&J %END;  
	/ eventcode=1 &options; 
	BASELINE covariates=_covs OUT=S cif=cif0;
	strata &strat;
run;
	proc sort data=s out=_s;
		by &time;
		run;
	data s_;
		merge
		%do is=1 %to &maxstrat; 
			_s(keep=&time cif0 &strat where=(&strat=&is) rename=(cif0=cif0&is)) 
			%end;;
		by &time;
		drop &strat ;
		run;

	data s_;
		set s_;
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

PROC LIFETEST DATA=RR outcif=t noprint;
 TIME  &TIME*STAT(0) / eventcode=1; 
run;
data t;
	set t;
	survival=1-cif; 
	run;


PROC LIFETEST DATA=RR OUTSURV=B NOPRINT ;  ******* reverse KM estimate;
 TIME  &TIME*STAT(1,2); 
run;

proc freq data=rr noprint;
	tables stat / out=_stats;
	run;
proc means data=_stats noprint;
	var stat;
	output out=_stats2 min=smin;
	run;
data _null_;
	set _stats2;
	call symput('smin',smin);
	run;

ods select all;

DATA U3(KEEP= &TIME  KM cif) ;
 SET T;
 KM=SURVIVAL;
 cif=cif;
 IF KM ne .;
run;

DATA U2(KEEP= &TIME  KM cif) ;
 SET U3;
 BY &TIME;
 IF LAST.&TIME;
run;


DATA U1(KEEP= &TIME BKM) ;
 SET B; 
 if survival ne . then BKM=SURVIVAL;
 RETAIN BKM;
run;


DATA U;
 MERGE U1 U2(IN=W);
 BY &TIME;
 IF W;
run;

DATA u5;
 MERGE S_(IN=W) U ;
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
 _zeitp_=&time;
  %do is=1 %to &maxstrat;
	_cif0&is._=cif0&is; * baseline CIF;
	%end;
 _cif_=cif;	  * uncond. CIF;
 _bkm_=bkm;   * reverse KM;
 _line_=_n_;
 keep _line_ _zeitp_ _cif0: _cif_ _bkm_;
run;

proc means data=rr noprint; 
	var &time;
	output out=_maxt max=maxt;
	run;

data u7;
merge rr u6;
by _line_; 
run;

data _time_;
x=time();
output;
run;

		PROC IML;

				start explvar(PAR,DATALINES,IOPARMS);
					*** code translated from original fortran SUREV.FOR by G. H. 2012 Apr 12;

					zeit=repeat(0,ioparms[1]);
					exbeta=repeat(0,ioparms[1]);
					istat=repeat(0,ioparms[1]);
					inde=repeat(0,ioparms[1]);
					zeitp=repeat(0,ioparms[3]);
					cif0=j(ioparms[3],ioparms[4],0);	* baseline CIF;
					cif=repeat(0,ioparms[3]);			* uncond. CIF;
					bkm=repeat(0,ioparms[3]);			* reverse KM;
					dup=repeat(0,ioparms[3]);
				    cov=j(ioparms[1],ioparms[2],0);

				   n=ioparms[1];  * # individuals;
				   nk=ioparms[2]; * # covariables;
				   np=ioparms[3]; * # distinct failure times for event of interest;
				   ns=ioparms[4]; * # strata;
				    
				   zeit=datalines[,1];
				   istat=datalines[,2];
				   cov=datalines[,3:(nk+2)];
				   strt=datalines[,nk+3];
				   zeitp=datalines[,nk+4];
				   cif0=datalines[,(nk+5):(nk+5+ns-1)];
				   cif=datalines[,(nk+5+ns)];
				   bkm=datalines[,(nk+5+ns+1)];

				   * end of all input;

				   	  * construction of pointer (inde) and of counter of tied failure times of interest/*!!!!!!*/ (dup);
				      IZ=0;                                                             
				      do j = 1 to np; * loop over distinct times to event of interest t(j);
				      	  dup[j]=0;                                                      
						  mark31:
					      if iz = n then goto mark32;                                            
						  iz=iz+1;                                                           
					      if istat[iz]=0 then do;                                        
						      inde[iz]=j;                                                        
						      goto mark31;                                                           
					      end;                                                              
					      if istat[iz]=1 & zeit[iz]=zeitp[j] then do; 
						      dup[j]=dup[j]+1;                                                   
						      inde[iz]=j;                                                        
						      goto mark31;                                                           
					      end;                                                       
					      if istat[iz]=2 then do; 
						      inde[iz]=j;                                                        
						      goto mark31;                                                           
					      end;                                                              
					      if zeit[iz] > zeitp[j] then iz=iz-1; 
				          mark30:
				      end;  
					  mark32:
					  d1u_s=0; * for dir. est. of V;                                                      
				      d1c_s=0;                                                       
				      wei_s=0;                                                        
					  dn=0; * variant 1;
		  			  ds=0; * variant 1;


					  do i=1 to n;                                                      
					      xbeta=0;                                                      
					      do i1=1 to nk; 
							 mark19:
							 xbeta=xbeta+cov[i,I1]*par[i1];                                    
					      end;
					   	  mark18:
						  exbeta[i]=exp(xbeta);
				      end; 
						                            
				      do j=1 to np; * BIG LOOP OVER DISTINCT EoI FAILURE TIMES; 
					      resc=0;                                                           
					      resk=0;     
						  q_sum=0; 
						  term_vw=0;
			  			  term_n=0; 
			  			  term_s=0; 
						  q_kl=0; q_gr=0;
					      resc_ind=0;                                                           
					      resk_ind=0;     

					      do i=1 to n; * loop over individuals;    
							* direct estimate of V;
						      fg=1-(1-cif0[j,strt[i]])**exbeta[i]; 

						      if zeit[i] > zeitp[j] then do; * alive at t(j);                                   
							      resc=resc+fg;                                                 
							      resk=resk+cif[j];
								  q=1; 
							      goto mark12;                                                          
						      end;                                                             
						      if zeit[i] <= zeitp[j] & istat[i] =1 then do; * event of interest at t(j);
							      resc=resc+(1-fg);                                                  
							      resk=resk+(1-cif[j]);      
								  q=1; 
							      goto mark12;                                                           
						      end;                                                             
						      if zeit[i] <= zeitp[j] & istat[i] = 0 then do; * censored at t(j);             
							      ife=inde[i];                                                  
							      fgc=1-(1-cif0[ife,strt[i]])**exbeta[i];
								  q=1; 
							      if fgc < 0.999999 then RESC=RESC+(fg*(1-fg)/(1-fgC) + (1-fg)*(1-(1-fg)/(1-fgC))); 
							      if cif[ife] > 0.000001 then  RESK=RESK+(cif[J]*(1-cif[J])/(1-cif[IFE]) + (1-cif[J])*(1-(1-cif[J])/(1-cif[IFE]))); 
						      end;                     
						      if zeit[i] <= zeitp[j] & istat[i] = 2 then do;
							      ife=inde[i];                                                  
							      q=bkm[j]/bkm[ife]; 
							      RESC=RESC +q*fg; 
							      RESK=RESK +q*cif[J]; 
						      end;  
					          mark12:

							* estimate of DN, DS;
							  if fg<cif[j] then do; /* >>> NECESSARY <<< */
							    q_kl=q_kl+q;
								if cif[j]<1 then 
									term_n=term_n + q*((fg-cif[j])/cif[j])**2;
								end;
							  if fg>cif[j] then do; /* >>> SUFFICIENT <<< */
							    q_gr=q_gr+q;
								if cif[j]>0 then 
									term_s=term_s + q*((cif[j]-fg)/(1-cif[j]))**2;
								end;

							  q_sum=q_sum+q; 


				          end; /* i-loop */   
						

						  D1C_S=D1C_S+(RESC/q_sum)*DUP[J]/BKM[J]; * numerator for dir. V;
					      D1U_S=D1U_S+(RESK/q_sum)*DUP[J]/BKM[J]; * denom. for dir. V;
					      WEI_S=WEI_S+DUP[J]/BKM[J];                                                                              

						  if q_kl>0 then 
					      	dn=dn+sqrt(term_n/q_kl)*DUP[J]/BKM[J];                                
						  if q_gr>0 then
					      	ds=ds+sqrt(term_s/q_gr)*DUP[J]/BKM[J];                                

					      mark8:
				      end; /* j-loop */                                                      

				      D=D1U_S/WEI_S;                                          
				      Dx=D1C_S/WEI_S;                                         
				      ioparms[4]=(D-Dx)/D; /* direct estimate of V */

				      ioparms[5]=dn/wei_s;                                          
				      ioparms[6]=ds/wei_s;
					  if ioparms[4]<0 then ioparms[4]=0; 
					  if ioparms[5]<0 then ioparms[5]=0; 
					  if ioparms[6]<0 then ioparms[6]=0; 

					  mark546:

				finish; /* explvar */



		use u7;
		 read all var("&time" || "stat"  %do j=1 %to &vz;  || "&&v&j" %end; || "&strat"	|| "_zeitp_" 
						%do is=1 %to &maxstrat; || "_cif0&is._" %end;
						|| "_cif_"
						||"_bkm_") 
		  into datlin;
		 close u7;
		use u6;
		 read all var("_zeitp_") into u6;
		 close u6;
		use v;
		 read all var("&v1" %if &vz>1 %then %do; 
		   %do j=2 %to &vz;  || "&&v&j" %end;
		  %end;) into par;
		 close v;

		ioparms=repeat(0,6,1);
		ioparms[1]=nrow(datlin);
		ioparms[2]=&vz;
		ioparms[3]=nrow(u6);
 		ioparms[4]=&maxstrat;
		if ioparms[3]<ioparms[1] then do;
		  do i=ioparms[3]+1 to ioparms[1];
		   do j=&vz+4 to &vz+4+&maxstrat+2;
		    datlin[i,j]=0;
		   end;
		  end;
		 end;

		x=time();

		call explvar(par,datlin,ioparms);

		y=time();
		V_dir=ioparms[4];
		DN=ioparms[5];
		DS=ioparms[6];
		outset=V_dir||DN||DS;

		z=y-x;
		create &out from outset [colname={'V_dir','DN','DS'}];
		append from outset;
		close &out;
		create _ttt_ from z [colname={'time'}];
		append from z;
		close _ttt_;
		quit; /* proc iml */

data _ttt_;
set _ttt_;
time=floor(time*100+0.5)/100;
put "NOTE: Execution of surev used " time "seconds.";
run;

data &out;
set &out;
run;

%if &print %then %do;
	proc print noobs label;
	var V_dir DN DS;
	run;
	%end;
%mend;

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

