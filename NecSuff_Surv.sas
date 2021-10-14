*******************************************************************;
*
* NECSUFF_SURV
* ============
*
* SAS-macro for computation of degrees of necessity and sufficiency
* for a survival outcome variable
*
* This macro has been adapted from the macro %relipmcr by Heinze & 
* Schemper (2003) that calculates the direct V measure of explained 
* variation as defined in Schemper & Henderson (2000)
*
* Output: marginal and partial measures are reported for:
*		  DN and DS variant 1 (see Gleiss & Schemper, 2019),
*		  EV, i.e. indirect V_W measure of explained variation
*			   as defined in Schemper & Henderson (2000)
*
* Gleiss, A. & Schemper, M., Quantifying degrees of necessity and 
* sufficiency in cause-effect relationships with dichotomous and 
* survival outcome.
* Statistics in Medicine 2019, 38:4733–4748.
*
* Schemper M & Henderson R., Predictive accuracy and explained 
* variation in Cox regression. 
* Biometrics 2000, 56(1):249-255.
*
* Heinze, G. & Schemper, M., Comparing the importance of 
* prognostic factors in Cox and logistic regression using SAS. 
* Computer Methods and Programs in Biomedicine 2003, 71:155-163.
*
*
* Author:  Andreas Gleiss
* Version: 1.0
* Date:    15 Jan 2020
*
* Macro parameters:
* =================
* 
* data		name of SAS data set
* time		name of time variable of survival outcome 
* cens		name of censoring variable of survival outcome 
* censval	value of censoring variable indicating censoring
*			(default=0)
* varlist	list of independent variables
* cand		list of candidate variables
* class		list of class variables
* nboot		number of bootstrap replicates (=1 no bootstrap)
*			(bootstrap comparison only relating to EV)
* shrink	indicator for shrinking parameter estimates (1=yes)
* firth		indicator for using Firth-type bias correction (1=yes)
* seed		seed for bootstrap random number generation
* correctp	indicator for Tukey correction of bootstrap p-values
* s0meth	method for estimating baseline survival function
* alpha		significance level (bootstrap only)
* clp		confidence level (bootstrap only)
* groups	specify which variables belong to which group
* gnames	list of group names (see example below)
* printclp	indicator to print CI for bootstrap p-values
* histogram indicator to show histogram of bootstrapped values
* useranks	if =1 replaces survival times by their ranks (default)
* varlen	maximum length of variable names in output
*
* (for option details see Heinze & Schemper, 2003)
*
*******************************************************************;


%macro NecSuff_Surv(data=_last_, time=time, cens=cens, censval=0, varlist=, cand=, class=, 
		  nboot=1, shrink=0, firth=0, seed=395748, 
          correctp=0, s0meth=pl, alpha=0.05, clp=0.99, groups=, gnames=, printclp=0, histogram=0, 
          useranks=1, varnlen=15);

			%let Version= 1_0;
			%let Build=   202001101600;


%macro srcore(srcdata=_srcwork, srctime=&time, srccens=&cens, srck=1, rby=_rby, srcout=_srcout, pevname=pev);

* core program to compute PEV for a data set in by-groups       ;
* covariate names must be _src1 to _src&srck                    ;
* by variable must be numeric integer, starting with 1          ;
* and no missing values, please                                 ;

%if &srck=0 %then %do;
	 data &srcout;
	 _sample_=0;
	 &pevname=0;
	 output;
	 run;
	%end;
%else %do;

	options nonotes nosource;
	*options nomprint nomacrogen;
	*filename SASCBTBL 'T:\SAS\sasdll\surev.def';

	data _ttt_;
	x=time();
	output;
	run;

	/*%if &firth=0 %then %do;*/
	 %if &shrink=1 %then %do;
	  ods select none;
	  ods output GlobalTests=_glob;
	 %end;
	 proc phreg data=&srcdata outest=_srcv   /* --> par */
	 %if &shrink=0 %then %do; 
	  noprint 
	 %end; ;
	 model &srctime*&srccens(0)=%do srj=1 %to &srck; _src&srj %end;  %if &firth=1 %then %do; /firth %end; ;
	 %if &shrink = 0 %then %do;
	  baseline out=_srcs survival=s0 /method=&s0meth;
	 %end;
	 by &rby;
	 run;
	 %if &shrink=1 %then %do;
	  ods select all;
	  data _glob;
	  set _glob; 
	   if Test="Likelihood Ratio";
	   _sf_=(ChiSq-&srck)/ChiSq;
	   keep &rby _sf_;
	  run;

	  data _srcv;
	  merge _srcv _glob;
	  by &rby;
	  %do srj=1 %to &srck;
	   _src&srj=_src&srj*_sf_;
	  %end;
	  run;
	 %end;

	/*%end;*/


	/*%if %eval(&shrink+&firth)>0 %then %do;*/
	 %if &shrink=1 %then %do;
	 data &srcdata;
	 merge &srcdata _srcv (rename=(
	  %do srj=1 %to &srck;
	   _src&srj=_b&srj
	  %end;
	  ));
	 by &rby;
	  _lp_=0;
	  %do srj=1 %to &srck;
	   _lp_=_lp_+_src&srj*_b&srj;
	  %end;
	 run;
	 proc phreg noprint data=&srcdata;
	 model &srctime*&srccens(0)=/offset=_lp_ maxiter=0;
	 baseline out=_srcs survival=s0 /method=&s0meth;
	 by &rby;
	 run;
	%end;


	DATA _srcU4;
	 MERGE _srcS(IN=W) _srcU ;
	 BY  &rby &srctime;
	 IF W;
	run;

	DATA _srcU5;
	 SET _srcU4;
	 RETAIN BKMLAG KMLAG;
	 IF BKM=. THEN BKM=BKMLAG;
	 IF  KM=. THEN BKM= KMLAG;
	 OUTPUT;
	 BKMLAG=BKM;
	 KMLAG= KM;
	RUN;

	DATA _srcu6;
	 SET _srcU5 ;
	 BY  &rby &srctime;
	 IF FIRST.&srctime;
	 IF &srctime > 0.0; 
	run;

	proc means data=_srcu6 noprint;
	var &srctime;
	output out=_nkm n=nkm;
	by &rby;
	run;


	proc iml;
	 
	 * read survival curves;
	start explvar(PAR,DATALINES,IOPARMS);
		*** code translated from original fortran SUREV.FOR by G. H. 2012 Apr 12;

		zeit=repeat(0,ioparms[1]);
		exbeta=repeat(0,ioparms[1]);
		istat=repeat(0,ioparms[1]);
		inde=repeat(0,ioparms[1]);
		zeitp=repeat(0,ioparms[3]);
		s0=repeat(0,ioparms[3]);
		skm=repeat(0,ioparms[3]);
		bkm=repeat(0,ioparms[3]);
		dup=repeat(0,ioparms[3]); * tied failure times;
	    cov=j(ioparms[1],ioparms[2],0);

	   n=ioparms[1]; 	* #obs;
	   nk=ioparms[2]; 	* #covar;
	   np=ioparms[3]; 	* #distinct failure times;
	    
	   zeit=datalines[,1];
	   istat=datalines[,2];
	   cov=datalines[,3:(nk+2)];
	   zeitp=datalines[,nk+3];
	   s0=datalines[,nk+4];
	   skm=datalines[,nk+5];
	   bkm=datalines[,nk+6];
	  
	* end of all input;
	* construction of pointer (inde) and of counter of tied failure times (dup);


	      IZ=0;                                                             
	      do j = 1 to np;*DO 30 J=1,NP;  
	      	  dup[j]=0; *DUP(J)=0.;                                                         
			  mark31:
		      if iz = n then goto mark32; *IF (IZ .EQ. N) GOTO 32;                                            
			  iz=iz+1; *31 IZ=IZ+1;                                                           
		      if istat[iz]=0 then do; *IF (ISTAT(IZ) .EQ. 0) THEN;                                        
			      inde[iz]=j; *INDE(IZ)=J;                                                        
			      goto mark31; *GOTO 31;                                                           
		      end; *ENDIF;                                                             
		      if zeit[iz]=zeitp[j] then do; *IF (ZEIT(IZ) .EQ. ZEITP(J)) THEN;                                  
			      dup[j]=dup[j]+1; *DUP(J)=DUP(J)+1.;                                                  
			      inde[iz]=j; *INDE(IZ)=J;                                                        
			      goto mark31; *GOTO 31;                                                           
		      end; *ENDIF;                                                             
		      if zeit[iz] > zeitp[j] then iz=iz-1; *IF (ZEIT(IZ) .GT. ZEITP(J)) IZ=IZ-1;
	          mark30:
	      end; *30 CONTINUE;                                                          
		  mark32:
		  d1u=0;                                                         
	      d1c=0;                                                          
	      wei=0; *WEI_S=0.;  
		  ev=0;
		  dn=0; * variant 1;
		  ds=0; * variant 1;
	      do i=1 to n; *DO 18 I=1,N;                                                       
		      xbeta=0; *XBETA=0.;                                                          
		      do i1=1 to nk; *DO 19 I1=1,NK;
				 mark19:
				 xbeta=xbeta+cov[i,I1]*par[i1]; * 19 XBETA=XBETA+COV(I,I1)*PAR(I1);                                     
		      end;
		   	  mark18:
			  exbeta[i]=exp(xbeta); *18 EXBETA(I)=DEXP(XBETA);
	      end; 

		  *C     BIG LOOP OVER DISTINCT FAILURE TIMES  ;
			          
	      do j=1 to np; *DO 8 J=1,NP;                                                       
		      resc=0; *RESC=0.;                                                           
		      resk=0; *RESK=0.; 
			  term_ev=0;
			  term_n=0; 
			  term_s=0; 

			  n_kl=0; n_gr=0;

		      do i=1 to n; *DO 12 I=1,N;   
				      cox=s0[j]**exbeta[i];*COX=S0(J)**EXBETA(I);   
					  resc=resc+2*cox*(1-cox);
					  resk=resk+2*skm[j]*(1-skm[j]);
					  if skm[j]>0.000001 & skm[j]<0.99999 then
						term_ev=term_ev + (1-(cox*(1-cox))/(skm[j]*(1-skm[j])));
					  if cox>skm[j] then do; /* >>> NECESSARY <<< */
					    n_kl=n_kl+1;
						if skm[j]<1 then /* falls skm[j]=1 => cox>1 in diesem Teil => gar nix */
							term_n=term_n + ((cox-skm[j])/(1-skm[j]))**2;
						end;
					  if cox<skm[j] then do; /* >>> SUFFICIENT <<< */
					    n_gr=n_gr+1;
						if skm[j]>0 then /* falls skm[j]=0 => cox=0 in diesem Teil => gar nix */
							term_s=term_s + ((skm[j]-cox)/(skm[j]))**2;
						end;
							/* * originally in %relimpcr: SH measure = direct V;
						      if zeit[i] > zeitp[j] then do; *IF (ZEIT(I) .GT. ZEITP(J)) THEN;                                   
							      resc=resc+(1-cox); *RESC=RESC+(1.-COX);                                                
							      resk=resk+(1-skm[j]);*RESK=RESK+(1.-SKM(J));                                             
							      goto mark12; *GOTO 12 ;                                                          
						      end; *ENDIF;                                                            
						      if zeit[i] <= zeitp[j] & istat[i] ^= 0 then do; *IF (ZEIT(I) .LE. ZEITP(J) .AND. ISTAT(I) .NE. 0) THEN;             
							      resc=resc+cox; *RESC=RESC+COX;                                                     
							      resk=resk+skm[j]; *RESK=RESK+SKM(J);                                                  
							      goto mark12; *GOTO 12;                                                           
						      end; *ENDIF;                                                             
						      if zeit[i] <= zeitp[j] & istat[i] = 0 then do; *IF (ZEIT(I) .LE. ZEITP(J) .AND. ISTAT(I) .EQ. 0) THEN;             
							      ife=inde[i]; *IFE=INDE(I);                                                       
							      coxc=s0[ife]**exbeta[i]; *COXC=S0(IFE)**EXBETA(I);                                           
							      if coxc > 0.000001 then RESC=RESC+((1-COX)*COX/COXC+COX*(1-(COX/COXC))); *IF (COXC .GT. 0.000001) RESC=RESC+((1.-COX)*COX/COXC+COX*(1.-(COX/COXC)));                 
							      if skm[ife] > 0.000001 then  RESK=RESK+((1-SKM[J])*SKM[J]/SKM[IFE]+SKM[J]*(1-(SKM[J]/SKM[IFE]))); *IF (SKM(IFE) .GT. 0.000001) RESK=RESK+((1.-SKM(J))*SKM(J)/SKM(IFE)+SKM(J)*(1.-(SKM(J)/SKM(IFE))));                                               
						      end; *ENDIF;    
							*/ 

			          mark12:

	          end; *12 CONTINUE;  * Ende i-Schleife;                                                      

		      *D1U=D1U+(RESK/N)*DUP[J]/BKM[J]; * Nenner für indirektes V; 
		      *D1C=D1C+(RESC/N)*DUP[J]/BKM[J]; * Zähler für indirektes V;
			  ev=ev+(term_ev/n)*DUP[J]/BKM[J]; * indirektes V_W;
		      WEI=WEI+DUP[J]/BKM[J];     * Gewichtersumme;  
			  if n_kl>0 then 
		      	dn=dn+sqrt(term_n/n_kl)*DUP[J]/BKM[J];                                
			  if n_gr>0 then
		      	ds=ds+sqrt(term_s/n_gr)*DUP[J]/BKM[J];                                

		      mark8:
	      end;*8 CONTINUE;   * Ende j-Schleife;                                                       

	*print wei;

	      *ioparms[4]=1 - D1C / D1U;   * indirect V (Gewichtersumme wei kürzt sich raus); 
	      ioparms[4]=ev/wei;   		  * indirect V_W;                                
	      ioparms[5]=dn/wei;                                          
	      ioparms[6]=ds/wei;
	  if ioparms[4]<0 then ioparms[4]=0; 
	  if ioparms[5]<0 then ioparms[5]=0; 
	  if ioparms[6]<0 then ioparms[6]=0; 
		  mark546:
		   *546   CONTINUE;

	finish; /* explvar */


	 use _srcu6;
	  read all var("&srctime" || "s0" || "km" ||"bkm") into kmall;
	 close _srcu6;
	 
	 * read data;

	 use &srcdata;
	  read all var("&srctime" || "&srccens" %do srj=1 %to &srck; || "_src&srj" %end;) into dlall;
	 close &srcdata;
	 
	 * read parameter estimates;

	 use _srcv;
	  read all var("_src1" %if &srck>1 %then %do; 
	   %do srj=2 %to &srck;  || "_src&srj" %end;
	  %end;) into parall;
	 close _srcv;

	 * read n_obs for all datasets;

	 use _ndat;
	  read all var("ndat") into ndat;
	 close _ndat;

	 * read n_t (number of distinct uncensored survival times) for all datasets;

	 use _nkm;
	  read all var("nkm") into nkm;
	 close _nkm;

	 x=time();
	 pev=repeat(0,nrow(nkm),3);
	 index=0;
	 indexkm=0;
	 do rby=1 to nrow(nkm);

	  datlin=repeat(0,ndat[rby],&srck+6);

	  do i=1 to ndat[rby];
	   index=index+1;
	   do j=1 to &srck+2;
	    datlin[i,j]=dlall[index,j];
	   end;
	   if i<=nkm[rby] then do;
	    indexkm=indexkm+1;
	    do j=1 to 4;
	     datlin[i,&srck+2+j]=kmall[indexkm,j];
	    end;
	   end;
	  end;

	  ioparms=repeat(0,6,1);
	  ioparms[1]=ndat[rby];    * number of observations;
	  ioparms[2]=&srck;          * number of covariates;
	  ioparms[3]=nkm[rby];       * number of distinct uncensored survival times;
	  if ioparms[3]<ioparms[1] then do;
	   do i=ioparms[3]+1 to ioparms[1];
	    do j=&srck+3 to &srck+6;
	     datlin[i,j]=0;
	    end;
	   end;
	  end;
	  par=parall[rby,];
	*  print rby ioparms, par, datlin;
	*  CALL modulei('*E','SUREV',par,datlin,ioparms);
	*  print par, datlin, ioparms;
	*  stop;
	  call explvar(par,datlin,ioparms);
	  ev=ioparms[4];
	  dn=ioparms[5];
	  ds=ioparms[6];
	  pev[rby,1]=ev; pev[rby,2]=dn; pev[rby,3]=ds;
	  *print "&title ","EV = " ev,
	        "DN = " dn,
	        "DS = " ds;
	 end;
	 create &srcout from pev[colname={"&pevname","&pevname._DN","&pevname._DS"}];
	 append from pev;
	 close &srcout;

	 create &srcout.par from parall; /*!!!!!!*/
	 append from parall;
	 close &srcout.par;

	quit; * Ende IML;

	data &srcout;
	set &srcout;
	_sample_=_n_-1;
	run;

	options notes source;
	*options mprint macrogen;
	%end; /* zu %else %do ganz oben */

data _ttt_;
set _ttt_;
time=time()-x;
time=floor(time*100+0.5)/100;
put "NOTE: Execution of srcore used " time "seconds.";
run;


%mend; * srcore;







%let clpval=%sysevalf(&clp*100);

%let nvar=0;
%if &varlist ne %then %do;
 %do %while(%scan(&varlist,&nvar+1)~=);
  %let nvar=%eval(&nvar+1);
  %let var&nvar=%scan(&varlist,&nvar);
 %end;
%end;

%let ngroup=0;
%if &groups ne %then %do;
 %do k=1 %to &nvar;
  %let group&k=%scan(&groups,&k);
  %if &&group&k gt &ngroup %then %do;
   %let ngroup=%eval(&&group&k);
  %end;
 %end;
 %do k=1 %to &ngroup;
  %if &gnames ne %then %do;
   %let gname&k=%scan(&gnames,&k);
  %end;
  %else %do;
   %let gname&k=&k;
  %end;
 %end;
%end;

data _xxx;
put "NOTE: &ngroup groups found.";
time=time();
output;
run;

%let ncand=0;
%if &cand ne %then %do;
 %do %while(%scan(&cand,&ncand+1)~=);
  %let ncand=%eval(&ncand+1);
  %let cand&ncand=%scan(&cand,&ncand);
 %end;
%end;

%if &ncand ne 0 %then %do;
 %let mode=CANDIDATE;
 %let nobj=&ncand;
 %do k=1 %to &nobj;
  %let obj&k=&&cand&k;
 %end;
%end;
%else %if &ngroup ne 0 %then %do;
 %let mode=GROUP;
 %let nobj=&ngroup;
 %do k=1 %to &ngroup;
  %let obj&k=&&gname&k;
 %end;
%end;
%else %if &varlist ne  %then %do;
 %let mode=TOTAL;
 %let nobj=&nvar;
  %do k=1 %to &nobj;
  %let obj&k=&&var&k;
 %end;
%end;

%if &nobj=1 %then %do;
 %let nboot=1;
 %let printclp=0;
 %let histogram=0;
%end;

data _worksr;
set &data;
if &time ne .;
if &cens ne .;
if &cens=&censval then &cens=0;
else &cens=1;
%if &varlist ne %then %do;
 %do j=1 %to &nvar;
  if &&var&j ne .;
 %end;
%end;
%if &cand ne %then %do;
 %do j=1 %to &ncand;
  if &&cand&j ne .;
 %end;
%end;
run;

%if &useranks=1 %then %do;
 proc rank data=_worksr out=_worksr;
 var &time;
 ranks _rtime_;
 run;
 %let time = _rtime_;
%end;

%let nclass=0;

%if &class ne %then %do; 

 %do %while(%scan(&class,&nclass+1)~=);
  %let nclass=%eval(&nclass+1);
  %let class&nclass=%scan(&class,&nclass);
 %end; 


 proc freq data=_worksr;
 %do j=1 %to &nclass;
  tables &&class&j / out=_clt&j noprint;
 %end;
 run;

 proc iml;
 use _worksr;
 read all var("&class1" 
  %if &nclass>1 %then %do;
   %do j=2 %to &nclass;
    || "&&class&j"
   %end;
  %end; ) into x;
 close _worksr;
 newvar=0;
 n=nrow(x);
 df=repeat(0,1,&nclass);
 %do j=1 %to &nclass;
  use _clt&j;
  read all var("&&class&j") into factor&j;
  close _clt&j;
  df[&j]=nrow(factor&j)-1;
  newvar=newvar+df[&j];
 %end;

 newx=repeat(0,n,newvar); 

 do i=1 to n;
  index=0;
  %do j=1 %to &nclass;
   do ii=1 to df[&j];
    index=index+1;
    newx[i,index]=(x[i,&j]=factor&j[ii+1]);
   end;
  %end;
 end;

 create _df from df[colname={%do j=1 %to &nclass; "df&j" %end;}];
 append from df;
 close _df; 

 create _dummies from newx;
 append from newx;
 close _dummies;
 

 quit; 

 data _df;
 set _df;
  %do j=1 %to &nclass;
   call symput("df&j", df&j);
  %end;
 run; 


 data _dummies;
 set _dummies;
  %let index=0;
  %do j=1 %to &nclass;
   %do jj=1 %to &&df&j;
    %let index=%eval(&index+1);
    rename col&index=_dum&j._&jj;
   %end;
  %end;
 run;

 data _worksr;
 merge _worksr _dummies;
 run;

%end;  *class;
 
 *** mapping of factors to dummies ***;

 %if &varlist ne %then %do;
  %do j=1 %to &nvar;
   %let dumv&j=0;
   %do k=1 %to &nclass;
    %if &&var&j=&&class&k %then %do;
	 %let dumv&j=&k;
	%end;
   %end;
  %end;
 %end;
 %if &cand ne %then %do;
  %do j=1 %to &ncand;
   %let dumc&j=0;
   %do k=1 %to &nclass;
    %if &&cand&j=&&class&k %then %do;
	 %let dumc&j=&k;
	%end;
   %end;
  %end;
 %end;


  ********* erzeuge grossen bootstrap-Datensatz **********;
  ********* balanced bootstrap                  **********;

%if &nboot=1 %then %do;
  data _bootsr;
  set _worksr;
   _sample_=0;
  run;
 %end;
				 %else %do;
				  proc iml; 
				  use _worksr; 
				  read all var("&time"||"&cens"
				   %if &varlist ne %then %do;
				    %do k=1 %to &nvar;
					 %if &&dumv&k=0 %then %do;
				      ||"&&var&k"
					 %end;
					 %else %do;
					  %let help=&&dumv&k;
					  %do kk=1 %to &&df&help;
					   || "_dum&help._&kk"
					  %end;
					 %end;
				    %end; 
				   %end;
				   %if &cand ne %then %do;
				    %do k=1 %to &ncand;
					 %if &&dumc&k = 0 %then %do;
				      ||"&&cand&k"
					 %end;
					 %else %do;
					  %let help=&&dumc&k;
					  %do kk=1 %to &&df&help;
					   || "_dum&help._&kk"
					  %end;
				     %end; 
				    %end; 
				   %end;
				   ) into data; 
				  close _worksr;
				  n=nrow(data);
				  bootdata=repeat(0,n*&nboot,ncol(data)+1);
				  longrow=repeat(0,n*&nboot,1);
				  index=1;
				  do iboot=1 to &nboot;
				   do i=1 to n;
				    longrow[index]=i;
					index=index+1;
				   end;
				  end;
				  m=n*&nboot;
				  do iboot=1 to &nboot;
				   do i=1 to n;
				    draw=floor(ranuni(&seed)*m)+1;
				    take=longrow[draw];
				    longrow[draw]=longrow[m];
					m=m-1;
				    bootdata[(iboot-1)*n+i,1]=iboot;
				    do j=1 to ncol(data);
				     bootdata[(iboot-1)*n+i,j+1]=data[take,j];
				    end;
				   end;
				  end;
				  create _bootsr from bootdata[colname={"_sample_","&time","&cens"
				   %if &varlist ne %then %do;
				    %do k=1 %to &nvar;
					 %if &&dumv&k=0 %then %do;
				      ,"&&var&k"
					 %end;
					 %else %do;
					  %let help=&&dumv&k;
					  %do kk=1 %to &&df&help;
					   ,"_dum&help._&kk"
					  %end;
					 %end;
				    %end; 
				   %end;
				   %if &cand ne %then %do;
				    %do k=1 %to &ncand;
					 %if &&dumc&k = 0 %then %do;
				      ,"&&cand&k"
					 %end;
					 %else %do;
					  %let help=&&dumc&k;
					  %do kk=1 %to &&df&help;
					   ,"_dum&help._&kk"
					  %end;
				     %end; 
				    %end; 
				   %end;
				   }];
				  append from bootdata;
				  close _bootsr;
				  quit;
				  data _worksr;
				  set _worksr;
				  _sample_=0;
				  run;
				  data _bootsr;
				  set _worksr _bootsr;
				  run;
				 %end;
  
  
  *************** sorge vor fuer srcore ****************;

  PROC SORT DATA=_bootsr;
  BY _sample_ &time;
  run;

 %let varlist2=; 
 %if &varlist ne %then %do;
  %do j=1 %to &nvar;
   %let jump=0;
   %if &nclass gt 0 %then %do;
    %do jj=1 %to &nclass;
     %if &&class&jj=&&var&j %then %do;
	  %let jump=1;
      %do jjj=1 %to &&df&jj;
       %let varlist2=&varlist2 _dum&jj._&jjj;
      %end;
     %end;
	%end;
   %end;
   %if &jump=0 %then %do;
 	%let varlist2=&varlist2 &&var&j;
   %end;
  %end;
 %end;
 %if &cand ne %then %do;
  %do j=1 %to &ncand;
   %let jump=0;
   %if &nclass gt 0 %then %do;
    %do jj=1 %to &nclass;
     %if &&class&jj=&&cand&j %then %do;
	  %let jump=1;
      %do jjj=1 %to &&df&jj;
       %let varlist2=&varlist2 _dum&jj._&jjj;
      %end;
     %end;
	%end;
   %end;
   %if &jump=0 %then %do;
 	%let varlist2=&varlist2 &&cand&j;
   %end;
  %end;
 %end;

  PROC STANDARD DATA=_bootsr OUT=_bootsr1 MEAN=0;
  VAR &varlist2;
  by _sample_;
  run;

  proc means data=_bootsr1 noprint;
  var &time;
  output out=_ndat n=ndat;
  by _sample_;
  run;
  
  ods graphics off;
  proc lifetest data=_bootsr1 outsurv=_srct noprint;
  time &time*&cens(0);
  by _sample_;
  run;

  proc lifetest data=_bootsr1 outsurv=_srcb noprint;
  time &time*&cens(1);
  by _sample_;
  run;
  ods graphics on;


  DATA _srcu3(KEEP=_sample_ &time  KM) ;
  SET _srct;
   KM=SURVIVAL;
   IF KM ne . AND _CENSOR_=0;
  run;

  DATA _srcu2(KEEP=_sample_ &time  KM) ;
  SET _srcU3;
   BY _sample_ &time;
   IF LAST.&time;
  run;

  DATA _srcu1(KEEP=_sample_ &time BKM) ;
  SET _srcb;
   if survival ne . then BKM=SURVIVAL;
*  bkm=survival;
   *BY &TIME;
   RETAIN BKM;
  run;

  DATA _srcu;
  MERGE _srcU1 _srcU2(IN=W);
   BY _sample_ &time;
   IF W;
  run;

%if &nboot=1 %then %do;
 data _ndat;
 set _ndat;
  call symput("n",ndat);
 run;
%end;


%macro po(a,c,b);
 _option_="&a";
 _meaning_="&c";
 _value_="&b";
 output;
%mend;


data _settings;
length _option_ $10 _meaning_ $42 _value_ $ 50;
label _option_="Macro option" _meaning_="Description" _value_="Value";
%po(version, Version, &version);
%po(build, Build time stamp, &build);
%po(---, ---, ---);
%po(data, data set, &data);
%po(time, time variable, &time);
%po(cens, censoring indicator, &cens);
%po(censval, censoring value, &censval);
%po(s0meth, Method for baseline survival function, &s0meth);
%po(varlist, Variable list, &varlist);
%po(class, Class variables, &class);
%po(groups, Group belongings, &groups);
%po(gnames, Group names, &gnames);
%po(cand, Candidate list, &cand);
%po(firth, Firth option, &firth);
%po(shrink, Shrinkage, &shrink);
%po(nboot, Number of bootstrap resamples, &nboot);
%po(seed, Seed for resampling, &seed);
%po(correctp, Correct p-values by Tukeys HSD, &correctp);
%po(clp, Confidence level for p-values, &clp);
%po(alpha, Significance level,&alpha);
%po(printclp, Print confidence intervals, &printclp);
%po(histgoram, Plot histogram of PEV values, &histogram);
%po(useranks, Rank-transform survival times, &useranks);
%po(varnlen, Maximum length of variable names in output, &varnlen);
run;

proc print noobs label;
title3 "RELIMPCR macro";
title4 "Overview of selected options";
run;
title3;

%if &cand ne %then %do;

**********************************************************************;
****************** candidate mode ************************************;
**********************************************************************;

  ******************* marginale PEVs *****************************;
  %do j=1 %to &ncand;
   data _take_;
   set _bootsr1;
    %let jump=0;
    %do jj=1 %to &nclass;
     %if &&class&jj=&&cand&j %then %do;
	  %let jump=&&df&jj;
      %do jjj=1 %to &&df&jj;
       _src&jjj=_dum&jj._&jjj;
      %end;
     %end;
	%end;
	%if &jump=0 %then %do;
	 _src1=&&cand&j;
	 %let jump=1;
	%end;
   run;

   %srcore(srcdata=_take_, srck=&jump, rby=_sample_, srcout=_pevmc&j, pevname=_pevmc&j);  

  %end;
  %if &varlist ne %then %do;
   %do j=1 %to &nvar;
    data _take_;
    set _bootsr1;
	if _sample_=0;
     %let jump=0;
     %do jj=1 %to &nclass;
      %if &&class&jj=&&var&j %then %do; 
 	   %let jump=&&df&jj;
       %do jjj=1 %to &&df&jj;
        _src&jjj=_dum&jj._&jjj;
       %end;
      %end;
  	 %end;
	 %if &jump=0 %then %do;
	  _src1=&&var&j;
	  %let jump=1;
	 %end;
    run;

    %srcore(srcdata=_take_, srck=&jump, rby=_sample_, srcout=_pevmv&j, pevname=_pevmv&j);  
   %end;
  %end;



  ******************* partielle PEVs *****************************;
  ******************* Leermodell (nur Variablen aus varlist) *****;
  
  %if &varlist ne %then %do;
   data _take_;
   set _bootsr1;
   if _sample_=0;
    %let srckj=0;
    %do j=1 %to &nvar;
     %let jump=0;
     %do jj=1 %to &nclass;
      %if &&class&jj=&&var&j %then %do;
	   %let jump=1;
       %do jjj=1 %to &&df&jj;
	    %let srckj=%eval(&srckj+1);
        _src&srckj=_dum&jj._&jjj;
       %end;
      %end;
	 %end;
	 %if &jump=0 %then %do;
	   %let srckj=%eval(&srckj+1);
 	  _src&srckj=&&var&j;
     %end;
    %end;
   run;


   %srcore(srcdata=_take_, srck=&srckj, rby=_sample_, srcout=_pevfull, pevname=_pevfull);  

  ************************ Partielle Modelle (mit jeweils einem Kandidaten dazu) ****;
   %do k=1 %to &ncand;
    data _take_;
    set _bootsr1;
     %let srckj=0;    
     %do j=1 %to &nvar;
      %let jump=0;
      %do jj=1 %to &nclass;
       %if &&class&jj=&&var&j %then %do; 
	    %let jump=1;
        %do jjj=1 %to &&df&jj;
	     %let srckj=%eval(&srckj+1);
         _src&srckj=_dum&jj._&jjj;
        %end; *do jjj ;
       %end; *if classjj ;

	  %end; * do jj;
	  %if &jump=0 %then %do;
	   %let srckj=%eval(&srckj+1);
 	   _src&srckj=&&var&j;
	  %end; *if jump;
	 %end; *do j;

     %let jump=0;
     %do jj=1 %to &nclass;
      %if &&class&jj=&&cand&k %then %do; 
       %let jump=1;
       %do jjj=1 %to &&df&jj;
        %let srckj=%eval(&srckj+1);
        _src&srckj=_dum&jj._&jjj;
       %end; *do jjj ;
      %end; *if classjj ;
	 %end; * do jj;
	 %if &jump=0 %then %do;
	  %let srckj=%eval(&srckj+1);
 	  _src&srckj=&&cand&k;
     %end; *if jump;
    run;
    %srcore(srcdata=_take_, srck=&srckj, rby=_sample_, srcout=_pevp&k, pevname=_pevp&k);  
   %end;
  %end;



********************                           *****************************;
  data _pev_;
  merge %if &varlist ne %then %do; _pevfull  %end;
   %do j=1 %to &ncand;
    _pevmc&j %if &varlist ne %then %do; _pevp&j %end;
   %end; 
   %if &varlist ne %then %do;
    %do j=1 %to &nvar;
     _pevmv&j
	%end;
   %end; ;
   by _sample_;
   %do j=1 %to &ncand;
    %if &varlist ne %then %do;
	 _asfull=arsin(sqrt(_pevfull));
     if _sample_=0 then do;
      _pevp&j=max(_pevp&j-_pevfull,0);
      _pevp&j._DN=max(_pevp&j._DN-_pevfull_DN,0);
      _pevp&j._DS=max(_pevp&j._DS-_pevfull_DS,0);
	 end;
	 _asp&j=arsin(sqrt(_pevp&j));
	%end;
    _asmc&j=arsin(sqrt(_pevmc&j));
   %end;
   %if &varlist ne %then %do;
    %do j=1 %to &nvar;
	 _asmv&j=arsin(sqrt(_pevmv&j));
	%end;
   %end;
   %do j=1 %to &ncand-1;
    %do k=&j+1 %to &ncand;
     %if &varlist ne %then %do;
      dp&j._&k=_asp&j-_asp&k;
	 %end;
	 dm&j._&k=_asmc&j-_asmc&k;
	%end;
   %end;
  run;

  proc univariate noprint data=_pev_;
  where _sample_>0;
  var %if &varlist ne %then %do;
/*       _asfull   */
	  %end;
      %do j=1 %to &ncand; 
       %if &varlist ne %then %do;
        _asp&j 
	   %end;
       _asmc&j 
      %end;
	  %if &varlist ne %then %do;
	   %do j=1 %to &nvar;
	    _asmv&j
	   %end;
	  %end;
   %do j=1 %to &ncand-1; 
    %do k=&j+1 %to &ncand;
     %if &varlist ne %then %do;
      dp&j._&k 
     %end;
     dm&j._&k
	%end;
   %end; ;
  output out=_relimp_ 

   std=    
   %if &varlist ne %then %do;
/*    sasfull   */
   %end;
   %do j=1 %to &ncand;
    %if &varlist ne %then %do;
	 sasp&j
	%end;
	sasmc&j
   %end;
   %if &varlist ne %then %do;
    %do j=1 %to &nvar;
     sasmv&j 
	%end;
   %end;
   %do j=1 %to &ncand-1; 
    %do k=&j+1 %to &ncand;
	 %if &varlist ne %then %do;
 	  sdp&j._&k 
     %end;
     sdm&j._&k
	%end;
   %end; 
  ;
  run;

 
  data _srcbbm_;
  set _pev_;
  if _sample_=0;
   %if &varlist ne %then %do;
    masfull=_asfull;
   %end;
   %do j=1 %to &ncand;
    %if &varlist ne %then %do;
	 masp&j=_asp&j;
	%end;
	masmc&j=_asmc&j;
   %end;
   %if &varlist ne %then %do;
    %do j=1 %to &nvar;
     masmv&j =_asmv&j;
	%end;
   %end;
   %do j=1 %to &ncand-1; 
    %do k=&j+1 %to &ncand;
	 %if &varlist ne %then %do;
 	  mdp&j._&k=dp&j._&k;
     %end;
     mdm&j._&k=dm&j._&k;
	%end;
   %end; 
   run;

  data _relimp_;
  merge _srcbbm_ _relimp_;
  file print;
   length name $ &varnlen;
   %if &nboot>1 %then %do;
    _chia=0.5*(probit(1-(1-&clp)/2)+sqrt(2*&nboot-3))**2;
	_chi1ma=0.5*(-probit(1-(1-&clp)/2)+sqrt(2*&nboot-3))**2;
    %do j=1 %to &ncand-1;
     %do k=&j+1 %to &ncand;
	  sdm&j._&k.up=sqrt((&nboot-1)*sdm&j._&k**2/_chia);
	  sdm&j._&k.lo=sqrt((&nboot-1)*sdm&j._&k**2/_chi1ma);
      tm&j._&k=abs(mdm&j._&k)/(sdm&j._&k);
      tm&j._&k.lo=abs(mdm&j._&k)/(sdm&j._&k.up);
      tm&j._&k.up=abs(mdm&j._&k)/(sdm&j._&k.lo);
      %if &varlist ne %then %do;
	   sdp&j._&k.up=sqrt((&nboot-1)*sdp&j._&k**2/_chia);
	   sdp&j._&k.lo=sqrt((&nboot-1)*sdp&j._&k**2/_chi1ma);
       tp&j._&k=abs(mdp&j._&k)/(sdp&j._&k);
       tp&j._&k.lo=abs(mdp&j._&k)/(sdp&j._&k.up);
       tp&j._&k.up=abs(mdp&j._&k)/(sdp&j._&k.lo);
      %end;

      %if %eval(&ncand*&correctp)>2 %then %do;
       pm&j._&k=floor((1-probmc('RANGE',tm&j._&k,.,.,&ncand))*10000+0.5)/10000;
	   pm&j._&k.lo=floor((1-probmc('RANGE',tm&j._&k.lo,.,.,&ncand))*10000+0.5)/10000;
	   pm&j._&k.up=floor((1-probmc('RANGE',tm&j._&k,up,.,.,&ncand))*10000+0.5)/10000;
       %if &varlist ne %then %do;
           pp&j._&k=floor((1-probmc('RANGE',tp&j._&k,.,.,&ncand))*10000+0.5)/10000;
           pp&j._&k.lo=floor((1-probmc('RANGE',tp&j._&k.lo,.,.,&ncand))*10000+0.5)/10000;
           pp&j._&k.up=floor((1-probmc('RANGE',tp&j._&k.up,.,.,&ncand))*10000+0.5)/10000;
       %end;
      %end;
      %else %do;
       pm&j._&k=floor((1-probnorm(tm&j._&k))*2*10000+0.5)/10000;
       pm&j._&k.lo=floor((1-probnorm(tm&j._&k.lo))*2*10000+0.5)/10000;
       pm&j._&k.up=floor((1-probnorm(tm&j._&k.up))*2*10000+0.5)/10000;
       %if &varlist ne %then %do;
        pp&j._&k=floor((1-probnorm(tp&j._&k))*2*10000+0.5)/10000;
        pp&j._&k.lo=floor((1-probnorm(tp&j._&k.lo))*2*10000+0.5)/10000;
        pp&j._&k.up=floor((1-probnorm(tp&j._&k.up))*2*10000+0.5)/10000;
       %end;
      %end;
        label %if &varlist ne %then %do;
             pp&j._&k="Partial P-value for &&cand&j &&cand&k"
                %end;
            pm&j._&k="Marginal P-value for &&cand&j &&cand&k";
     %end;
    %end; ;
    hsd=probmc('RANGE',.,0.95,.,3);
   %end;
   %if &varlist ne %then %do;
    mpfull=floor(1000*sin(masfull)**2+0.5)/10;
	label mpfull="PEV for model without candidates";
	%do j=1 %to &nvar;
	 mmv&j=floor(1000*sin(masmv&j)**2+0.5)/10;
	 label mmv&j="Marginal PEV for &&var&j";
	%end;
   %end;
   %do j=1 %to &ncand;
    %if &varlist ne %then %do;
     mpc&j=floor(1000*sin(masp&j)**2+0.5)/10;
	%end;
    mmc&j=floor(1000*sin(masmc&j)**2+0.5)/10;
	label    
     %if &varlist ne %then %do;
      mpc&j="Partial PEV for &&cand&j"
	 %end;
     mmc&j="Marginal PEV for &&cand&j";
   %end;
   run;

   data _relimptab1;
   set _relimp_;
   length _NAME_ $&varnlen _ROLE_ $15;
   %if &varlist ne %then %do;
     %do j=1 %to &nvar;
       _NAME_="&&var&j";
	   _ROLE_="Basic model";
	   EV_marg=mmv&j;
	   DN_marg=round(_pevmv&j._DN*100,0.1);
	   DS_marg=round(_pevmv&j._DS*100,0.1);
	   EV_part=.;
	   DN_part=.;
	   DS_part=.;
	   output;
	 %end;
	 _NAME_="Model"; 
     _ROLE_="w/o candidates";
	 EV_marg=mpfull;
	 DN_marg=round(_pevfull_DN*100,0.1);
	 DS_marg=round(_pevfull_DS*100,0.1);
	 output;
   %end;
   %do j=1 %to &ncand;
     _name_="&&cand&j";
	 _ROLE_="Candidate";
	 EV_marg=mmc&j;
	 DN_marg=round(_pevmc&j._DN*100,0.1);
	 DS_marg=round(_pevmc&j._DS*100,0.1);
	 EV_part=mpc&j;
	 DN_part=round(_pevp&j._DN*100,0.1);
	 DS_part=round(_pevp&j._DS*100,0.1);
	 output;
   %end;
   label _NAME_="Name of variable" _role_="Role" marginal="Marginal PEV" partial="Partial PEV";
   keep _name_ _role_ EV_marg DN_marg DS_marg EV_part DN_part DS_part;
   run;
   title3;
   proc print noobs label;
   title3 "Marginal and Partial EV, DN and DS (%)";
   var _name_ _role_ EV_marg DN_marg DS_marg EV_part DN_part DS_part;
   run;
				   %if &nboot>1 %then %do;
						data _marginal_p;
				   		set _relimp_;
						%do j=1 %to &ncand;
				    	 _name_="&&cand&j";
						 %do jj=&j %to &ncand;
						  &&cand&jj=pm&j._&jj;
						  if pm&j._&jj.lo<&alpha & pm&j._&jj.up>&alpha then _m&jj="*";
						  else _m&jj="";
						  label _m&jj="M";
						 %end;
						 output;
						%end;
						label _name_="Variable";
						keep _name_ %do j=1 %to &ncand; &&cand&j _m&j %end;;
						run;
						proc print noobs label;
						%if &correctp=1 %then %do;
						 title3 "Tukey(HSD)-corrected P-values for comparison of marginal PEV";
						%end;
						%else %do;
						 title3 "P-values for comparison of marginal PEV";
						%end;
						title4 "* denotes: 99% CI for p-value contains the critical value &alpha";
						var _name_ %do j=1 %to &ncand; &&cand&j _m&j %end;;
						run;
						data _partial_p;
				   		set _relimp_;
						%do j=1 %to &ncand;
				    	 _name_="&&cand&j";
						 %do jj=&j %to &ncand;
						  &&cand&jj=pp&j._&jj;
						  if pp&j._&jj.lo<&alpha & pp&j._&jj.up>&alpha then _m&jj="*";
						  else _m&jj="";
						  label _m&jj="M";
						 %end;
						 output;
						%end;
						label _name_="Variable";
						keep _name_ %do j=1 %to &ncand; &&cand&j _m&j %end;;
						run;
						proc print noobs label;
						%if &correctp=1 %then %do;
						 title3 "Tukey(HSD)-corrected P-values for comparison of partial PEV";
						%end;
						%else %do;
						 title3 "P-values for comparison of partial PEV";
						%end;
						title4 "* denotes: 99% CI for p-value contains the critical value &alpha";
						var _name_ %do j=1 %to &ncand; &&cand&j _m&j %end;;
						run;
				   %end;


%end;

****************************************************************************;
********************** total mode ******************************************;
********************** no groups  ******************************************;
****************************************************************************;
%if &ngroup=0 %then %do;
%if &varlist ne %then %do;
 %if &cand = %then %do;

 *****************************  partielle und marginale PEVs, alle mit allen vergleichen ***;
  ******************* marginale PEVs *****************************;
  %do j=1 %to &nvar;
   data _take_;
   set _bootsr1;
    %let jump=0;
    %do jj=1 %to &nclass;
     %if &&class&jj=&&var&j %then %do;
	  %let jump=&&df&jj;
      %do jjj=1 %to &&df&jj;
       _src&jjj=_dum&jj._&jjj;
      %end;
     %end;
	%end;
	%if &jump=0 %then %do;
	 _src1=&&var&j;
	 %let jump=1;
	%end;
   run;

   %srcore(srcdata=_take_, srck=&jump, rby=_sample_, srcout=_pevm&j, pevname=_pevm&j);  
  %end;
  ******************* partielle PEVs *****************************;
  ******************* Vollmodell      *****************************;
   data _take_;
   set _bootsr1;
    %let srckj=0;
    %do j=1 %to &nvar;
	 %let jump=0;
     %do jj=1 %to &nclass;
      %if &&class&jj=&&var&j %then %do;
	   %let jump=1;
       %do jjj=1 %to &&df&jj;
	    %let srckj=%eval(&srckj+1);
        _src&srckj=_dum&jj._&jjj;
       %end;
      %end;
	 %end;
	 %if &jump=0 %then %do;
	   %let srckj=%eval(&srckj+1);
 	  _src&srckj=&&var&j;
     %end;
    %end;

   run;

   %srcore(srcdata=_take_, srck=&srckj, rby=_sample_, srcout=_pevfull, pevname=_pevfull);  

  ************************ Partielle Modelle *************************************;
  %do k=1 %to &nvar;
   data _take_;
   set _bootsr1;
    %let srckj=0;    
    %do j=1 %to &nvar;
	 %if &j ne &k %then %do;
      %let jump=0;
      %do jj=1 %to &nclass;
       %if &&class&jj=&&var&j %then %do;
	    %let jump=1;
        %do jjj=1 %to &&df&jj;
	     %let srckj=%eval(&srckj+1);
         _src&srckj=_dum&jj._&jjj;
        %end;
       %end;
	  %end;
	  %if &jump=0 %then %do;
	   %let srckj=%eval(&srckj+1);
 	   _src&srckj=&&var&j;
	  %end;
	 %end;
	%end;
   run;

   %srcore(srcdata=_take_, srck=&srckj, rby=_sample_, srcout=_pevp&k, pevname=_pevp&k);  

  %end;
  
  data _pev_;
  merge _pevfull 
   %do j=1 %to &nvar;
    _pevm&j _pevp&j
   %end; ;
   by _sample_;
   %do j=1 %to &nvar;
*    _pevp&j=max(_pevfull-_pevp&j,0);
	_asp&j=arsin(sqrt(_pevp&j));
	_asm&j=arsin(sqrt(_pevm&j));
   %end;
   %do j=1 %to &nvar-1;
    %do k=&j+1 %to &nvar;
	 dp&j._&k=_asp&k-_asp&j;
	 dm&j._&k=_asm&j-_asm&k;
	%end;
   %end;
   %do j=1 %to &nvar;
    _pevp&j=max(_pevfull-_pevp&j,0);
    _pevp&j._DN=max(_pevfull_DN-_pevp&j._DN,0);
    _pevp&j._DS=max(_pevfull_DS-_pevp&j._DS,0);
	_asp&j=arsin(sqrt(_pevp&j));
	_asm&j=arsin(sqrt(_pevm&j));
   %end;
  run;


  proc univariate noprint data=_pev_;
  where _sample_>0;
  var %do j=1 %to &nvar; _asp&j _asm&j %end;
   %do j=1 %to &nvar-1; 
    %do k=&j+1 %to &nvar;
	 dp&j._&k dm&j._&k
	%end;
   %end; ;
  output out=_relimp_ 
   std=   %do j=1 %to &nvar; sasp&j sasm&j %end;
   %do j=1 %to &nvar-1; 
    %do k=&j+1 %to &nvar;
	 sdp&j._&k sdm&j._&k
	%end;
   %end; ;
  run;

  data _srcbbm_;
  set _pev_;
  if _sample_=0;
   %do j=1 %to &nvar; 
    masp&j=_asp&j;
    masm&j=_asm&j;
   %end;
   %do j=1 %to &nvar-1; 
    %do k=&j+1 %to &nvar;
	 mdp&j._&k=dp&j._&k;
     mdm&j._&k=dm&j._&k;
	%end;
   %end; 
  run;


  data _relimp_;
  merge _srcbbm_ _relimp_;
  file print;
   length name $ &varnlen;
   %if &nboot>1 %then %do;
    _chia=0.5*(probit(1-(1-&clp)/2)+sqrt(2*&nboot-3))**2;
	_chi1ma=0.5*(-probit(1-(1-&clp)/2)+sqrt(2*&nboot-3))**2;
    %do j=1 %to &nvar-1;
     %do k=&j+1 %to &nvar;
	  sdm&j._&k.up=sqrt((&nboot-1)*sdm&j._&k**2/_chia);
	  sdm&j._&k.lo=sqrt((&nboot-1)*sdm&j._&k**2/_chi1ma);
	  sdp&j._&k.up=sqrt((&nboot-1)*sdp&j._&k**2/_chia);
	  sdp&j._&k.lo=sqrt((&nboot-1)*sdp&j._&k**2/_chi1ma);

      tm&j._&k=abs(mdm&j._&k)/(sdm&j._&k*(&nboot-1)/(&nboot));
      tp&j._&k=abs(mdp&j._&k)/(sdp&j._&k*(&nboot-1)/(&nboot));
      tm&j._&k.lo=abs(mdm&j._&k)/(sdm&j._&k.lo*(&nboot-1)/(&nboot));
      tp&j._&k.lo=abs(mdp&j._&k)/(sdp&j._&k.lo*(&nboot-1)/(&nboot));
      tm&j._&k.up=abs(mdm&j._&k)/(sdm&j._&k.up*(&nboot-1)/(&nboot));
      tp&j._&k.up=abs(mdp&j._&k)/(sdp&j._&k.up*(&nboot-1)/(&nboot));
      %if %eval(&nvar*&correctp)>3 %then %do;
       pm&j._&k=floor((1-probmc('RANGE',tm&j._&k,.,.,&nvar))*10000+0.5)/10000;
       pp&j._&k=floor((1-probmc('RANGE',tp&j._&k,.,.,&nvar))*10000+0.5)/10000;
       pm&j._&k.up=floor((1-probmc('RANGE',tm&j._&k.lo,.,.,&nvar))*10000+0.5)/10000;
       pp&j._&k.up=floor((1-probmc('RANGE',tp&j._&k.lo,.,.,&nvar))*10000+0.5)/10000;
       pm&j._&k.lo=floor((1-probmc('RANGE',tm&j._&k.up,.,.,&nvar))*10000+0.5)/10000;
       pp&j._&k.lo=floor((1-probmc('RANGE',tp&j._&k.up,.,.,&nvar))*10000+0.5)/10000;
      %end;
      %else %do;
       pm&j._&k=floor((1-probnorm(tm&j._&k))*2*10000+0.5)/10000;
       pp&j._&k=floor((1-probnorm(tp&j._&k))*2*10000+0.5)/10000;
       pm&j._&k.up=floor((1-probnorm(tm&j._&k.lo))*2*10000+0.5)/10000;
       pp&j._&k.up=floor((1-probnorm(tp&j._&k.lo))*2*10000+0.5)/10000;
       pm&j._&k.lo=floor((1-probnorm(tm&j._&k.up))*2*10000+0.5)/10000;
       pp&j._&k.lo=floor((1-probnorm(tp&j._&k.up))*2*10000+0.5)/10000;
      %end;
      label pp&j._&k="Partial P-value for &&var&j &&var&k"
            pm&j._&k="Marginal P-value for &&var&j &&var&k";
     %end;
    %end; ;
    hsd=probmc('RANGE',.,0.95,.,3);
   %end;
   %do j=1 %to &nvar;
    mp&j=floor(1000*sin(masp&j)**2+0.5)/10;
    mm&j=floor(1000*sin(masm&j)**2+0.5)/10;
	label mp&j="Partial PEV for &&var&j"
          mm&j="Marginal PEV for &&var&j";
   %end;
   run;

   data _relimptab1;
   set _relimp_;
   length _NAME_ $&varnlen ;
   %if &varlist ne %then %do;
     %do j=1 %to &nvar;
       _NAME_="&&var&j";
	   EV_marg=mm&j;
	   DN_marg=round(_pevm&j._DN*100,0.1);
	   DS_marg=round(_pevm&j._DS*100,0.1);
	   EV_part=mp&j;
	   DN_part=round(_pevp&j._DN*100,0.1);
	   DS_part=round(_pevp&j._DS*100,0.1);
	   output;
	 %end;
	 _NAME_="Model";
	 EV_marg=round(_pevfull*100,0.1);
	 DN_marg=round(_pevfull_DN*100,0.1);
	 DS_marg=round(_pevfull_DS*100,0.1);
	 EV_part=.;
	 DN_part=.;
	 DS_part=.;
	 output;
   %end;
   label _NAME_="Name of variable"  marginal="Marginal PEV" partial="Partial PEV";
   keep _name_ EV_marg DN_marg DS_marg EV_part DN_part DS_part;
   run;
   title3;
   proc print noobs label;
   title3 "Marginal and Partial EV, DN and DS (%)";
   var _name_ EV_marg DN_marg DS_marg EV_part DN_part DS_part;
   run;
			   %if &nboot>1 %then %do;
					data _marginal_p;
			   		set _relimp_;
					%do j=1 %to &nvar;
			    	 _name_="&&var&j";
					 %do jj=&j %to &nvar;
					  &&var&jj=pm&j._&jj;
					  if pm&j._&jj.lo<&alpha & pm&j._&jj.up>&alpha then _m&jj="*";
					  else _m&jj="";
					  label _m&jj="M";
					 %end;
					 output;
					%end;
					label _name_="Variable";
					keep _name_ %do j=1 %to &nvar; &&var&j _m&j %end;; 
					run;
					proc print noobs label;
					%if &correctp=1 %then %do;
					 title3 "Tukey(HSD)-corrected P-values for comparison of marginal PEV";
					%end;
					%else %do;
					 title3 "P-values for comparison of marginal PEV";
					%end;
					title4 "* denotes: 99% CI for p-value contains the critical value &alpha";
					var _name_ %do j=1 %to &nvar; &&var&j _m&j %end;;
					run;
					data _partial_p;
			   		set _relimp_;
					%do j=1 %to &nvar;
			    	 _name_="&&var&j";
					 %do jj=&j %to &nvar;
					  &&var&jj=pp&j._&jj;
					  if pp&j._&jj.lo<&alpha & pp&j._&jj.up>&alpha then _m&jj="*";
					  else _m&jj="";
					  label _m&jj="M";
					 %end;
					 output;
					%end;
					label _name_="Variable";
					keep _name_ %do j=1 %to &nvar; &&var&j _m&j %end;; 
					run;
					proc print noobs label;
					%if &correctp=1 %then %do;
					 title3 "Tukey(HSD)-corrected P-values for comparison of partial PEV";
					%end;
					%else %do;
					 title3 "P-values for comparison of partial PEV";
					%end;
					title4 "* denotes: 99% CI for p-value contains the critical value &alpha";
					var _name_ %do j=1 %to &nvar; &&var&j _m&j %end;;
					run;
			   %end;

%end;
%end;
%end;

****************************************************************************;
********************** total mode ******************************************;
********************** groups     ******************************************;
****************************************************************************;

%if &ngroup>0 %then %do;
 %if &varlist ne %then %do;
  %if &cand = %then %do;

  *****************************  partielle und marginale PEVs, alle Gruppen mit allen vergleichen ***;
  ******************* marginale PEVs *****************************;
  %do g=1 %to &ngroup;
   data _take_;
   set _bootsr1;
    %let jump=0;
    %do j=1 %to &nvar;
	 %if &&group&j=&g %then %do;
	  %let classok=0;
      %do jj=1 %to &nclass;
       %if &&class&jj=&&var&j %then %do; 
	    %let classok=1;
        %do jjj=1 %to &&df&jj;
		 %let jump=%eval(&jump+1);
         _src&jump=_dum&jj._&jjj;
        %end; *jjj;
       %end; *if;
	  %end; *jj;
	  %if &classok=0 %then %do;
	   %let jump=%eval(&jump+1);
	   _src&jump=&&var&j;
	  %end; *if classok;
	 %end; *if groupj=g;
	%end; *j;
   run;

   %srcore(srcdata=_take_, srck=&jump, rby=_sample_, srcout=_pevm&g, pevname=_pevm&g);  

  %end;
  ******************* partielle PEVs *****************************;
  ******************* Vollmodell      *****************************;
   data _take_;
   set _bootsr1;
    %let srckj=0;
    %do j=1 %to &nvar;
	 %let jump=0;
     %do jj=1 %to &nclass;
      %if &&class&jj=&&var&j %then %do;
	   %let jump=1;
       %do jjj=1 %to &&df&jj;
	    %let srckj=%eval(&srckj+1);
        _src&srckj=_dum&jj._&jjj;
       %end;
      %end;
	 %end;
	 %if &jump=0 %then %do;
	   %let srckj=%eval(&srckj+1);
 	  _src&srckj=&&var&j;
     %end;
    %end;

   run;

   %srcore(srcdata=_take_, srck=&srckj, rby=_sample_, srcout=_pevfull, pevname=_pevfull);  

  ************************ Partielle Modelle *************************************;
  %do k=1 %to &ngroup;
   data _take_;
   set _bootsr1;
    %let srckj=0;    
    %do j=1 %to &nvar;
	 %if &&group&j ne &k %then %do;
	  %let classok=0;
      %do jj=1 %to &nclass;
       %if &&class&jj=&&var&j %then %do;
	    %let classok=1;
        %do jjj=1 %to &&df&jj;
	     %let srckj=%eval(&srckj+1);
         _src&srckj=_dum&jj._&jjj;
        %end;
       %end;
	  %end;
	  %if &classok=0 %then %do;
      %let srckj=%eval(&srckj+1);
       _src&srckj=&&var&j;
	  %end;
	 %end;
    %end;
   run;

   %srcore(srcdata=_take_, srck=&srckj, rby=_sample_, srcout=_pevp&k, pevname=_pevp&k);  

  %end;
  
  data _pev_;
  merge _pevfull 
   %do j=1 %to &ngroup;
    _pevm&j _pevp&j
   %end; ;
   by _sample_;
   %do j=1 %to &ngroup;
*    _pevp&j=max(_pevfull-_pevp&j,0);
	_asp&j=arsin(sqrt(_pevp&j));
	_asm&j=arsin(sqrt(_pevm&j));
   %end;
   %do j=1 %to &ngroup-1;
    %do k=&j+1 %to &ngroup;
	 dp&j._&k=_asp&k-_asp&j;
	 dm&j._&k=_asm&j-_asm&k;
	%end;
   %end;
   %do j=1 %to &ngroup;
    _pevp&j=max(_pevfull-_pevp&j,0);
    _pevp&j._DN=max(_pevfull_DN-_pevp&j._DN,0);
    _pevp&j._DS=max(_pevfull_DS-_pevp&j._DS,0);
	_asp&j=arsin(sqrt(_pevp&j));
	_asm&j=arsin(sqrt(_pevm&j));
   %end;
  run;

/*proc print data=_pev_ noobs;
	var _pevm:;
	run;*/

  proc univariate noprint data=_pev_;
  where _sample_>0;
  var %do j=1 %to &ngroup; _asp&j _asm&j %end;
   %do j=1 %to &ngroup-1; 
    %do k=&j+1 %to &ngroup;
	 dp&j._&k dm&j._&k
	%end;
   %end; ;
  output out=_relimp_ 
   std=   %do j=1 %to &ngroup; sasp&j sasm&j %end;
   %do j=1 %to &ngroup-1; 
    %do k=&j+1 %to &ngroup;
	 sdp&j._&k sdm&j._&k
	%end;
   %end; ;
  run;

  data _srcbbm_;
  set _pev_;
  if _sample_=0;
   %do j=1 %to &ngroup; 
    masp&j=_asp&j;
    masm&j=_asm&j;
   %end;
   %do j=1 %to &ngroup-1; 
    %do k=&j+1 %to &ngroup;
	 mdp&j._&k=dp&j._&k;
     mdm&j._&k=dm&j._&k;
	%end;
   %end; 
  run;


  data _relimp_;
  merge _srcbbm_ _relimp_;
			   %if &nboot>1 %then %do;
			    _chia=0.5*(probit(1-(1-&clp)/2)+sqrt(2*&nboot-3))**2;
				_chi1ma=0.5*(-probit(1-(1-&clp)/2)+sqrt(2*&nboot-3))**2;
			    %do j=1 %to &ngroup-1;
			     %do k=&j+1 %to &ngroup;
				  sdm&j._&k.up=sqrt((&nboot-1)*sdm&j._&k**2/_chia);
				  sdm&j._&k.lo=sqrt((&nboot-1)*sdm&j._&k**2/_chi1ma);
				  sdp&j._&k.up=sqrt((&nboot-1)*sdp&j._&k**2/_chia);
				  sdp&j._&k.lo=sqrt((&nboot-1)*sdp&j._&k**2/_chi1ma);

			      tm&j._&k=abs(mdm&j._&k)/(sdm&j._&k*(&nboot-1)/(&nboot));
			      tp&j._&k=abs(mdp&j._&k)/(sdp&j._&k*(&nboot-1)/(&nboot));
			      tm&j._&k.lo=abs(mdm&j._&k)/(sdm&j._&k.lo*(&nboot-1)/(&nboot));
			      tp&j._&k.lo=abs(mdp&j._&k)/(sdp&j._&k.lo*(&nboot-1)/(&nboot));
			      tm&j._&k.up=abs(mdm&j._&k)/(sdm&j._&k.up*(&nboot-1)/(&nboot));
			      tp&j._&k.up=abs(mdp&j._&k)/(sdp&j._&k.up*(&nboot-1)/(&nboot));
			      %if %eval(&ngroup*&correctp)>2 %then %do;
			       pm&j._&k=floor((1-probmc('RANGE',tm&j._&k,.,.,&ngroup))*10000+0.5)/10000;
			       pp&j._&k=floor((1-probmc('RANGE',tp&j._&k,.,.,&ngroup))*10000+0.5)/10000;
			       pm&j._&k.up=floor((1-probmc('RANGE',tm&j._&k.lo,.,.,&ngroup))*10000+0.5)/10000;
			       pp&j._&k.up=floor((1-probmc('RANGE',tp&j._&k.lo,.,.,&ngroup))*10000+0.5)/10000;
			       pm&j._&k.lo=floor((1-probmc('RANGE',tm&j._&k.up,.,.,&ngroup))*10000+0.5)/10000;
			       pp&j._&k.lo=floor((1-probmc('RANGE',tp&j._&k.up,.,.,&ngroup))*10000+0.5)/10000;
			      %end;
			      %else %do;
			       pm&j._&k=floor((1-probnorm(tm&j._&k))*2*10000+0.5)/10000;
			       pp&j._&k=floor((1-probnorm(tp&j._&k))*2*10000+0.5)/10000;
			       pm&j._&k.up=floor((1-probnorm(tm&j._&k.lo))*2*10000+0.5)/10000;
			       pp&j._&k.up=floor((1-probnorm(tp&j._&k.lo))*2*10000+0.5)/10000;
			       pm&j._&k.lo=floor((1-probnorm(tm&j._&k.up))*2*10000+0.5)/10000;
			       pp&j._&k.lo=floor((1-probnorm(tp&j._&k.up))*2*10000+0.5)/10000;
			      %end;
			      label pp&j._&k="Partial P-value for &&gname&j &&gname&k"
			            pm&j._&k="Marginal P-value for &&gname&j &&gname&k";
			     %end;
			    %end; ;
			    hsd=probmc('RANGE',.,0.95,.,3);
			   %end;
   %do j=1 %to &ngroup;
    mp&j=floor(1000*sin(masp&j)**2+0.5)/10;
    mm&j=floor(1000*sin(masm&j)**2+0.5)/10;
	label mp&j="Partial PEV for &&gname&j"
          mm&j="Marginal PEV for &&gname&j";
   %end;
   run;


  %let maxlen=1;
  %let maxname=1;
  %do j=1 %to &ngroup;
    %let vars&j=;
	%if %length(&&gname&j)>&maxname %then %let maxname=%length(&&gname&j);
    %do k=1 %to &nvar; 
     %if &&group&k=&j %then %do; 
      %let vars&j=&&vars&j &&var&k ; 
	  %if %length(&&vars&j)>&maxlen %then %let maxlen=%length(&&vars&j);
     %end; 
    %end;
   %end;


data _relimptab1;
   set _relimp_;
   length _NAME_ $&maxname _Vars_ $&maxlen;
     %do j=1 %to &ngroup;
       _NAME_="&&gname&j";
	   _VARS_="&&vars&j";
	   EV_marg=mm&j;
	   DN_marg=round(_pevm&j._DN*100,0.1);
	   DS_marg=round(_pevm&j._DS*100,0.1);
	   EV_part=mp&j;
	   DN_part=round(_pevp&j._DN*100,0.1);
	   DS_part=round(_pevp&j._DS*100,0.1);
	   output;
	 %end;
	 _NAME_="Model"; 
	 _vars_="All";
	 EV_marg=round(_pevfull*100,0.1);
	 DN_marg=round(_pevfull_DN*100,0.1);
	 DS_marg=round(_pevfull_DS*100,0.1);
	 EV_part=.;
	 DN_part=.;
	 DS_part=.;
	 output;
   label _NAME_="Name of group" _vars_="Variables";
   keep _name_ _vars_ EV_marg DN_marg DS_marg EV_part DN_part DS_part;
   run;
   title3;
   proc print noobs label;
   title3 "Marginal and Partial EV, DN and DS (%)";
   var _name_ _vars_ EV_marg DN_marg DS_marg EV_part DN_part DS_part;
   run;
			   %if &nboot>1 %then %do;
					data _marginal_p;
			   		set _relimp_;
					length _name_ $ &maxname;
					%do j=1 %to &ngroup;
			    	 _name_="&&gname&j";
					 %do jj=&j %to &ngroup;
					  %qcmpres(&&gname&jj)=pm&j._&jj;
					  if pm&j._&jj.lo<&alpha & pm&j._&jj.up>&alpha then _m&jj="*";
					  else _m&jj="";
					  label _m&jj="M";
					 %end;
					 output;
					%end;
					label _name_="Group";
					keep _name_ %do j=1 %to &ngroup; &&gname&j _m&j %end;;
					run;
					proc print noobs label;
					%if &correctp=1 %then %do;
					 title3 "Tukey(HSD)-corrected P-values for comparison of marginal PEV";
					%end;
					%else %do;
					 title3 "P-values for comparison of marginal PEV";
					%end;
					title4 "* denotes: 99% CI for p-value contains the critical value &alpha";
					var _name_ %do j=1 %to &ngroup; &&gname&j _m&j %end;;
					run;
					data _partial_p;
			   		set _relimp_;
					length _name_ $ &maxname;
					%do j=1 %to &ngroup;
			    	 _name_="&&gname&j";
					 %do jj=&j %to &ngroup;
					  %qcmpres(&&gname&jj)=pp&j._&jj;
					  if pp&j._&jj.lo<&alpha & pp&j._&jj.up>&alpha then _m&jj="*";
					  else _m&jj="";
					  label _m&jj="M";
					 %end;
					 output;
					%end;
					label _name_="Group";
					keep _name_ %do j=1 %to &ngroup; &&gname&j _m&j %end;
					run;
					proc print noobs label;
					%if &correctp=1 %then %do;
					 title3 "Tukey(HSD)-corrected P-values for comparison of partial PEV";
					%end;
					%else %do;
					 title3 "P-values for comparison of partial PEV";
					%end;
					title4 "* denotes: 99% CI for p-value contains the critical value &alpha";
					var _name_ %do j=1 %to &ngroup; &&gname&j _m&j %end;;
					run;
			   %end;

 %end;
%end;
%end;

%if &nboot*&printclp>1 %then %do;
 %let clpp=%sysevalf(&clp*100);
 %let clpp=&clpp.%;
 data _CLPtable_;
 set _relimp_;
 %if &ngroup>0 %then %do;
  length na1 $ &maxname na2 $ &maxname;
 %end;
 %else %do;
  length na1 $ &varnlen na2 $ &varnlen;
 %end;
 %if &ngroup>0 %then %let bis=&ngroup;
 %else %if &ncand>1 %then %let bis=&ncand;
 %else %if &varlist ne %then %let bis=&nvar;
 %do j=1 %to %eval(&bis-1);
  %do k=%eval(1+&j) %to &bis;
   %if &ngroup>0 %then %do; 
    na1="&&gname&j";
	na2="&&gname&k";
   %end;
   %else %if &ncand>1 %then %do;
    na1="&&cand&j";
	na2="&&cand&k";
   %end;
   %else %if &varlist ne %then %do;
    na1="&&var&j";
	na2="&&var&k";
   %end;
   pm_lo=pm&&j._&k.lo;
   pm_up=pm&&j._&k.up;
   pp_lo=pp&&j._&k.lo;
   pp_up=pp&&j._&k.up;
   output;
  %end;
 %end;
 label na1="Term" na2="with Term" pm_lo="Marginal PEV p-value lower &clpp CL" pm_up="Marginal PEV p-value upper &clpp CL" pp_lo="Partial PEV p-value lower &clpp CL" pp_up="Partial PEV p-value upper &clpp CL";
run;
proc print noobs label;
title3 "&clpp confidence intervals for p-values for comparison of PEV";
title4 "based on &nboot bootstrap resamples";
var na1 na2 pm_lo pm_up pp_lo pp_up;
run;

%end;


%if &histogram=1 %then %do;

 proc capability  noprint data=_pev_;
 var
 %if &ngroup>0 %then %let bis=&ngroup;
 %else %if &ncand>1 %then %let bis=&ncand;
 %else %if &varlist ne %then %let bis=&nvar;
 %do j=1 %to %eval(&bis-1);
  %do k=%eval(1+&j) %to &bis;
   dm&j._&k dp&j._&k
  %end;
 %end;
 ;
 histogram
 %if &ngroup>0 %then %let bis=&ngroup;
 %else %if &ncand>1 %then %let bis=&ncand;
 %else %if &varlist ne %then %let bis=&nvar;
 %do j=1 %to %eval(&bis-1);
  %do k=%eval(1+&j) %to &bis;
   dm&j._&k dp&j._&k
  %end;
 %end;
 ;
 run;

%end;


data _xxx;
set _xxx;
time=time()-time;
time=floor(time*100+0.5)/100;
minutes=floor(time/60);
seconds=time-minutes*60;
*if minutes>0 then put "NOTE: Execution of macro used " minutes "minutes " seconds "seconds.";
*else put "NOTE: Execution of macro used " seconds "seconds.";
call symput('minute',compress(minutes));
call symput('second',compress(seconds));
run;

options nonotes nosource;
proc datasets noprint;
delete _bootsr _bootsr1 _ndat _nkm _pevfull _srcb _srcbbm_ _srcs _srct
       _srcu _srcu1 _srcu2 _srcu3 _srcu4 _srcu5 _srcu6 _srcv _take_ _ttt_ _worksr _xxx
	   %if &cand ne %then %do; 
        %do c=1 %to &ncand; 
         _pevmc&c _pevp&c
        %end; 
	    %if &varlist ne %then %do;
		 %do v=1 %to &nvar;
		  _pevmv&v
		 %end;
		%end;
	   %end;
	   %else %do;
	    %do v=1 %to &nvar;
	     _pevm&v _pevp&v
    	%end;
	   %end;
	   %if &class ne %then %do;
	    %do c=1 %to &nclass;
		 _clt&c
		%end;
		_df _dummies
	   %end;
	   ;
run;
quit;

options notes source;
title; title2; title3; title4;
%put NOTE: Execution of macro used &minute minutes &second seconds.;
%mend;


/*
%NecSuff_Surv(data=sh_tab3.pbcneu, time=zeit, cens=status, 
			varlist=age edema logbil albu logpro,
			class=edema);
%NecSuff_Surv(data=sh_tab3.pbcneu, time=zeit, cens=status, 
			varlist=age  albu edema logpro,
			class=edema, cand=logbil);
%NecSuff_Surv(data=sh_tab3.pbcneu, time=zeit, cens=status, 
			varlist=age edema logbil albu logpro,
			class=edema, nboot=100);
%NecSuff_Surv(data=sh_tab3.pbcneu, time=zeit, cens=status, 
			varlist=age edema logbil albu logpro,
			class=edema, groups=1 2 3 3 3, gnames=age edema labor);
*/
