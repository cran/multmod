
// includes from the plugin

#include <Rcpp.h>


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes

using namespace Rcpp;
double pnorm(const double x)
{
  const double c0 = 0.2316419;
  const double c1 = 1.330274429;
  const double c2 = 1.821255978;
  const double c3 = 1.781477937;
  const double c4 = 0.356563782;
  const double c5 = 0.319381530;
  const double c6 = 0.398942280401;
  const double negative = (x < 0 ? 1.0 : 0.0);
  const double xPos = (x < 0.0 ? -x : x);
  const double k = 1.0 / ( 1.0 + (c0 * xPos));
  const double y1 = (((((((c1*k-c2)*k)+c3)*k)-c4)*k)+c5)*k;
  const double y2 = 1.0 - (c6*std::exp(-0.5*xPos*xPos)*y1);
  return ((1.0-negative)*y2) + (negative*(1.0-y2));
}


using namespace Rcpp;
double bvndfun(const double dh, const double dk, const double r)
{
double twopi=2*3.14159265358979323846;

double Wm [10][3];
double Xm [10][3];

Wm [0][0] = 0.1713244923791705;
Wm [1][0] = 0.3607615730481384;
Wm [2][0] = 0.4679139345726904;

Xm [0][0] = -0.9324695142031522;
Xm [1][0] = -0.6612093864662647;
Xm [2][0] = -0.2386191860831970;

Wm [0][1] = 0.4717533638651177 / 10;
Wm [1][1] = 0.1069393259953183;
Wm [2][1] = 0.1600783285433464;
Wm [3][1] = 0.2031674267230659;
Wm [4][1] = 0.2334925365383547;
Wm [5][1] = 0.2491470458134029;

Xm [0][1] = -0.9815606342467191;
Xm [1][1] = -0.9041172563704750;
Xm [2][1] = -0.7699026741943050;
Xm [3][1] = -0.5873179542866171;
Xm [4][1] = -0.3678314989981802;
Xm [5][1] = -0.1252334085114692;

Wm [0][2] = 0.1761400713915212 / 10;
Wm [1][2] = 0.4060142980038694 / 10;
Wm [2][2] = 0.6267204833410906 / 10;
Wm [3][2] = 0.8327674157670475 / 10;
Wm [4][2] = 0.1019301198172404;
Wm [5][2] = 0.1181945319615184;
Wm [6][2] = 0.1316886384491766;
Wm [7][2] = 0.1420961093183821;
Wm [8][2] = 0.1491729864726037;
Wm [9][2] = 0.1527533871307259;

Xm [0][2] = -0.9931285991850949;
Xm [1][2] = -0.9639719272779138;
Xm [2][2] = -0.9122344282513259;
Xm [3][2] = -0.8391169718222188;
Xm [4][2] = -0.7463319064601508;
Xm [5][2] = -0.6360536807265150;
Xm [6][2] = -0.5108670019508271;
Xm [7][2] = -0.3737060887154196;
Xm [8][2] = -0.2277858511416451;
Xm [9][2] = -0.7652652113349733;

int ng;
int lg;
if (fabs(r) < 0.3)
    {
        ng = 1;
        lg = 3;
    } else {
        if (fabs(r) < 0.75)
        {
            ng = 2;
            lg = 6;
        } else {
            ng = 3;
            lg = 10;
        } }       
   
double H = dh;
double K = dk;
double HK = H*K;
    

double bvn=0;
    if (fabs(r) < 0.925)
    {
        if (fabs(r) > 0)
        {
            double HS = (H*H + K*K) / 2;
            double ASR = asin(r);
            for (int i=0; i < lg; i++)
            {
                for (int j=-1; j<2; j+=2)
                {
                    double SN = sin(ASR*(j*Xm[i] [ng-1] + 1)/2);
                    bvn = bvn + Wm[i][ng-1] * exp((SN*HK-HS)/(1-SN*SN));
                }
            }
            bvn = bvn*ASR/(2*twopi);
        }
        bvn = bvn + pnorm(-H) * pnorm(-K);
    }
 else {
        if (r < 0)
        {
            K = -K;
            HK = -HK;
        }
        if (fabs(r) < 1)
        {
            double AS = (1-r)*(1+r);
            double A = sqrt(AS);
            double BS = (H-K) * (H-K);
            double CC = (4-HK) / 8;
            double DD = (12-HK) / 16;
            double ASR = -(BS/AS + HK) / 2;
            if (ASR > -100)
            {
                bvn = A*exp(ASR)*(1 - CC*(BS-AS)*(1-DD*BS/5)/3 + CC*DD*AS*AS/5);
            }
            if (-HK < 100)
            {
                double B = sqrt(BS);
                double bvn = bvn - exp( -HK/2)*sqrt(twopi)*pnorm(-B/A)*B*(1 - CC*BS*(1 - DD*BS/5)/3);
            }
            A = A/2;
            for (int i=0; i < lg; i++)
            {
                for (int j=-1; j<2; j+=2)
                {
                    double XS = (A*( j*Xm[i][ng-1] + 1))*(A*( j*Xm[i][ng-1] + 1));
                    double RS = sqrt(1-XS);
                    double ASR = -(BS/XS+HK)/2;
                    if (ASR > -100)
                    {
                        bvn = bvn + A*Wm[i][ng-1]*exp(ASR)*(exp(-HK*(1-RS)/(2*(1+RS)))/RS - (1+CC*XS*(1+DD*XS)));
                    }
                }
            }
            bvn = -bvn/twopi;
        }
        if (r > 0)
        {
  double maxHK=H;
	if(K>H) maxHK=K; 
            bvn = bvn + pnorm(-maxHK);
        } else {
            bvn = -bvn;
            if (K > H)
            {
                bvn = bvn + pnorm(K) - pnorm(H);
            }
        }
    }

return bvn;
}


using namespace Rcpp;
double probfun(const double xc, const double xcor)
{
return (pnorm(xc)-pnorm(-xc)-(bvndfun(xc,xc,xcor)-bvndfun(-xc,xc,xcor)-bvndfun(xc,-xc,xcor)+bvndfun(-xc,-xc,xcor)));
}


// declarations
extern "C" {
SEXP file5d5320db( SEXP iidresp, SEXP numobs, SEXP nummod, SEXP zvec) ;
SEXP file640c223a( SEXP iidresp, SEXP numobs, SEXP zvec) ;	
}



// definition


SEXP file640c223a( SEXP iidresp, SEXP numobs, SEXP zvec ){
	BEGIN_RCPP
	
	Rcpp::NumericMatrix Am(iidresp); 
	int Xnumobs = as<int>(numobs);
	Rcpp::NumericVector xzvec(zvec);
	int nrows = Am.nrow(); 
	int ncolumns = Am.ncol(); 
	Rcpp::NumericVector Xvardiag(nrows); 
	
	for (int i = 0; i < nrows; i++) { 
		for (int k = 0; k < ncolumns; k++) {
			Xvardiag(i) += Am(i,k)*Am(i,k)/Xnumobs;       
		}
	}
	
	for (int i = 0; i < nrows; i++){
		xzvec(i) *= 1/sqrt(Xvardiag(i));
	}
	
	return xzvec;
	
	
	
	
	END_RCPP
}


SEXP file5d5320db( SEXP iidresp, SEXP numobs, SEXP nummod, SEXP zvec ){
BEGIN_RCPP

Rcpp::NumericMatrix Am(iidresp); 
int Xnumobs = as<int>(numobs);
int xnummod = as<int>(nummod);
Rcpp::NumericVector xzvec(zvec); 
int nrows = Am.nrow(); 
int ncolumns = Am.ncol(); 
Rcpp::NumericVector Xvardiag(nrows); 
Rcpp::NumericVector ycorvec(nrows*nrows); 
int n_ycorvec=ycorvec.size();
Rcpp::NumericVector ycorlist(xnummod-1);


for (int i = 0; i < nrows; i++) { 
for (int k = 0; k < ncolumns; k++) {
Xvardiag(i) += Am(i,k)*Am(i,k)/Xnumobs;       
}
}
for (int i = 0; i < nrows; i++) { 
for (int j = 0; j < nrows; j++) {         
for (int k = 0; k < ncolumns; k++) {
ycorvec(i*nrows + j) += (Am(i,k)*Am(j,k)/Xnumobs)/(sqrt(Xvardiag(i))*sqrt(Xvardiag(j)));  	
						// B vektoren = 1. rÊkke i correlationsmatricen, 
				     		// 2. rÊkke, 3. rÊkke osv. 
}
} 
} 

ycorvec=abs(ycorvec);

// IndsÊtter -10 i diagonalen sÂ der ikke vÊlges et element her

for(int i = 0; i<xnummod; i++)
ycorvec(i+i*xnummod)=-10;

// Finder f¯rste max og noterer vÊrdi og placering 

double maxVal = ycorvec(0);   	// Initialisering af maxVal
int maxPos=-1;			// Initialisering af maxPos
for (int i=1; i<n_ycorvec; i++) {
        if (ycorvec(i) > maxVal) {
            maxVal = ycorvec(i);
	    maxPos = i;
        }
    }

// F¯rste max-vÊrdi indsÊttes pÂ f¯rste plads i listen

ycorlist(0)=maxVal;

// RÊkke og s¯jle nummer for placering af max findes og der 
// indsÊttes -10 pÂ de to pladser (symmetrisk matrix) 
// sÂ den fundne vÊrdi ikke vil blive brugt igen

int pos [2];
pos[0] = floor((maxPos)/xnummod);
pos[1] =(maxPos % xnummod);

ycorvec(pos[0]*xnummod+pos[1])=-10;
ycorvec(pos[1]*xnummod+pos[0])=-10;



// Finder nÊste max i rÊkker der svarer til rÊkker fra tidligere max

for(int k=1; k<(xnummod-1);k++){

maxVal = ycorvec[pos[0]*xnummod];	// Initialisering af maxVal til f¯rste element 
					// i den ene af de to rÊkker der kan s¯ges i
maxPos = pos[0]*xnummod;		// Initialisering af maxPos til f¯rste position der ledes i

for(int j=0; j<2; j++ ){
for (int i=0; i<xnummod; i++) {
        if (ycorvec(pos[j]*xnummod+i) > maxVal) {
            maxVal = ycorvec(pos[j]*xnummod+i);
	    maxPos = pos[j]*xnummod+i;    //
        }
    }
}

ycorlist(k)=maxVal;

int pos2 [2];

pos2[0] = floor((maxPos)/xnummod);
pos2[1] =(maxPos % xnummod);

ycorvec(pos2[0]*xnummod+pos2[1])=-10;
ycorvec(pos2[1]*xnummod+pos2[0])=-10;


if(pos[0]==pos2[0]){
for(int i=0; i<xnummod; i++){
ycorvec(pos[0]*xnummod+i)=-10;
ycorvec(i*xnummod+pos[0])=-10;
}
ycorvec(pos[1]*xnummod+pos2[1])=-10;
ycorvec(pos2[1]*xnummod+pos[1])=-10;
pos[0]=pos[1];
pos[1]=pos2[1];
}

else if(pos[1]==pos2[1]){
for(int i=0; i<xnummod; i++){
ycorvec(pos[1]*xnummod+i)=-10;
ycorvec(i*xnummod+pos[1])=-10;
}
ycorvec(pos[0]*xnummod+pos2[0])=-10;
ycorvec(pos2[0]*xnummod+pos[0])=-10;
pos[0]=pos[0];
pos[1]=pos2[0];
}

else if(pos[0]==pos2[1]){
for(int i=0; i<xnummod; i++){
ycorvec(pos[0]*xnummod+i)=-10;
ycorvec(i*xnummod+pos[0])=-10;
}
ycorvec(pos[1]*xnummod+pos2[0])=-10;
ycorvec(pos2[0]*xnummod+pos[1])=-10;
pos[0]=pos[1];
pos[1]=pos2[0];
}

else if(pos[1]==pos2[0]){
for(int i=0; i<xnummod; i++){
ycorvec(pos[1]*xnummod+i)=-10;
ycorvec(i*xnummod+pos[1])=-10;}

ycorvec(pos[0]*xnummod+pos2[1])=-10;
ycorvec(pos2[1]*xnummod+pos[0])=-10;
pos[0]=pos[0];
pos[1]=pos2[1];
}


}


Rcpp::NumericVector xsum(xnummod);

for(int i=0; i<xnummod; i++ ){
xsum(i)=0;
for(int j=0; j<(xnummod-1); j++){ 
xsum(i)+= probfun(xzvec(i),ycorlist(j));
}
}

return xsum;


END_RCPP
}



