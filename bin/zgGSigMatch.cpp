// mex -v CXXFLAGS='-ffast-math -O3' SPMatch4c.cpp

#include <mex.h>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <climits>
#include <vector>
#include <algorithm>
#include <functional>
using namespace std;

// #if 0
//     #define MAX_CRITERION
// #endif

const unsigned idxSamples       = 0;
const unsigned idxSamLocations	= 1;
const unsigned idxSamScales 	= 2;
const unsigned idxExemplars     = 3;
const unsigned idxExeLocations	= 4;
const unsigned idxExeScales 	= 5;
const unsigned idx_p			= 6;
const unsigned idx_p_t			= 7;
const unsigned idx_p_v			= 8;
const unsigned idx_kernel_type	= 9;
const unsigned idx_nearest_k	= 10;


const double   infinity         = 1e300;

const double   cyclen           = 1;

void
mexFunction(int nout, mxArray *out[],
            int nin, const mxArray *in[])
{
	/* ------------------------------------------------------------------
	**                                                Check the arguments
	** --------------------------------------------------------------- */

	if (nin == 0) {
		mexPrintf ( "h = SPMatch4( Samples,   SamLocations, SamScales, ...\n"
		            "              Exemplars, ExeLocations, ExeScales, ...\n" 
					"              p, p_t, p_v, kernel_type )\n"
					"   kernel_type: 0 - Cauchy (default), 1 - Gaussian\n" );
		return;
	}

	if (nin <8 ) {
		mexErrMsgTxt("At least 8 input arguments.");
    } else if (nin > 12 ) {
		mexErrMsgTxt("At most 11 input arguments.");
	} else if (nout != 1) {
		mexErrMsgTxt("Exact 1 output argument.");
	}

    for (int i=0; i<nin; i++) {
    	if(!mxIsDouble(in[i])) 
            mexErrMsgTxt("Every argument should be double.");
    }
    
    double p, p_t, p_v;
    unsigned nearest_k;
	int kernel_type;
    
	if (nin<=idx_kernel_type) {
		kernel_type = 0;
	} else {
		kernel_type = (int)mxGetScalar( in[idx_kernel_type] );
	}
    
    if (nin<=idx_p_v) {
        p_v =  1;
    } else {
        p_v = *mxGetPr(in[idx_p_v]);
    }
    if (nin<=idx_p_t) {
        p_t = 8;
    } else {
        p_t = *mxGetPr(in[idx_p_t]);
    }
    if (nin<=idx_p) {
        p = 0.08;
    } else {
        p = *mxGetPr(in[idx_p]);
    }
    if (nin<=idx_nearest_k) {
        nearest_k = UINT_MAX;
    } else {
        double nk = mxGetScalar(in[idx_nearest_k]);
        if ( nk > double(UINT_MAX) || nk<1 ) {
            nearest_k = UINT_MAX;
        } else {
            nearest_k = int( nk );
        }
    }

    
//     printf("%lf,%lf,%lf", p, p_t, p_v);
    
    const mxArray * Samples,   * SamLocations, * SamScales, 
                  * Exemplars, * ExeLocations, * ExeScales;
    Samples         = in[idxSamples];
	SamLocations    = in[idxSamLocations];
	SamScales       = in[idxSamScales];
	Exemplars       = in[idxExemplars];
	ExeLocations    = in[idxExeLocations];
	ExeScales       = in[idxExeScales];
    
    const double  * dSamples,   * dSamLocations, * dSamScales, 
                  * dExemplars, * dExeLocations, * dExeScales;
    dSamples        = mxGetPr(Samples);
	dSamLocations   = mxGetPr(SamLocations);
	dSamScales      = mxGetPr(SamScales);
	dExemplars      = mxGetPr(Exemplars);
	dExeLocations   = mxGetPr(ExeLocations);
	dExeScales      = mxGetPr(ExeScales);
	/* ------------------------------------------------------------------
	**                                                         Likelihood
	** --------------------------------------------------------------- */
	
	mxArray* _h = mxCreateNumericMatrix( 1, 1, mxDOUBLE_CLASS, mxREAL );
	out[0] = _h;
    
    double&  h = *(mxGetPr(_h));
    h = 0;

    const double pi = 3.1415926535897931159979634685441851615906;
    const double   lowest          = -infinity;
        
    double detaNomi = 4*p/pi;
//     double detaNomi = 1/(sqrt(2*pi)*p);
    double logdetaNomi = log(detaNomi);
    double p42 = 4*p*p;
    const double log2 = log(2.0);
    double gauLOGnormalizer = log(sqrt(2*pi)*p);
    
    unsigned exeNum, samNum, FeatureLen;
    exeNum = (unsigned)mxGetM(Exemplars);
    samNum = (unsigned)mxGetM(Samples);
    FeatureLen = (unsigned)mxGetN(Samples);

   
    for (unsigned k=0; k<samNum; k++) {
        const double* theSample = dSamples + k;
        double  theSamLoc = dSamLocations[k], 
                theSamSca = dSamScales[k];
        unsigned nonzeroNum = 0;

        double  mh = 0;
      
        vector<double> all_like(exeNum);
        for (unsigned i=0; i<exeNum; i++) {
            
			if ( p_t>0 ) {
				double deltaLoc;
                deltaLoc = abs(dExeLocations[i] - theSamLoc);
				if ( deltaLoc > cyclen/2 )
					deltaLoc = cyclen - deltaLoc;
				if ( deltaLoc >= p_t )
                    continue;
			}
			if ( p_v>0 ) {
				double deltaSca;
				deltaSca  = abs(log(dExeScales[i] / theSamSca)/log2);
				if ( deltaSca >= p_v )
					continue;
			}
            
            nonzeroNum++;
            
            double like = 0;
            const double* theExemplar = dExemplars + i;
			double factor_const = 0;
			switch ( kernel_type ) {
				case 0:		//Cauchy
				{
					double t = 1.0;
					const unsigned t_total = 16;
					factor_const = logdetaNomi;
					for (unsigned j=0; j<FeatureLen; j++) {
						double fd;
						fd = theSample[j*samNum] - theExemplar[j*exeNum];
						t*=fd*fd+p42;
						if ( !(j%t_total) ) {
							like += -log(t);
							t = 1.0;
						}
						//like += - log(fd*fd+p42);
					}
					like += -log(t);
					break;
				}
				case 1:		//Gaussian
				{
					factor_const = -gauLOGnormalizer;
					for (unsigned j=0; j<FeatureLen; j++) {
						double fd;
						fd = theSample[j*samNum] - theExemplar[j*exeNum];
						like += -(fd*fd);
					}
					like = like / (p42/2);
					break;
				}
			}
			like += factor_const*FeatureLen;

            all_like[i] = like;
        }
        
        sort( all_like.begin(), all_like.end(), std::greater<double>() );
        
        unsigned usedExeNum = std::min(exeNum, nearest_k);
        for ( unsigned i=0; i<usedExeNum; ++i ) {
            mh += exp(all_like[i]);
        }

        if (mh>0) {
            mh = log(mh/(double)usedExeNum);
        }
        
        if ( mh > lowest ) {
            h += mh;
        }
    }
    
    h = h/(samNum);

}