#ifndef _splinesDUNE_h_
#define _splinesDUNE_h_

#include "splines/splineFDBase.h"

class splinesDUNE : virtual public splineFDBase
{
  public:
	splinesDUNE(const char *spline, int nutype, int nevents, int DetID, covarianceXsec* xsec_cov = NULL); // constructor for etrue-var1 splines
	splinesDUNE(const char *spline, int nutype, int nevents, double BinningOpt, int DetID, covarianceXsec* xsec_cov = NULL); // constructor for etrue-var1-var2 splines
	virtual ~splinesDUNE(){};
	void SetupSplines();
	void SetupSplines(int BinningOpt);

	void FindUniqueModes();
};

#endif
