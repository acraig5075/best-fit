#pragma once

#include <boost/numeric/ublas/banded.hpp>

namespace ublas = boost::numeric::ublas;

#define ROUNDOFF       (0.000000001)
#define PARMGE(a,b)    ((a) >= ((b) - (ROUNDOFF + ROUNDOFF * fabs((b)))))
#define PARMZERO(a)    (fabs((a)) <= ROUNDOFF)
#define PARMLT(a,b)    ((a) < ((b) - (ROUNDOFF + ROUNDOFF * fabs((b)))))
#define PARMLTZERO(a)  ((a) < -ROUNDOFF)
#define PARMEQ(a,b)    (fabs((a) - (b)) <= (ROUNDOFF + ROUNDOFF * fabs((b))))
#define PARMGT(a,b)    ((a) > ((b) + (ROUNDOFF + ROUNDOFF * fabs((b)))))
#define PARMLE(a,b)    ((a) <= ((b) + (ROUNDOFF + ROUNDOFF * fabs((b)))))
#define PARMGTZERO(a)  ((a) > ROUNDOFF)
#define PARMLEZERO(a)  ((a) <= ROUNDOFF)
#define PARMGEZERO(a)  ((a) >= -ROUNDOFF)
#define PARMNOTZERO(a) (fabs((a)) > ROUNDOFF)
#define PARMNE(a,b)    (fabs((a) - (b)) > (ROUNDOFF + ROUNDOFF * fabs((b))))


class BestFit
{
	// Construction
public:
	BestFit(int observations, int unknowns);
	virtual ~BestFit();
private:
protected:

	// Operators
public:
private:
protected:

	// Implementation
public:
	void SetVerbosity(int verbosity);
	void AddObservation(double x, double y);
	void Compute();

private:
	virtual void GenerateProvisionals() = 0;
	virtual void FormulateMatrices() = 0;
	virtual double SolveAt(double x, double y) const = 0;
	virtual void EvaluateFinalResiduals(int point, double &vxi, double &vyi) const;

	bool HasConverged() const;
	bool IsDegenerate(int iteration) const;
	bool InvertMatrix(const ublas::matrix<double> &input, ublas::matrix<double> &inverse);
	bool EvaluateUnknowns();
	void EvaluateResiduals();
	void EvaluateAdjustedUnknowns();
	void EvaluateAdjustedObservations();
	void GlobalCheck();
	void ErrorAnalysis(int iterations);
	void OutputSimpleSolution() const;
protected:

	// Variables
public:
private:
	ublas::matrix<double> m_solution;

	int m_verbosity;
	int m_count;
protected:
	ublas::matrix<double> m_residuals;
	ublas::matrix<double> m_design;
	ublas::matrix<double> m_l;
	ublas::banded_matrix<double> m_qweight;
	ublas::matrix<double> m_observations;
	ublas::matrix<double> m_provisionals;
	ublas::matrix<double> m_b;

	int m_numObs;
	int m_numUnknowns;
	double m_minx;
	double m_maxx;
	double m_miny;
	double m_maxy;
};
