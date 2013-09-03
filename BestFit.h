#pragma once

#include <boost/numeric/ublas/banded.hpp>
#include <iostream>

namespace ublas = boost::numeric::ublas;


/***********************************************************
************************************************************
**********************  BestFitIO **************************
************************************************************
***********************************************************/

struct BestFitIO
{
	int numPoints;
	double *points;
	int verbosity;
	int numOutputFields;
	double outputFields[5];
	bool wantAdjustedObs;
	bool wantResiduals;
	double *residuals;

	enum { LineGradient, LineYIntercept };
	enum { CircleCentreX, CircleCentreY, CircleRadius };
	enum { EllipseCentreX, EllipseCentreY, EllipseMajor, EllipseMinor, EllipseRotation };

	BestFitIO()
	 : numPoints(0)
	 , points(NULL)
	 , verbosity(1)
	 , numOutputFields(0)
	 , wantAdjustedObs(false)
	 , wantResiduals(false)
	 , residuals(NULL)
	 {}
};

/***********************************************************
************************************************************
**********************  BestFit ****************************
************************************************************
***********************************************************/

class BestFit
{
	// Construction
public:
	BestFit(int unknowns, std::ostream &oStream);
	virtual ~BestFit();
private:
protected:

	// Operators
public:
private:
protected:

	// Implementation
public:
	void Compute(BestFitIO &in, BestFitIO &out);
	virtual double SolveAt(double x, double y) const = 0;

private:
	virtual void GenerateProvisionals() = 0;
	virtual void FormulateMatrices() = 0;
	virtual void EvaluateFinalResiduals(int point, double &vxi, double &vyi) const;
	virtual void OutputAdjustedUnknowns(std::ostream &oStream) const = 0;

	void SetVerbosity(int verbosity);
	void Compute();
	void ResizeMatrices();
	void AddObservation(int count, double x, double y);
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
	void FillOutput(BestFitIO &out) const;
protected:

	// Variables
public:
private:
	int m_verbosity;
	std::ostream &m_oStream;
	ublas::matrix<double> m_solution;
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

/***********************************************************
************************************************************
**********************  BestFitLine ************************
************************************************************
***********************************************************/

class BestFitLine : public BestFit
{
public:
	BestFitLine(std::ostream &oStream);

private:
	void GenerateProvisionals();
	void FormulateMatrices();
	double SolveAt(double x, double y) const;
	void EvaluateFinalResiduals(int point, double &vxi, double &vyi) const;
	void OutputAdjustedUnknowns(std::ostream &oStream) const;
};

/***********************************************************
************************************************************
**********************  BestFitCircle **********************
************************************************************
***********************************************************/

class BestFitCircle : public BestFit
{
public:
	BestFitCircle(std::ostream &oStream);

private:
	void GenerateProvisionals();
	void FormulateMatrices();
	double SolveAt(double x, double y) const;
	void OutputAdjustedUnknowns(std::ostream &oStream) const;
};

/***********************************************************
************************************************************
**********************  BestFitEllipse *********************
************************************************************
***********************************************************/

class BestFitEllipse : public BestFit
{
public:
	BestFitEllipse(std::ostream &oStream);

private:
	void GenerateProvisionals();
	void FormulateMatrices();
	double SolveAt(double x, double y) const;
	void OutputAdjustedUnknowns(std::ostream &oStream) const;
};

/***********************************************************
************************************************************
**********************  BestFitFactory *********************
************************************************************
***********************************************************/

struct BestFitFactory
{
	static BestFit *Create(int type, std::ostream &oStream);
};

/***********************************************************
************************************************************
**********************  Double *****************************
************************************************************
***********************************************************/

class Double
{
	static const double Accuracy;

public:
	static bool IsZero(double a);
	static bool IsNotZero(double a);
	static bool AreEqual(double a, double b);
	static bool AreNotEqual(double a, double b);
	static bool IsLessThan(double a, double b);
	static bool IsGreaterThan(double a, double b);
	static bool IsLessThanOrEquals(double a, double b);
	static bool IsGreaterThanOrEquals(double a, double b);
	static bool IsLessThanZero(double a);
	static bool IsGreaterThanZero(double a);
	static bool IsLessThanOrEqualsZero(double a);
	static bool IsGreaterThanOrEqualsZero(double a);
};