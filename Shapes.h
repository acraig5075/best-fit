#pragma once

#include "BestFit.h"

/***********************************************************
************************************************************
**********************  BestFitLine ************************
************************************************************
***********************************************************/

class BestFitLine : public BestFit
{
public:
	BestFitLine(int observations);

private:
	void GenerateProvisionals();
	void FormulateMatrices();
	double SolveAt(double x, double y) const;
	void EvaluateFinalResiduals(int point, double &vxi, double &vyi) const;
};

/***********************************************************
************************************************************
**********************  BestFitCircle **********************
************************************************************
***********************************************************/

class BestFitCircle : public BestFit
{
public:
	BestFitCircle(int observations);

private:
	void GenerateProvisionals();
	void FormulateMatrices();
	double SolveAt(double x, double y) const;
};

/***********************************************************
************************************************************
**********************  BestFitEllipse *********************
************************************************************
***********************************************************/

class BestFitEllipse : public BestFit
{
public:
	BestFitEllipse(int observations);

private:
	void GenerateProvisionals();
	void FormulateMatrices();
	double SolveAt(double x, double y) const;
};

/***********************************************************
************************************************************
**********************  BestFitFactory *********************
************************************************************
***********************************************************/

class BestFitFactory
{
public:
	BestFitFactory();
	static BestFit *Create(int type, int observations);
};
