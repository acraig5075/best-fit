#include "BestFit.h"
#include <limits>
#include <cassert>
#include <stdexcept>
#include <iomanip>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#define CONVERGENCE_CRITERIA 0.000000001
#define MAX_ITERATIONS 50

BestFit::BestFit(int observations, int unknowns)
	: m_solution(unknowns,1),
		m_verbosity(1),
		m_count(0),
		m_residuals(observations, 1),
		m_design(observations, unknowns),
		m_l(observations, 1),
		m_qweight(observations, observations, 0, 0),
		m_observations(observations, 2),
		m_provisionals(unknowns,1),
		m_b(observations, observations * 2),
		m_numObs(observations),
		m_numUnknowns(unknowns),
		m_minx(std::numeric_limits<double>::max()),
		m_maxx(std::numeric_limits<double>::min()),
		m_miny(std::numeric_limits<double>::max()),
		m_maxy(std::numeric_limits<double>::min())
{
}

BestFit::~BestFit()
{
}

// Set the level of detail required for output to stdout
void BestFit::SetVerbosity(int verbosity)
{
	m_verbosity = verbosity;
}

// Another another observable (coordinate) to the computation object
void BestFit::AddObservation(double x, double y)
{
	if (m_count >= m_numObs)
		throw std::out_of_range("Adding more than expected");

	m_observations(m_count, 0) = x;
	m_observations(m_count, 1) = y;

	m_count++;
	m_minx = std::min<double>(m_minx, x);
	m_maxx = std::max<double>(m_maxx, x);
	m_miny = std::min<double>(m_miny, y);
	m_maxy = std::max<double>(m_maxy, y);
}

// Do the least-squares adjustment
void BestFit::Compute()
{
	assert(m_count == m_numObs);
	if (m_count != m_numObs)
		return;

	if (m_verbosity > 1)
		std::cout << "Observations:    " << m_observations << std::endl;

	GenerateProvisionals();

	int it = 0;
	while (true)
		{
		FormulateMatrices();

		if (m_verbosity > 1)
			{
			std::cout << "Provisionals:    " << m_provisionals << std::endl;
			std::cout << "l-matrix:        " << m_l << std::endl;
			}

		// evaluate the unknowns - small corrections to be applied to the provisional unknowns
		if (EvaluateUnknowns())
			{
			it++;

			if (HasConverged())
				break;
			if (IsDegenerate(it))
				break;

			// add the unknowns to the provisional values
			EvaluateAdjustedUnknowns();
			}
		else
			break;

		}

	bool successful = (it > 0 && it < MAX_ITERATIONS);
	if (successful)
		{
		if (0 == m_verbosity)
			{
			OutputSimpleSolution();
			}
		else
			{
			// evaluate the residuals
			EvaluateResiduals();

			// add the residuals to the provisional observations
			EvaluateAdjustedObservations();

			// Check that the adjusted unknowns and adjusted observations satisfy
			// the original line/circle/ellipse equation.
			GlobalCheck();

			// Subsequent error analysis and statistical output
			ErrorAnalysis(it);
			}
		}
}

// Matrix inversion routine using LU decomposition
bool BestFit::InvertMatrix(const ublas::matrix<double> &input, ublas::matrix<double> &inverse)
{
	// create a working copy of the input
	ublas::matrix<double> A(input);

	// create a permutation matrix for the LU-factorization
	ublas::permutation_matrix<double> pm(A.size1());

	// perform LU-factorization
	int res = ublas::lu_factorize(A, pm);
	if (res != 0)
		return false;

	// create identity matrix of "inverse"
	inverse.assign(ublas::identity_matrix<double>(A.size1()));

	// back-substitute to get the inverse
	ublas::lu_substitute(A, pm, inverse);

	return true;
}

// Determine solution vector
bool BestFit::EvaluateUnknowns()
{
	ublas::matrix<double> pa = ublas::prod(m_qweight, m_design);
	ublas::matrix<double> atpa = ublas::prod(ublas::trans(m_design), pa);

	ublas::matrix<double> inverse(atpa.size1(), atpa.size2());
	if (InvertMatrix(atpa, inverse))
		{
		ublas::matrix<double> pl = ublas::prod(m_qweight, m_l);
		ublas::matrix<double> atpl = ublas::prod(ublas::trans(m_design), pl);

		m_solution = ublas::prod(inverse, atpl);

		if (m_verbosity > 2)
			std::cout << "Iter. solution   " << m_solution << std::endl;
		return true;
		}
	else
		std::cout << "No solution. Cannot invert matrix." << std::endl;
	return false;
}


// Add solution to the provisionals
void BestFit::EvaluateAdjustedUnknowns()
{
	m_provisionals += m_solution;

	if (m_verbosity > 2)
		std::cout << "Iter adj unknowns" << m_provisionals << std::endl;
}

// Residuals to the initial observations
void BestFit::EvaluateResiduals()
{
	// v = Ax-l, quasi-residuals in the case of circle and ellipse
	m_residuals = ublas::prod(m_design, m_solution) - m_l;

	if (m_verbosity > 1)
		std::cout << "Residuals:       " << m_residuals << std::endl;
}

void BestFit::EvaluateFinalResiduals(int point, double &vxi, double &vyi) const
{
	// In the case of the quasi-parametric case (circle and ellipse) the actual
	// residuals for both the x and y coordinate need to be extracted.
	vxi = m_design(point, 0) * m_qweight(point, point) * m_residuals(point, 0);
	vyi = m_design(point, 1) * m_qweight(point, point) * m_residuals(point, 0);

	//	ublas::matrix<double> bt = ublas::trans(m_b);
	//	ublas::matrix<double> pv = ublas::prod(m_qweight, m_residuals);
	//	ublas::matrix<double> btpv = ublas::prod(bt, pv);
	//	vxi = -btpv(point * 2 + 0,0);
	//	vyi = -btpv(point * 2 + 1,0);
	//	std::cout << avxi - vxi << "," << avyi - vyi << std::endl;
}

// Add residuals to initial observations
void BestFit::EvaluateAdjustedObservations()
{
	for (int i = 0; i < m_numObs; i++)
		{
		double vxi = 0.0;
		double vyi = 0.0;
		EvaluateFinalResiduals(i, vxi, vyi); // overridden for normal (non-quasi-parametric) BestFitLine

		m_observations(i, 0) += vxi;
		m_observations(i, 1) += vyi;
		}

	if (m_verbosity > 1)
		std::cout << "Adj. observations" << m_observations << std::endl;
}

// Global check on the quality of the adjustment
void BestFit::GlobalCheck()
{
	ublas::vector<double> global(m_numObs);
	bool pass = (m_numObs > 0);

	for (int i = 0; i < m_numObs; i++)
		{
		double x = m_observations(i, 0);
		double y = m_observations(i, 1);

		global (i) = SolveAt(x, y); // should be zero
		if (pass)
			//pass = PARMZERO(global(i));
			pass = fabs(global(i)) < 0.01; // TODO: Is this too lax?
		}

	if (m_verbosity > 1)
		std::cout << "Global check:    " << global << std::endl;
	if (m_verbosity > 0)
		std::cout << "Global check of the adjustment     ***" << (pass ? "PASSES***" : "FAILS***") << std::endl;

	// Just for fun let's check that aTPv is zero too, really only neccessary
	// if the above global check fails.
	ublas::matrix<double> atp = ublas::prod(ublas::trans(m_design), m_qweight);
	ublas::matrix<double> atpv = ublas::prod(atp, m_residuals);

	pass = PARMZERO(atpv(0,0));
	if (!pass && m_verbosity > 1)
		std::cout << "aTPv             " << atpv << std::endl;
	if (m_verbosity > 0)
		std::cout << "Check of the evaluated unknowns    ***" << (pass ? "PASSES***" : "FAILS***") << std::endl;
}

// Predicate for seeing whether the solution has become sufficiently small
bool BestFit::HasConverged() const
{
	for (int i = 0; i < m_numUnknowns; i++)
		{
		if (fabs(m_solution(i, 0)) > CONVERGENCE_CRITERIA)
			return false;
		}
	return true;
}

// Predicate for seeing whether the solution is diverging?
bool BestFit::IsDegenerate(int iteration) const
{
	bool degenerate = (iteration >= MAX_ITERATIONS);
	if (degenerate)
		std::cout << "No solution. Does not converge." << std::endl;

	return degenerate;
}

// Output stats
// TODO: Variance-covariance matrices
void BestFit::ErrorAnalysis(int iterations)
{
	ublas::matrix<double> pv = ublas::prod(m_qweight, m_residuals);
	ublas::matrix<double> vtpv = ublas::prod(ublas::trans(m_residuals), pv);
	int dfreedom = m_numObs - m_numUnknowns;
	double variance = vtpv(0, 0) / dfreedom;
	double stddev = sqrt(variance);

	if (m_verbosity > 0)
		{
		std::cout << std::setprecision(6)

			<< std::setiosflags(std::ios::fixed)
			<< "Number of observations             " << m_numObs << std::endl
			<< "Number of unknowns                 " << m_numUnknowns << std::endl
			<< "Degrees of freedom                 " << dfreedom << std::endl
			<< "Iterations until convergence       " << iterations << std::endl
			<< "Variance                           " << variance << std::endl
			<< "Std. dev. observation unit weight  " << stddev << std::endl;

		std::cout << "***********************************" << std::endl;
		for (int i = 0; i < m_numUnknowns; i++)
			{
			std::cout << "Adj. unknown " << i + 1
				<< "                     "
				<< m_provisionals(i, 0) << std::endl;
			}
		}
}

// Output at it's very simplest
void BestFit::OutputSimpleSolution() const
{
	for (int i = 0; i < m_numUnknowns; i++)
		std::cout << m_provisionals(i, 0) << std::endl;
}
