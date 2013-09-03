#include <complex>
#include "BestFit.h"

const double Double::Accuracy = 0.000000001;
	
bool Double::IsZero(double a)
{
	return fabs(a) <= Accuracy;
}
bool Double::IsNotZero(double a)
{
	return fabs(a) > Accuracy;
}
bool Double::AreEqual(double a, double b)
{
	return fabs(a - b) <= (Accuracy + Accuracy * fabs(b));
}
bool Double::AreNotEqual(double a, double b)
{
	return fabs(a - b) > (Accuracy + Accuracy * fabs(b));
}
bool Double::IsLessThan(double a, double b)
{
	return a < (b - (Accuracy + Accuracy * fabs(b)));
}
bool Double::IsGreaterThan(double a, double b)
{
	return a > (b + (Accuracy + Accuracy * fabs(b)));
}
bool Double::IsLessThanOrEquals(double a, double b)
{
	return a <= (b + (Accuracy + Accuracy * fabs(b)));
}
bool Double::IsGreaterThanOrEquals(double a, double b)
{
	return a >= (b - (Accuracy + Accuracy * fabs(b)));
}
bool Double::IsLessThanZero(double a)
{
	return a < -Accuracy;
}
bool Double::IsGreaterThanZero(double a)
{
	return a > Accuracy;
}
bool Double::IsLessThanOrEqualsZero(double a)
{
	return a <= Accuracy;
}
bool Double::IsGreaterThanOrEqualsZero(double a)
{
	return a >= -Accuracy;
}
