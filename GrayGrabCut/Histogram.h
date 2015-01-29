
#include "stdafx.h"

class Histogram{

public:
	Histogram(Mat _hist);

	void createHist(const Mat& img, const Mat& mask);

	double probability(const uchar gray) const;

private:
	Mat hist;

};