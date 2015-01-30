
#include "stdafx.h"

class Histogram{

public:
	Histogram(Mat _hist);
	void createHist(const Mat& img, Mat& mask);
	double probability(const uchar gray) const;

private:
	Mat hist;
};