
#include "stdafx.h"

class Histogram{

public:
	Histogram(Mat _hist);
	void createHist(const Mat& img, Mat& mask);
	void createSuperHist(vector<Vec3d>& centers,vector<vector<int> >& contains, vector<int>& SPMask);
	double probability(const uchar gray) const;

private:
	Mat hist;
};