#include "Histogram.h"

Histogram::Histogram(Mat _hist)
{
	hist = _hist;
}

void Histogram::createHist(const Mat& img, const Mat& mask)
{
	int histSize = 255;
	float range[] = {0, 255};
	const float* histRange = {range};
	bool uniform = false;	//是否归一化/均衡化
	bool accumulate = false;	//是否累计之前记录

	calcHist( &img, 1, 0, mask, hist, 1, &histSize, &histRange, uniform, accumulate );

}

double Histogram::probability(const uchar gray) const
{
	return hist.at<double>(gray);
}