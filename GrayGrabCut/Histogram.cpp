#include "Histogram.h"

Histogram::Histogram(Mat _hist)
{
	hist = _hist;
}

void Histogram::createHist(const Mat& img, Mat& mask)
{
	

	int histSize = 256;

	hist.create(1, histSize, CV_32FC1); 
	int count = img.cols * img.rows;
	int maskCount = 0;

	vector<int> histCount(histSize, 0);
	for( int y = 0; y < img.rows; y++ )
    {
        for( int x = 0; x < img.cols; x++ )
        {
			if(mask.at<uchar>(y,x) > 0)
			{
				uchar value = img.at<uchar>(y,x);
				histCount[value]++;
				maskCount++;
			}
		}
	}

	for(int i = 0; i < histCount.size(); i++)
	{
		hist.at<float>(0,i) = (float)histCount[i] / maskCount;
	}


}

double Histogram::probability(const uchar gray) const
{
	return hist.at<float>(gray);
}