#include "Histogram.h"

Histogram::Histogram(Mat _hist)
{
	hist = _hist;
}

void Histogram::createHist(const Mat& img, Mat& mask)
{
	//int histSize = 256;
	//float range[] = {0, 255};
	//const float* histRange = {range};
	//bool uniform = false;	//是否归一化/均衡化
	//bool accumulate = false;	//是否累计之前记录
	//calcHist( &img, 1, 0, mask, hist, 1, &histSize, &histRange, uniform, accumulate );

	int histSize = 256;

	hist.create(1, histSize, CV_32FC1); 
	int count = img.cols * img.rows;

	vector<int> histTotal(histSize, 0);
	for( int y = 0; y < img.rows; y++ )
    {
        for( int x = 0; x < img.cols; x++ )
        {
			if(mask.at<uchar>(y,x) > 0)
			{
				uchar value = img.at<uchar>(y,x);
				histTotal[value]++;
			}
		}
	}

	for(int i = 0; i < histTotal.size(); i++)
	{
		hist.at<float>(0,i) = (float)histTotal[i] / count;
	}

}

double Histogram::probability(const uchar gray) const
{
	return hist.at<float>(gray);
}