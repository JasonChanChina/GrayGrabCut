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

void Histogram::createSuperHist(vector<Vec3d>& centers, vector<vector<int> >& contains, vector<int>& SPMask)
{

	int histSize = 256;

	hist.create(1, histSize, CV_32FC1); 
	int count = centers.size();
	int maskCount = 0;

	vector<int> histCount(histSize, 0);

	for(int i = 0; i < count; i++)
	{
		if(SPMask[i] > 0)
		{
			int gray = (int)(centers[i][0]+0.5);
			if(gray > 255) gray = 255;
			if(gray < 0) gray = 0;
			int singleGrayCount = contains[i].size();
			histCount[gray] += singleGrayCount;
			maskCount += singleGrayCount;
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

	////直方图概率改进
	//// 假设：
	////领域： -2   -1  0   +1  +2
	////概率： 0.1 0.2 0.4 0.2 0.1

	//vector<double> histPro(5, 0);
	//histPro[0] = 0.1;
	//histPro[1] = 0.2;
	//histPro[2] = 0.4;
	//histPro[3] = 0.2;
	//histPro[4] = 0.1;

	//int start = gray - 2;
	//int end = gray + 3;
	//double pro = 0;
	//int j = 0;
	//for(int i = start; i < end; i++)
	//{
	//	if(i >= 0 && i <= 255)
	//	{
	//		pro += histPro[j]*hist.at<float>(i);
	//	}
	//	j++;
	//}

	//return pro;


}