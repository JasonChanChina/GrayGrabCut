
#pragma once

#include "stdafx.h"

#include "Histogram.h"
#include "gcgraph.hpp"


 

class GrayGrabCut
{
public:
	double calcBeta(const Mat& img);
	void calcNWeights( const Mat& img, Mat& leftW, Mat& upleftW, Mat& upW, Mat& uprightW, double beta, double gamma );
	void checkMask( const Mat& img, const Mat& mask );
	void initMaskWithRect( Mat& mask, Size imgSize, Rect rect );
	void initHists( const Mat& img, const Mat& mask, Histogram& bgdHist, Histogram& fgdHist );
	void constructGCGraph( const Mat& img, const Mat& mask, const Histogram& bgdHist, const Histogram& fgdHist, double lambda,
                       const Mat& leftW, const Mat& upleftW, const Mat& upW, const Mat& uprightW,
                       GCGraph<double>& graph );
	void estimateSegmentation( GCGraph<double>& graph, Mat& mask );
	void graygrabCut( InputArray _img, InputOutputArray _mask, Rect rect,
                  InputOutputArray _bgdModel, InputOutputArray _fgdModel,
                  int iterCount, int mode=GC_EVAL );

private:
	void editMask(Mat& mask, bool isFgd);

};
