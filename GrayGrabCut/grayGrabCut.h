
#pragma once

#include "stdafx.h"

#include "Histogram.h"
#include "Superpixel.h"
#include "gcgraph.hpp"




class GrayGrabCut
{
public:
	double calcBeta(const Mat& img);
	void checkMask( const Mat& img, const Mat& mask );
	void initMaskWithRect( Mat& mask, Size imgSize, Rect rect );

	//pixel cut
	void calcNWeights( const Mat& img, Mat& leftW, Mat& upleftW, Mat& upW, Mat& uprightW, double beta, double gamma );
	void initHists( const Mat& img, const Mat& mask, Histogram& bgdHist, Histogram& fgdHist );
	void constructGCGraph( const Mat& img, const Mat& mask, const Histogram& bgdHist, const Histogram& fgdHist, double lambda,
		const Mat& leftW, const Mat& upleftW, const Mat& upW, const Mat& uprightW, GCGraph<double>& graph );
	void estimateSegmentation( GCGraph<double>& graph, Mat& mask );
	void graygrabCut( InputArray _img, InputOutputArray _mask, Rect rect,
		InputOutputArray _bgdModel, InputOutputArray _fgdModel, int iterCount, int mode=GC_EVAL );


	//superpixel cut
	void calcSuperNWeights(vector<Vec3d>& centers, vector<vector<int> >& arcs, vector<vector<int> >& contains,
		vector<vector<double> >& nweights, int& sideCount, double beta, double gamma );
	void initSuperHists(vector<Vec3d>& centers, vector<int>& SPMask, Histogram& bgdHist, Histogram& fgdHist );
	void constructSuperGCGraph(vector<Vec3d> centers, vector<vector<int> > arcs, vector<int>& SPMask,vector<vector<double> >& nweights, 
		Histogram& fgdHist ,Histogram& bgdHist , double lambda, int sideCount, vector<int>& vtxs,  GCGraph<double>& graph);
	void estimateSuperSegmentation( GCGraph<double>& graph, vector<int>& SPMask , vector<int> vtxs);
	void graySupergrabCut(vector<int> pixelLabels, int kindOfLabels, vector<vector<int> > arcs, vector<Vec3d> centers, 
		vector<vector<int> >& contains, vector<int>& SPMask, Mat& bgdModel, Mat& fgdModel, int iterCount, double beta, bool firstTime);


	void SPMask2Mask(vector<int>& SPMask, Mat& mask, vector<vector<int> >& contains);
	void Mask2SPMask( Mat& mask, vector<int>& SPMask, vector<vector<int> >& contains, vector<int> pixelLabels );
private:
	void editMask(Mat& mask, bool isFgd);
	void editSPMask(vector<int>& SPMask, bool isFgd);

};
