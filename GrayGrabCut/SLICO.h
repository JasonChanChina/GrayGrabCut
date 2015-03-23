#pragma once

#include "stdafx.h"


class SLICO
{
public:
	SLICO();
	virtual ~SLICO();


	void DoSuperpixelSegmentation_ForGivenMat( const Mat& img, vector<int>& pixelLabels, int& kindOfLabels, int _step);

	void PerformSLICO_ForGivenStepSize(
		const unsigned int*			ubuff,
		const int					width,
		const int					height,
		int*						pixelLabels,
		int&						kindOfLabels,
		const int&					STEP,
		const double&				m);

	void PerformSLICO_ForGivenK(
		const unsigned int*			ubuff,
		const int					width,
		const int					height,
		int*						pixelLabels,
		int&						kindOfLabels,
		const int&					K,
		const double&				m);

	void GetArcAndCenterOfSuperpixels( const Mat& img, vector<int>& pixelLabels, int& kindOfLabels, vector<vector<int> >& arcs, vector<Vec3d>& centers, vector<vector<int> >& spLabels);

	void DrawContoursAroundSegments(Mat& mat, vector<int>& pixelLabels, Scalar color);
	//void DrawAverageColor(Mat &mat, vector<int>& pixelLabels, vector<Vec6d>& centers);

private:

	void PerformSuperpixelSegmentation_VariableSandM(
		vector<double>&				kseedsg,
		vector<double>&				kseedsx,
		vector<double>&				kseedsy,
		int*						pixelLabels,
		const int&					STEP,
		const int&					NUMITR);

	void GetLABXYSeeds_ForGivenStepSize(
		vector<double>&				kseedsg,
		vector<double>&				kseedsx,
		vector<double>&				kseedsy,
		const int&					STEP,
		const bool&					perturbseeds,
		const vector<double>&		edgemag);

	void GetLABXYSeeds_ForGivenK(
		vector<double>&				kseedsg,
		vector<double>&				kseedsx,
		vector<double>&				kseedsy,
		const int&					STEP,
		const bool&					perturbseeds,
		const vector<double>&		edges);


	void PerturbSeeds(
		vector<double>&				kseedsg,
		vector<double>&				kseedsx,
		vector<double>&				kseedsy,
		const vector<double>&		edges);

	void DetectLabEdges(
		const double*				gvec,
		const int&					width,
		const int&					height,
		vector<double>&				edges);

	void DoRGBtoLABConversion(
		const unsigned int*&		ubuff,
		double*&					gvec);



	void EnforceLabelConnectivity(
		const int*					labels,
		const int&					width,
		const int&					height,
		int*						nlabels,//input labels that need to be corrected to remove stray labels
		int&						kindOfLabels,//the number of labels changes in the end if segments are removed
		const int&					K); //the number of superpixels desired by the user


private:
	int										m_width;
	int										m_height;
	int										m_depth;

	int										m_size;

	double*									m_gvec;
	//double*									m_lvec;
	//double*									m_avec;
	//double*									m_bvec;

	double**								m_gvecvec;
	//double**								m_lvecvec;
	//double**								m_avecvec;
	//double**								m_bvecvec;

};

