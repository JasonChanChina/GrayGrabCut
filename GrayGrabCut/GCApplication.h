#include "GrayGrabcut.h"
#include "SLICO.h"
#include <iostream>

using namespace std;


class GCApplication
{
public:
    enum{ STATE_UP = 0, STATE_DOWN = 1};
    static const int radius = 2;
    static const int thickness = -1;

    void reset();
	void clearMarkInImage();
    void setImageAndWinName( const Mat& _image, const string& _winName );

	void show();
    void showImage();
	void showSuperImage();
    void mouseClick( int event, int x, int y, int flags, void* param );

	int nextIter();
	int nextSuperIter();
    int getIterCount() const;

	void updateMaskFromRectsAndPixels();

	void getBinMask(Mat& comMask, Mat& binMask );

	void SPMask2Mask();
	void mask2SPMask();

	int superpixelSegmentation(int _step);

	
	void dyeInvalidRegion();

	

private:
    void setCircleInMask( int flags, Point p, bool isPr );

    string winName;
	string superName;
	bool isFirst;

	uchar leftState, rightState;

	vector<Point> fgdPxls, bgdPxls, prFgdPxls, prBgdPxls;
	vector<Rect> rects;
	int iterCount;

	GrayGrabCut grab;
	SLICO slico;

	///////////////////////
	Mat bgdModel, fgdModel;

	// pixel
	Mat originImage;
    Mat image;
    Mat mask;
    
	// connection
	vector<int> pixelLabels;
	int kindOfLabels; 
	vector<vector<int> > contains;

	// superpixel
    vector<Vec3d> centers;
	vector<vector<int> > arcs;
	vector<int> SPMask;

	Mat originSuperImage;
	Mat superImage;
	double beta;

	Point beginPoint, middlePoint, endPoint;

};
