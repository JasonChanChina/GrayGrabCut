//#include "opencv2/highgui/highgui.hpp"
//#include "opencv2/imgproc/imgproc.hpp"
#include "GrayGrabcut.h"
#include <iostream>

using namespace std;
//using namespace cv;

static void help()
{
    cout << "\nThis program demonstrates GrabCut segmentation -- select an object in a region\n"
            "and then grabcut will attempt to segment it out.\n"
            "Call:\n"
            "./grabcut <image_name>\n"
        "\nSelect a rectangular area around the object you want to segment\n" <<
        "\nHot keys: \n"
        "\tESC - quit the program\n"
        "\tr - restore the original image\n"
		"\tp - clear invalid region\n"
        "\tn - next iteration\n"
        "\n"
        "\tleft mouse button - set rectangle\n"
        "\n"
        "\tCTRL+left mouse button - set GC_BGD pixels\n"
        "\tSHIFT+left mouse button - set CG_FGD pixels\n"
        "\n"
        "\tCTRL+right mouse button - set GC_PR_BGD pixels\n"
        "\tSHIFT+right mouse button - set CG_PR_FGD pixels\n" << endl;
}

const Scalar RED = Scalar(0,0,255);
const Scalar PINK = Scalar(230,130,255);
const Scalar BLUE = Scalar(255,0,0);
const Scalar LIGHTBLUE = Scalar(255,255,160);
const Scalar GREEN = Scalar(0,255,0);

const int BGD_KEY = CV_EVENT_FLAG_CTRLKEY;
const int FGD_KEY = CV_EVENT_FLAG_SHIFTKEY;

const Scalar COLOR_RECT = Scalar(255,255,255);
const Scalar COLOR_FGD = Scalar(128,128,128);
const Scalar COLOR_BGD = Scalar(0,0,0);
const Scalar COLOR_PR_FGD = Scalar(192,192,192);
const Scalar COLOR_PR_BGD = Scalar(64,64,64);


static void getBinMask( const Mat& comMask, Mat& binMask )
{
    if( comMask.empty() || comMask.type()!=CV_8UC1 )
        CV_Error( CV_StsBadArg, "comMask is empty or has incorrect type (not CV_8UC1)" );
    if( binMask.empty() || binMask.rows!=comMask.rows || binMask.cols!=comMask.cols )
        binMask.create( comMask.size(), CV_8UC1 );
    binMask = comMask & 1;
}

class GCApplication
{
public:
    enum{ NOT_SET = 0, IN_PROCESS = 1, SET = 2 };
    static const int radius = 2;
    static const int thickness = -1;

    void reset();
    void setImageAndWinName(Mat& _image, const string& _winName );
    void showImage() const;
    void mouseClick( int event, int x, int y, int flags, void* param );
    int nextIter();
    int getIterCount() const { return iterCount; }

	void dyeInvalidRegion();

private:
    void setRectInMask();
    void setLblsInMask( int flags, Point p, bool isPr );

    const string* winName;
    Mat* image;
    Mat mask;
    Mat bgdModel, fgdModel;

    uchar rectState, lblsState, prLblsState;
    bool isInitialized;

    Rect rect;
    vector<Point> fgdPxls, bgdPxls, prFgdPxls, prBgdPxls;
    int iterCount;
	GrayGrabCut grab;
};

void GCApplication::reset()
{
    if( !mask.empty() )
        mask.setTo(Scalar::all(GC_BGD));
    bgdPxls.clear(); fgdPxls.clear();
    prBgdPxls.clear();  prFgdPxls.clear();

    isInitialized = false;
    rectState = NOT_SET;
    lblsState = NOT_SET;
    prLblsState = NOT_SET;
    iterCount = 0;
}

void GCApplication::setImageAndWinName( Mat& _image, const string& _winName  )
{
    if( _image.empty() || _winName.empty() )
        return;
    image = &_image;
    winName = &_winName;
    mask.create( image->size(), CV_8UC1);
    reset();
}

void GCApplication::showImage() const
{
    if( image->empty() || winName->empty() )
        return;

    Mat res;
	res.create(image->size(), image->type());
	res.setTo(Vec3b(255,255,255));				//���ÿհ�������ɫ


    Mat binMask;
    if( !isInitialized )
        image->copyTo( res );
    else
    {
        getBinMask( mask, binMask );
        image->copyTo( res, binMask );
    }

    vector<Point>::const_iterator it;
    for( it = bgdPxls.begin(); it != bgdPxls.end(); ++it )
        circle( res, *it, radius, COLOR_BGD, thickness );
    for( it = fgdPxls.begin(); it != fgdPxls.end(); ++it )
        circle( res, *it, radius, COLOR_FGD, thickness );
    for( it = prBgdPxls.begin(); it != prBgdPxls.end(); ++it )
        circle( res, *it, radius, COLOR_PR_BGD, thickness );
    for( it = prFgdPxls.begin(); it != prFgdPxls.end(); ++it )
        circle( res, *it, radius, COLOR_PR_FGD, thickness );

    if( rectState == IN_PROCESS || rectState == SET )
        rectangle( res, Point( rect.x, rect.y ), Point(rect.x + rect.width, rect.y + rect.height ), COLOR_RECT, 2);

    imshow( *winName, res );
	string path = "D:\\aa_";
	path += "_a.png";
	imwrite(path, res);


}

void GCApplication::setRectInMask()
{
    assert( !mask.empty() );
    mask.setTo( GC_BGD );
    rect.x = max(0, rect.x);
    rect.y = max(0, rect.y);
    rect.width = min(rect.width, image->cols-rect.x);
    rect.height = min(rect.height, image->rows-rect.y);
    (mask(rect)).setTo( Scalar(GC_PR_FGD) );
}

void GCApplication::setLblsInMask( int flags, Point p, bool isPr )
{
    vector<Point> *bpxls, *fpxls;
    uchar bvalue, fvalue;
    if( !isPr )
    {
        bpxls = &bgdPxls;
        fpxls = &fgdPxls;
        bvalue = GC_BGD;
        fvalue = GC_FGD;
    }
    else
    {
        bpxls = &prBgdPxls;
        fpxls = &prFgdPxls;
        bvalue = GC_PR_BGD;
        fvalue = GC_PR_FGD;
    }
    if( flags & BGD_KEY )
    {
        bpxls->push_back(p);
        circle( mask, p, radius, bvalue, thickness );
    }
    if( flags & FGD_KEY )
    {
        fpxls->push_back(p);
        circle( mask, p, radius, fvalue, thickness );
    }
}

void GCApplication::mouseClick( int event, int x, int y, int flags, void* )
{
    // TODO add bad args check
    switch( event )
    {
    case CV_EVENT_LBUTTONDOWN: // set rect or GC_BGD(GC_FGD) labels
        {
            bool isb = (flags & BGD_KEY) != 0,
                 isf = (flags & FGD_KEY) != 0;
            if( rectState == NOT_SET && !isb && !isf )
            {
                rectState = IN_PROCESS;
                rect = Rect( x, y, 1, 1 );
            }
            if ( (isb || isf) && rectState == SET )
                lblsState = IN_PROCESS;
        }
        break;
    case CV_EVENT_RBUTTONDOWN: // set GC_PR_BGD(GC_PR_FGD) labels
        {
            bool isb = (flags & BGD_KEY) != 0,
                 isf = (flags & FGD_KEY) != 0;
            if ( (isb || isf) && rectState == SET )
                prLblsState = IN_PROCESS;
        }
        break;
    case CV_EVENT_LBUTTONUP:
        if( rectState == IN_PROCESS )
        {
            rect = Rect( Point(rect.x, rect.y), Point(x,y) );

			/////////////
			//rect = Rect(29,62,262,417);


			cout<<rect.x<<","<<rect.y<<","<<rect.width<<","<<rect.height<<endl;

            rectState = SET;
            setRectInMask();
            assert( bgdPxls.empty() && fgdPxls.empty() && prBgdPxls.empty() && prFgdPxls.empty() );
            showImage();
        }
        if( lblsState == IN_PROCESS )
        {
            setLblsInMask(flags, Point(x,y), false);
            lblsState = SET;
            showImage();
        }
        break;
    case CV_EVENT_RBUTTONUP:
        if( prLblsState == IN_PROCESS )
        {
            setLblsInMask(flags, Point(x,y), true);
            prLblsState = SET;
            showImage();
        }
        break;
    case CV_EVENT_MOUSEMOVE:
        if( rectState == IN_PROCESS )
        {
            rect = Rect( Point(rect.x, rect.y), Point(x,y) );
            assert( bgdPxls.empty() && fgdPxls.empty() && prBgdPxls.empty() && prFgdPxls.empty() );
            showImage();
        }
        else if( lblsState == IN_PROCESS )
        {
            setLblsInMask(flags, Point(x,y), false);
            showImage();
        }
        else if( prLblsState == IN_PROCESS )
        {
            setLblsInMask(flags, Point(x,y), true);
            showImage();
        }
        break;
    }
}

int GCApplication::nextIter()
{
	double begin = (double)getTickCount();
    if( isInitialized )
		grab.graygrabCut( *image, mask, rect, bgdModel, fgdModel, 1 );
    else
    {
        if( rectState != SET )
            return iterCount;

        if( lblsState == SET || prLblsState == SET )
            grab.graygrabCut( *image, mask, rect, bgdModel, fgdModel, 1, GC_INIT_WITH_MASK );
        else
            grab.graygrabCut( *image, mask, rect, bgdModel, fgdModel, 1, GC_INIT_WITH_RECT );

        isInitialized = true;
    }

	double end = (double)getTickCount();
	double time = (end-begin)/getTickFrequency();
	cout<<"Time="<<time<<endl;

    iterCount++;

    bgdPxls.clear(); fgdPxls.clear();
    prBgdPxls.clear(); prFgdPxls.clear();


    return iterCount;
}


void GCApplication::dyeInvalidRegion()
{
	int width = image->cols;
	int height = image->rows;
	

	Mat dyeMask;	//binMask
	dyeMask.create( image->size(), CV_8UC1);
	dyeMask.setTo(Scalar::all(0));

	if(width % 2 != 0 && height % 2 != 0)
	{
		Point point((int)(width+1)/2, (int)(height+1)/2);
		int dyeRadius = point.x;
		circle( dyeMask, point, dyeRadius, 1, thickness );
	}else
	{

		Point pointLU( (int)(width-1)/2, (int)(height-1)/2  );		//left + up
		Point pointRU( (int)(width+1)/2, (int)(height-1)/2  );		//right + up
		Point pointLD( (int)(width-1)/2, (int)(height+1)/2  );		//left + down
		Point pointRD( (int)(width+1)/2, (int)(height+1)/2  );		//right + down

		int dyeRadius = pointLU.x;
		circle( dyeMask, pointLU, dyeRadius, 1, thickness );
		circle( dyeMask, pointRU, dyeRadius, 1, thickness );
		circle( dyeMask, pointLD, dyeRadius, 1, thickness );
		circle( dyeMask, pointRD, dyeRadius, 1, thickness );
	}

	//Mat binMask;
	//getBinMask( dyeMask, binMask );
	Mat imageClone = image->clone();
	image->setTo(Vec3b(255,255,255));

	//Mat res;
	//res.create(image->size(), image->type());
	//res.setTo(Vec3b(255,255,255));				//���ÿհ�������ɫ
	
	imageClone.copyTo( *image, dyeMask );		//���Բ������



	//����һȦ��ɫԲ�α߽�
	if(width % 2 != 0 && height % 2 != 0)
	{
		Point point((int)(width+1)/2, (int)(height+1)/2);
		int dyeRadius = point.x;
		circle( *image, point, dyeRadius, Scalar(0,0,0), 1 );
	}else
	{

		Point pointLU( (int)(width-1)/2, (int)(height-1)/2  );		//left + up
		Point pointRU( (int)(width+1)/2, (int)(height-1)/2  );		//right + up
		Point pointLD( (int)(width-1)/2, (int)(height+1)/2  );		//left + down
		Point pointRD( (int)(width+1)/2, (int)(height+1)/2  );		//right + down

		int dyeRadius = pointLU.x;
		circle( *image, pointLU, dyeRadius, Scalar(0,0,0), 1 );
		circle( *image, pointRU, dyeRadius, Scalar(0,0,0), 1 );
		circle( *image, pointLD, dyeRadius, Scalar(0,0,0), 1 );
		circle( *image, pointRD, dyeRadius, Scalar(0,0,0), 1 );
	}




}


GCApplication gcapp;

static void on_mouse( int event, int x, int y, int flags, void* param )
{
    gcapp.mouseClick( event, x, y, flags, param );
}

int main()
{

	string filename;// = "test.jpg";
	cout<<"input your image name : ";
	cin>>filename;


    Mat image = imread( filename, 1);
    if( image.empty() )
    {
        cout << "\n Durn, couldn't read image filename " << filename << endl;
        return 1;
    }

    help();

    const string winName = "image";
    namedWindow( winName, WINDOW_AUTOSIZE );
    setMouseCallback( winName, on_mouse, 0 );

    gcapp.setImageAndWinName( image, winName );
    gcapp.showImage();

    for(;;)
    {
        int c = waitKey(0);
        switch( (char) c )
        {
        case '\x1b':
            cout << "Exiting ..." << endl;
            goto exit_main;
        case 'r':
            cout << endl;
            gcapp.reset();
            gcapp.showImage();
            break;
		case 'p':
			gcapp.dyeInvalidRegion();
			cout<<"clear invalid region"<<endl;
			gcapp.showImage();
			break;
        case 'n':
            int iterCount = gcapp.getIterCount();
            cout << "<" << iterCount << "... ";
            int newIterCount = gcapp.nextIter();
            if( newIterCount > iterCount )
            {
                gcapp.showImage();
                cout << iterCount << ">" << endl;
            }
            else
                cout << "rect must be determined>" << endl;
            break;
        }
    }

exit_main:
    destroyWindow( winName );
    return 0;
}
