//#include "opencv2/highgui/highgui.hpp"
//#include "opencv2/imgproc/imgproc.hpp"

#include "GCApplication.h"

using namespace std;
//using namespace cv;

static void help()
{
    cout << "\nThis program demonstrates GrabCut segmentation -- select an object in a region\n"
            "and then grabcut will attempt to segment it out.\n"
            "Call:\n"
            "./graySuperpixelGrabcut <image_name>\n"
        "\nSelect a rectangular area around the object you want to segment\n" <<
        "\nHot keys: \n"
        "\tESC - quit the program\n"
        "\tr - restore the original image\n"
		"\ts - Superpixel preProcess\n"
		"\tm - next Super iteration\n"
        "\tn - next iteration\n"
        "\n"
        "\tCtrl+left mouse button - set rectangle\n"
        "\n"
		"\tSHIFT+Left mouse button - set GC_PR_FGD pixels\n"
		"\tSHIFT+Right mouse button - set GC_PR_BGD pixels\n"
		"\n"
		"\tALT+Left mouse button - set GC_FGD pixels\n"
		"\tALT+Right mouse button - set GC_BGD pixels\n" << endl;
}





GCApplication gcapp;

static void on_mouse( int event, int x, int y, int flags, void* param )
{
    gcapp.mouseClick( event, x, y, flags, param );
}

int main()
{
	cout<<"\t\tHello, welcome to SLICO_iterGraphCut!"<<endl;
	cout<<"\t(Notice:program and images must be in same place.)"<<endl; 
	cout<<endl;
	cout<<"Please input the gray image name(no space)."<<endl;
	cout<<" -> ";
	string filename;
	cin>>filename;
	cout<<"======================================="<<endl;
	help();
    //if( argc!=2 )
    //{
    //    help();
    //    return 1;
    //}
    //string filename = argv[1];
    //if( filename.empty() )
    //{
    //    cout << "\nDurn, couldn't read in " << argv[1] << endl;
    //    return 1;
    //}
	//filename = "test.jpg";
    //Mat image = imread( filename, 1);		//RGB
	Mat grayImage=imread(filename,0);		//Gray


	//Mat grayImage;
	//cvtColor(image,grayImage, CV_BGR2GRAY);


    help();

    const string winName = "Main";
    namedWindow( winName, WINDOW_AUTOSIZE );
    setMouseCallback( winName, on_mouse, 0 );

    gcapp.setImageAndWinName( grayImage, winName );

    gcapp.show();

	int iterCount = 0;
	int newIterCount = 0;


    for(;;)
    {
		int c = waitKey(0);
		if(c == '\x1b')
		{
			cout << "Exiting ..." << endl;
			goto exit_main;
		}else if(c == 'r')
		{
			cout << endl;
			gcapp.reset();
			gcapp.show();
		}else if(c == 's')//superpixels
		{
			cout << ">superpixels start" << endl;
	
			gcapp.superpixelSegmentation();

			gcapp.showSuperImage();

			cout << ">superpixels end" << endl;

		}else if(c == 'n')//pixel cut
		{
			cout << ">pixel cut start" << endl;
			iterCount = gcapp.getIterCount();
			cout << "<" << iterCount << "... ";
			newIterCount = gcapp.nextIter();

			gcapp.clearMarkInImage();
			gcapp.show();
			cout << iterCount << ">" << endl;

			cout << ">pixel cut end" << endl;
		}else if(c == 'm')//superpixel cut
		{
			cout << ">superpixel cut start" << endl;
			iterCount = gcapp.getIterCount();
			cout << "<" << iterCount << "... ";

			newIterCount = gcapp.nextSuperIter();

			gcapp.clearMarkInImage();
			gcapp.show();
			cout << iterCount << ">" << endl;
			cout << ">superpixel cut end" << endl;
		}


    }

exit_main:
    destroyWindow( winName );
    return 0;
}
