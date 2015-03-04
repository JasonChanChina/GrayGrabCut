
#include "GrayGrabCut.h"


double GrayGrabCut::calcBeta( const Mat& img )
{
    double beta = 0;
    for( int y = 0; y < img.rows; y++ )
    {
        for( int x = 0; x < img.cols; x++ )
        {
            double color = img.at<uchar>(y,x);
            if( x>0 ) // left
            {
                double diff = color - (double)img.at<uchar>(y,x-1);
                beta += diff*diff;
            }
            if( y>0 && x>0 ) // upleft
            {
                double diff = color - (double)img.at<uchar>(y-1,x-1);
                beta += diff*diff;
            }
            if( y>0 ) // up
            {
                double diff = color - (double)img.at<uchar>(y-1,x);
                beta += diff*diff;
            }
            if( y>0 && x<img.cols-1) // upright
            {
                double diff = color - (double)img.at<uchar>(y-1,x+1);
                beta += diff*diff;
            }
        }
    }
    if( beta <= std::numeric_limits<double>::epsilon() )
        beta = 0;
    else
        beta = 1.f / (2 * beta/(4*img.cols*img.rows - 3*img.cols - 3*img.rows + 2) );

    return beta;
}

/*
  Calculate weights of noterminal vertices of graph.
  beta and gamma - parameters of GrabCut algorithm.
 */
void GrayGrabCut::calcNWeights( const Mat& img, Mat& leftW, Mat& upleftW, Mat& upW, Mat& uprightW, double beta, double gamma )
{
    const double gammaDivSqrt2 = gamma / std::sqrt(2.0f);
    leftW.create( img.rows, img.cols, CV_64FC1 );
    upleftW.create( img.rows, img.cols, CV_64FC1 );
    upW.create( img.rows, img.cols, CV_64FC1 );
    uprightW.create( img.rows, img.cols, CV_64FC1 );
    for( int y = 0; y < img.rows; y++ )
    {
        for( int x = 0; x < img.cols; x++ )
        {
            double color = img.at<uchar>(y,x);
            if( x-1>=0 ) // left
            {
                double diff = color - (double)img.at<uchar>(y,x-1);
                leftW.at<double>(y,x) = gamma * exp(-beta*diff*diff);
            }
            else
                leftW.at<double>(y,x) = 0;
            if( x-1>=0 && y-1>=0 ) // upleft
            {
                double diff = color - (double)img.at<uchar>(y-1,x-1);
                upleftW.at<double>(y,x) = gammaDivSqrt2 * exp(-beta*diff*diff);
            }
            else
                upleftW.at<double>(y,x) = 0;
            if( y-1>=0 ) // up
            {
                double diff = color - (double)img.at<uchar>(y-1,x);
                upW.at<double>(y,x) = gamma * exp(-beta*diff*diff);
            }
            else
                upW.at<double>(y,x) = 0;
            if( x+1<img.cols && y-1>=0 ) // upright
            {
                double diff = color - (double)img.at<uchar>(y-1,x+1);
                uprightW.at<double>(y,x) = gammaDivSqrt2 * exp(-beta*diff*diff);
            }
            else
                uprightW.at<double>(y,x) = 0;
        }
    }
}

/*
  Check size, type and element values of mask matrix.
 */
void GrayGrabCut::checkMask( const Mat& img, const Mat& mask )
{
    if( mask.empty() )
        CV_Error( CV_StsBadArg, "mask is empty" );
    if( mask.type() != CV_8UC1 )
        CV_Error( CV_StsBadArg, "mask must have CV_8UC1 type" );
    if( mask.cols != img.cols || mask.rows != img.rows )
        CV_Error( CV_StsBadArg, "mask must have as many rows and cols as img" );
    for( int y = 0; y < mask.rows; y++ )
    {
        for( int x = 0; x < mask.cols; x++ )
        {
            uchar val = mask.at<uchar>(y,x);
            if( val!=GC_BGD && val!=GC_FGD && val!=GC_PR_BGD && val!=GC_PR_FGD )
                CV_Error( CV_StsBadArg, "mask element value must be equel"
                    "GC_BGD or GC_FGD or GC_PR_BGD or GC_PR_FGD" );
        }
    }
}

/*
  Initialize mask using rectangular.
*/
void GrayGrabCut::initMaskWithRect( Mat& mask, Size imgSize, Rect rect )
{
    mask.create( imgSize, CV_8UC1 );
    mask.setTo( GC_BGD );

    rect.x = max(0, rect.x);
    rect.y = max(0, rect.y);
    rect.width = min(rect.width, imgSize.width-rect.x);
    rect.height = min(rect.height, imgSize.height-rect.y);

    (mask(rect)).setTo( Scalar(GC_PR_FGD) );
}


void GrayGrabCut::editMask(Mat& mask, bool isFgd)
{
	int value1 = GC_BGD;
	int value2 = GC_PR_BGD;
	if(!isFgd)
	{
		value1 = GC_FGD;
		value2 = GC_PR_FGD;
	}

    for( int y = 0; y < mask.rows; y++ )
    {
        for( int x = 0; x < mask.cols; x++ )
        {
            uchar value = mask.at<uchar>(y,x);
			if(value == value1 || value == value2)
			{
				mask.at<uchar>(y,x) = 0;
			}else
			{
				mask.at<uchar>(y,x) = 5;
			}
			
		}
	}

}

/*
  Initialize Hist background and foreground models using kmeans algorithm.
*/
void GrayGrabCut::initHists( const Mat& img, const Mat& mask, Histogram& bgdHist, Histogram& fgdHist )
{
	//mask中非0代表要被计算

	Mat bgdMask = mask.clone();
	editMask(bgdMask, false);
	bgdHist.createHist(img, bgdMask);


	Mat fgdMask = mask.clone();
	editMask(fgdMask, true);
	fgdHist.createHist(img, fgdMask);
}

/*
  Construct GCGraph
*/
void GrayGrabCut::constructGCGraph( const Mat& img, const Mat& mask, const Histogram& bgdHist, const Histogram& fgdHist, double lambda,
                       const Mat& leftW, const Mat& upleftW, const Mat& upW, const Mat& uprightW,
                       GCGraph<double>& graph )
{
    int vtxCount = img.cols*img.rows,
        edgeCount = 2*(4*img.cols*img.rows - 3*(img.cols + img.rows) + 2);	//8方向  未考虑S / T 点
		//			2（有向图）* （ 4（左上4方向）*cols*rows - 3（倍）*（cols(向上) + rows（向左）） + 2（补上重复减掉的））
    graph.create(vtxCount, edgeCount);
    Point p;
    for( p.y = 0; p.y < img.rows; p.y++ )
    {
        for( p.x = 0; p.x < img.cols; p.x++)
        {
            // add node
            int vtxIdx = graph.addVtx();
            uchar color = img.at<uchar>(p);

            // set t-weights
            double fromSource, toSink;
            if( mask.at<uchar>(p) == GC_PR_BGD || mask.at<uchar>(p) == GC_PR_FGD )
            {
				double kk = bgdHist.probability(color);
				double kk2 = fgdHist.probability(color);
				fromSource = -log( bgdHist.probability(color) );
				toSink = -log( fgdHist.probability(color) );
            }
            else if( mask.at<uchar>(p) == GC_BGD )
            {
                fromSource = 0;
                toSink = lambda;
            }
            else // GC_FGD
            {
                fromSource = lambda;
                toSink = 0;
            }
            graph.addTermWeights( vtxIdx, fromSource, toSink );

            // set n-weights
            if( p.x>0 )
            {
                double w = leftW.at<double>(p);
                graph.addEdges( vtxIdx, vtxIdx-1, w, w );
            }
            if( p.x>0 && p.y>0 )
            {
                double w = upleftW.at<double>(p);
                graph.addEdges( vtxIdx, vtxIdx-img.cols-1, w, w );
            }
            if( p.y>0 )
            {
                double w = upW.at<double>(p);
                graph.addEdges( vtxIdx, vtxIdx-img.cols, w, w );
            }
            if( p.x<img.cols-1 && p.y>0 )
            {
                double w = uprightW.at<double>(p);
                graph.addEdges( vtxIdx, vtxIdx-img.cols+1, w, w );
            }
        }
    }
}

/*
  Estimate segmentation using MaxFlow algorithm
*/
void GrayGrabCut::estimateSegmentation( GCGraph<double>& graph, Mat& mask )
{
    graph.maxFlow();
    Point p;
    for( p.y = 0; p.y < mask.rows; p.y++ )
    {
        for( p.x = 0; p.x < mask.cols; p.x++ )
        {
            if( mask.at<uchar>(p) == GC_PR_BGD || mask.at<uchar>(p) == GC_PR_FGD )
            {
                if( graph.inSourceSegment( p.y*mask.cols+p.x /*vertex index*/ ) )
                    mask.at<uchar>(p) = GC_PR_FGD;
                else
                    mask.at<uchar>(p) = GC_PR_BGD;
            }
        }
    }
}

void GrayGrabCut::graygrabCut( InputArray _img, InputOutputArray _mask, Rect rect,
                  InputOutputArray _bgdModel, InputOutputArray _fgdModel,
                  int iterCount, int mode )
{
    Mat img = _img.getMat();
    Mat& mask = _mask.getMatRef();
	Mat bgdModel = _bgdModel.getMat();
	Mat fgdModel = _fgdModel.getMat();

    if( img.empty() )
        CV_Error( CV_StsBadArg, "image is empty" );
    if( img.type() != CV_8UC1 )
        CV_Error( CV_StsBadArg, "image mush have CV_8UC1 type" );

	Histogram bgdHist(bgdModel), fgdHist(fgdModel);

    if( mode == GC_INIT_WITH_RECT || mode == GC_INIT_WITH_MASK )
    {
        if( mode == GC_INIT_WITH_RECT )
            initMaskWithRect( mask, img.size(), rect );
        checkMask( img, mask );
    }

    if( iterCount <= 0)
        return;

    if( mode == GC_EVAL )
        checkMask( img, mask );

    const double gamma = 50;
    const double lambda = 9*gamma;
    const double beta = calcBeta( img );
	
    Mat leftW, upleftW, upW, uprightW;
    calcNWeights( img, leftW, upleftW, upW, uprightW, beta, gamma );

    for( int i = 0; i < iterCount; i++ )
    {
        GCGraph<double> graph;
		initHists( img, mask, bgdHist, fgdHist);
        constructGCGraph(img, mask, bgdHist, fgdHist, lambda, leftW, upleftW, upW, uprightW, graph );
        estimateSegmentation( graph, mask );
    }
}
