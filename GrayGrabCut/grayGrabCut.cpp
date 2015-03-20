
#include "GrayGrabCut.h"



////////////////////////////////// pixel cut /////////////////////////////////
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
  Initialize Hist background and foreground models using kmeans algorithm.
*/
void GrayGrabCut::updateHists( const Mat& img, const Mat& mask, Histogram& bgdHist, Histogram& fgdHist )
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

void GrayGrabCut::graygrabCut( Mat& img, Mat& mask, Mat& bgdModel, Mat& fgdModel, int iterCount, double beta )
{

    if( img.empty() )
        CV_Error( CV_StsBadArg, "image is empty" );
    if( img.type() != CV_8UC1 )
        CV_Error( CV_StsBadArg, "image mush have CV_8UC1 type" );

	Histogram bgdHist(bgdModel), fgdHist(fgdModel);


    if( iterCount <= 0)
        return;


    const double gamma = 50;
    const double lambda = 9*gamma;

	
    Mat leftW, upleftW, upW, uprightW;
    calcNWeights( img, leftW, upleftW, upW, uprightW, beta, gamma );

    for( int i = 0; i < iterCount; i++ )
    {
        GCGraph<double> graph;
		updateHists( img, mask, bgdHist, fgdHist);
        constructGCGraph(img, mask, bgdHist, fgdHist, lambda, leftW, upleftW, upW, uprightW, graph );
        estimateSegmentation( graph, mask );
    }
}



/////////////////////////////////// superpixel cut /////////////////////////////////////////


void GrayGrabCut::calcSuperNWeights(vector<Vec3d>& centers, vector<vector<int> >& arcs,
									vector<vector<int> >& contains, vector<vector<double> >& nweights, int& sideCount,
									double beta, double gamma )
{
	int num = arcs.size();
	nweights.resize(num);
	for(int i = 0; i < num; i++) nweights[i].resize(num, 0);
	sideCount = 0;

	for(int i = 0; i < num; i++)
	{
		for(int j = i+1; j < num; j++)
		{
			if(arcs[i][j] == 0) continue;

			double diffG = centers[i][0] - centers[j][0];
			double diffX = centers[i][1] - centers[j][1];
			double diffY = centers[i][2] - centers[j][2];
			double distXY = std::sqrt(diffX*diffX + diffY*diffY);

			int countOfi = contains[i].size();
			int countOfj = contains[j].size();
			double countRate = (double)countOfi/countOfj;
			double countRate2 = (double)countOfj/countOfi;
			nweights[i][j] = exp(-beta * diffG*diffG) * (gamma * countRate / distXY);
			nweights[j][i] = exp(-beta * diffG*diffG) * (gamma * countRate2 / distXY);
			sideCount+=2;
		}
	}

}


void GrayGrabCut::updateSuperHists(vector<Vec3d>& centers, vector<vector<int> >& contains, vector<int>&  SPMask, Histogram& bgdHist, Histogram& fgdHist )
{
	//Mask: GC_PR_BGD  GC_PR_FGD GC_BGD GC_FGD
	//mask中非0代表要被计算

	vector<int> bgdSPMask(SPMask.size(), 0);
	editSPMask(bgdSPMask, false);
	bgdHist.createSuperHist(centers, contains, bgdSPMask);

	vector<int> fgdSPMask(SPMask.size(), 0);
	editSPMask(fgdSPMask, true);
	fgdHist.createSuperHist(centers, contains, fgdSPMask);
}


/*
  Construct GCGraph
*/
void GrayGrabCut::constructSuperGCGraph(vector<Vec3d>& centers, vector<vector<int> >& arcs, 
										vector<int>& SPMask, vector<vector<double> >& nweights, 
										Histogram& fgdHist ,Histogram& bgdHist , double lambda, 
										int sideCount, vector<int>& vtxs, GCGraph<double>& graph)
{

	int vtxCount = centers.size();
	int edgeCount = sideCount;			//正确计算了来回两条边

    graph.create(vtxCount, edgeCount);
	vtxs.resize(vtxCount, -1);

	for(int i = 0; i < vtxCount; i++)
	{
		// add node
		int vtxIdx = graph.addVtx();
		vtxs[i] = vtxIdx;
		int gray = centers[i][0];

		// set t-weights
		double fromSource, toSink;
		if(SPMask[i] == GC_PR_BGD || SPMask[i] == GC_PR_FGD)
		{
			fromSource = -log(bgdHist.probability(gray));
			toSink = -log(fgdHist.probability(gray));
		}
		else if( SPMask[i] == GC_BGD )
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

	}
    
	// set n-weights

	for(int i = 0; i < nweights.size(); i++)
	{
		for(int j = i+1; j < nweights[i].size(); j++)
		{
			if(arcs[i][j] != 0 )
				graph.addEdges(i, j, nweights[i][j], nweights[j][i]);
		}
	}
}

/*
  Estimate segmentation using MaxFlow algorithm
*/
void GrayGrabCut::estimateSuperSegmentation( GCGraph<double>& graph, vector<int>& SPMask , vector<int> vtxs)
{
    graph.maxFlow();

	for(int i = 0; i < SPMask.size(); i++)
	{
		if(SPMask[i] == GC_PR_FGD || SPMask[i] == GC_PR_BGD)
		{
			if(graph.inSourceSegment(vtxs[i]))
				SPMask[i] = GC_PR_FGD;
			else
				SPMask[i] = GC_PR_BGD;
		}
	}
}



//Mask: GC_PR_BGD  GC_PR_FGD GC_BGD GC_FGD

void GrayGrabCut::graySupergrabCut(vector<int>& pixelLabels, int& kindOfLabels, 
	vector<vector<int> >& arcs, vector<Vec3d>& centers, vector<vector<int> >& contains, vector<int>& SPMask,
	Mat& bgdModel, Mat& fgdModel, int iterCount, double _beta)
{
	Histogram bgdHist(bgdModel), fgdHist(fgdModel);

	//if(firstTime)
	//{
	//	vector<int> fmask(SPMask);
	//	vector<int> bmask(SPMask);
	//	editSPMask(fmask, true);
	//	editSPMask(bmask, false);
	//	bgdHist.createSuperHist(centers, fmask);
	//	fgdHist.createSuperHist(centers, bmask);
	//}

    if( iterCount <= 0) return;

    const double gamma = 50;
    const double lambda = 9*gamma;
    const double beta = _beta;
	
	vector<vector<double> > nweights;
	vector<int> vtxs;		//记录超像素和graph里顶点对应关系
	int sideCount = 0;
	calcSuperNWeights(centers, arcs, contains, nweights, sideCount, beta, gamma );

    for( int i = 0; i < iterCount; i++ )
    {
        GCGraph<double> graph;
		updateSuperHists(centers, contains, SPMask, bgdHist, fgdHist );
		constructSuperGCGraph(centers, arcs,  SPMask, nweights, fgdHist, bgdHist , lambda, sideCount, vtxs, graph);
        estimateSuperSegmentation( graph, SPMask, vtxs );
    }
}




/////////////////////////////////// Other /////////////////////////////////////////


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


////////////////////////////////////// private ///////////////////////////////////////

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


void GrayGrabCut::editSPMask(vector<int>& SPMask, bool isFgd)
{
	int value1 = GC_BGD;
	int value2 = GC_PR_BGD;
	if(!isFgd)
	{
		value1 = GC_FGD;
		value2 = GC_PR_FGD;
	}

	for(int i = 0; i < SPMask.size(); i++)
	{
		int value = SPMask[i];
		if(value == value1 || value == value2)
		{
			SPMask[i] = 0;
		}else{
			SPMask[i] = 5;
		}
	}

}