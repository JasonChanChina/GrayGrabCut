

#include "SLICO.h"

const int dx4[4] = {-1,  0,  1,  0};
const int dy4[4] = { 0, -1,  0,  1};
const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};



SLICO::SLICO()
{
	m_gvec = NULL;

	m_gvecvec = NULL;
}

SLICO::~SLICO()
{
	if(m_gvec) delete [] m_gvec;


	if(m_gvecvec)
	{
		for( int d = 0; d < m_depth; d++ ) delete [] m_gvecvec[d];
		delete [] m_gvecvec;
	}
	
}



void SLICO::DoRGBtoLABConversion(
	const unsigned int*&		ubuff,
	double*&					gvec)
{
	int sz = m_width*m_height;
	gvec = new double[sz];

	for( int j = 0; j < sz; j++ )
	{
		int g = (ubuff[j]      ) & 0xFF;
		gvec[j] = g;
	}
}



//绘制边界
void SLICO::DrawContoursAroundSegments(Mat& mat, vector<int>& pixelLabels, Scalar color)
{

	int lineThiner = 1;   // change to 2 or 3 for thinner lines

	int channels = mat.channels();
	int nRows = mat.rows;
	int nCols = mat.cols * channels;
	int height = mat.rows;
	int width = mat.cols * channels;

	//假设最多三通道
	Vec3d cvColor;
	cvColor[0] = color[0];
	cvColor[1] = color[1];
	cvColor[2] = color[2];

	
	int sz = mat.rows * mat.cols;
	vector<bool> istaken(sz, false);


	if(mat.isContinuous())
	{
		nCols *= nRows;
		nRows = 1;

		uchar* data = mat.ptr<uchar>(0);
		for(int i = 0; i < nCols; i+=channels)
		{
			int r = i / width;
			int c = i % width;

			int np(0);
			for(int k = 0; k < 8; k++)
			{
				int newr = r + dy8[k];
				int newc = c + dx8[k]*channels;
				if( (newc >= 0 && newc < width) && (newr >= 0 && newr < height) )
				{
					int newIndex = newr * width + newc;
					if(false == istaken[newIndex/channels])
					{
						if(pixelLabels[i/channels] != pixelLabels[newIndex/channels])
						{
							np++;
						}
					}

				}
			}
			if(np > lineThiner) //change to 2 or 3 for thinner lines
			{
				for(int chan = 0; chan < channels; chan++)
				{
					*(data+i+chan) = cvColor[chan];
				}
				istaken[i/channels] = true;
			}

		}

	}else
	{
		uchar* data;
		for(int r = 0; r < nRows; r++)
		{
			data = mat.ptr<uchar>(r);
			for(int c = 0; c < nCols; c+=channels)
			{
				int index = r * width + c;

				int np(0);
				for(int k = 0; k < 8; k++)
				{
					int newr = r + dy8[k];
					int newc = c + dx8[k]*channels;
					if( (newc >= 0 && newc < width) && (newr >= 0 && newr < height) )
					{
						int newIndex = newr * width + newc;
						if(false == istaken[newIndex/channels])
						{
							if(pixelLabels[index/channels] != pixelLabels[newIndex/channels])
							{
								np++;
							}
						}

					}
				}
				if(np > lineThiner) //change to 2 or 3 for thinner lines
				{
					for(int chan = 0; chan < channels; chan++)
					{
						*(data+c+chan) = cvColor[chan];
					}
					istaken[index/channels] = true;
				}
			}
		}
	}


	/*
	if(mat.isContinuous())
	{
		nCols *= nRows;
		nRows = 1;

		uchar* data = mat.ptr<uchar>(0);
		for(int i = 0; i < nCols; i+=channels)
		{
			int r = i / width;
			int c = i % width;
			bool isSidePixel = false;
			for(int k = 0; k < 8; k++)
			{
				int newr = r + dy8[k];
				int newc = c + dx8[k]*channels;
				if( (newc >= 0 && newc < width) && (newr >= 0 && newr < height) )
				{
					int newIndex = newr * width + newc;
					if(pixelLabels[i/channels] != pixelLabels[newIndex/channels])
					{
						isSidePixel = true;
						break;
					}
				}
			}
			if(isSidePixel)
			{
				*(data+i) = cvColor[0];
				*(data+i+1) = cvColor[1];
				*(data+i+2) = cvColor[2];
			}

		}


	}else
	{
		uchar* data;
		for(int r = 0; r < nRows; r++)
		{
			data = mat.ptr<uchar>(r);
			for(int c = 0; c < nCols; c+=channels)
			{
				int index = r * width + c;
				bool isSidePixel = false;
				for(int k = 0; k < 8; k++)
				{
					int newr = r + dy8[k];
					int newc = c + dx8[k]*channels;
					if( (newc >= 0 && newc < width) && (newr >= 0 && newr < height) )
					{
						int newIndex = newr * width + newc;
						if(pixelLabels[index/channels] != pixelLabels[newIndex/channels])
						{
							isSidePixel = true;
							break;
						}
					}
				}
				if(isSidePixel)
				{
					*(data+c) = cvColor[0];
					*(data+c+1) = cvColor[1];
					*(data+c+2) = cvColor[2];
				}
			}
		}
	}
	*/




}


void SLICO::DetectLabEdges(
	const double*				gvec,
	const int&					width,
	const int&					height,
	vector<double>&				edges)
{
	int sz = width*height;

	edges.resize(sz,0);
	for( int j = 1; j < height-1; j++ )
	{
		for( int k = 1; k < width-1; k++ )
		{
			int i = j*width+k;

			double dx = (gvec[i-1]-gvec[i+1])*(gvec[i-1]-gvec[i+1]);
			double dy = (gvec[i-width]-gvec[i+width])*(gvec[i-width]-gvec[i+width]);

			//edges[i] = (sqrt(dx) + sqrt(dy));
			edges[i] = (dx + dy);
		}
	}
}


void SLICO::PerturbSeeds(
	vector<double>&				kseedsg,
	vector<double>&				kseedsx,
	vector<double>&				kseedsy,
	const vector<double>&		edges)
{
	const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};
	
	int numseeds = kseedsg.size();

	for( int n = 0; n < numseeds; n++ )
	{
		int ox = kseedsx[n];//original x
		int oy = kseedsy[n];//original y
		int oind = oy*m_width + ox;

		int storeind = oind;
		for( int i = 0; i < 8; i++ )
		{
			int nx = ox+dx8[i];//new x
			int ny = oy+dy8[i];//new y

			if( nx >= 0 && nx < m_width && ny >= 0 && ny < m_height)
			{
				int nind = ny*m_width + nx;
				if( edges[nind] < edges[storeind])
				{
					storeind = nind;
				}
			}
		}
		if(storeind != oind)
		{
			kseedsx[n] = storeind%m_width;
			kseedsy[n] = storeind/m_width;
			kseedsg[n] = m_gvec[storeind];
		}
	}
}


void SLICO::GetLABXYSeeds_ForGivenStepSize(
	vector<double>&				kseedsg,
	vector<double>&				kseedsx,
	vector<double>&				kseedsy,
	const int&					STEP,
	const bool&					perturbseeds,
	const vector<double>&		edgemag)
{
	int numseeds(0);
	int n(0);

	//int xstrips = m_width/STEP;
	//int ystrips = m_height/STEP;
	int xstrips = (0.5+double(m_width)/double(STEP));
	int ystrips = (0.5+double(m_height)/double(STEP));

	int xerr = m_width  - STEP*xstrips;
	int yerr = m_height - STEP*ystrips;

	double xerrperstrip = double(xerr)/double(xstrips);
	double yerrperstrip = double(yerr)/double(ystrips);

	int xoff = STEP/2;
	int yoff = STEP/2;
	//-------------------------
	numseeds = xstrips*ystrips;
	//-------------------------
	kseedsg.resize(numseeds);
	kseedsx.resize(numseeds);
	kseedsy.resize(numseeds);

	for( int y = 0; y < ystrips; y++ )
	{
		int ye = y*yerrperstrip;
		for( int x = 0; x < xstrips; x++ )
		{
			int xe = x*xerrperstrip;
			int i = (y*STEP+yoff+ye)*m_width + (x*STEP+xoff+xe);
			
			kseedsg[n] = m_gvec[i];
			kseedsx[n] = (x*STEP+xoff+xe);
			kseedsy[n] = (y*STEP+yoff+ye);
			n++;
		}
	}

	
	if(perturbseeds)
	{
		PerturbSeeds(kseedsg, kseedsx, kseedsy, edgemag);
	}
}


void SLICO::GetLABXYSeeds_ForGivenK(
	vector<double>&				kseedsg,
	vector<double>&				kseedsx,
	vector<double>&				kseedsy,
	const int&					K,
	const bool&					perturbseeds,
	const vector<double>&		edgemag)
{
	int sz = m_width*m_height;
	double step = sqrt(double(sz)/double(K));
	int T = step;
	int xoff = step/2;
	int yoff = step/2;
	
	int n(0);int r(0);
	for( int y = 0; y < m_height; y++ )
	{
		int Y = y*step + yoff;
		if( Y > m_height-1 ) break;

		for( int x = 0; x < m_width; x++ )
		{
			//int X = x*step + xoff;//square grid
			int X = x*step + (xoff<<(r&0x1));//hex grid
			if(X > m_width-1) break;

			int i = Y*m_width + X;

			//_ASSERT(n < K);
			
			//kseedsl[n] = m_lvec[i];
			//kseedsa[n] = m_avec[i];
			//kseedsb[n] = m_bvec[i];
			//kseedsx[n] = X;
			//kseedsy[n] = Y;
			kseedsg.push_back(m_gvec[i]);
			kseedsx.push_back(X);
			kseedsy.push_back(Y);
			n++;
		}
		r++;
	}

	if(perturbseeds)
	{
		PerturbSeeds(kseedsg, kseedsx, kseedsy, edgemag);
	}
}



void SLICO::PerformSuperpixelSegmentation_VariableSandM(
	vector<double>&				kseedsg,
	vector<double>&				kseedsx,
	vector<double>&				kseedsy,
	int*						pixelLabels,
	const int&					STEP,
	const int&					NUMITR)
{
	int sz = m_width*m_height;
	const int numk = kseedsg.size();
	//double cumerr(99999.9);
	int numitr(0);

	//----------------
	int offset = STEP;
	m_size = STEP;
	if(STEP < 10) offset = STEP*1.5;
	//----------------

	vector<double> sigmag(numk, 0);
	vector<double> sigmax(numk, 0);
	vector<double> sigmay(numk, 0);
	vector<int> clustersize(numk, 0);
	vector<double> inv(numk, 0);//to store 1/clustersize[k] values
	vector<double> distxy(sz, DBL_MAX);
	vector<double> distg(sz, DBL_MAX);
	vector<double> distvec(sz, DBL_MAX);
	vector<double> maxg(numk, 10*10);//THIS IS THE VARIABLE VALUE OF M, just start with 10
	vector<double> maxxy(numk, STEP*STEP);//THIS IS THE VARIABLE VALUE OF M, just start with 10

	double invxywt = 1.0/(STEP*STEP);//NOTE: this is different from how usual SLIC/LKM works

	while( numitr < NUMITR )
	{
		//------
		//cumerr = 0;
		numitr++;
		//------

		distvec.assign(sz, DBL_MAX);
		for( int n = 0; n < numk; n++ )
		{
			int y1 = max((double)0,			kseedsy[n]-offset);
			int y2 = min((double)m_height,	kseedsy[n]+offset);
			int x1 = max((double)0,			kseedsx[n]-offset);
			int x2 = min((double)m_width,	kseedsx[n]+offset);

			for( int y = y1; y < y2; y++ )
			{
				for( int x = x1; x < x2; x++ )
				{
					int i = y*m_width + x;
					_ASSERT( y < m_height && x < m_width && y >= 0 && x >= 0 );

					double g = m_gvec[i];
					distg[i] = (g-kseedsg[n])*(g-kseedsg[n]);

					distxy[i] =		(x - kseedsx[n])*(x - kseedsx[n]) +
									(y - kseedsy[n])*(y - kseedsy[n]);

					//------------------------------------------------------------------------
					double dist = distg[i]/maxg[n] + distxy[i]*invxywt;//only varying m, prettier superpixels
					//double dist = distlab[i]/maxlab[n] + distxy[i]/maxxy[n];//varying both m and S
					//------------------------------------------------------------------------
					
					if( dist < distvec[i] )
					{
						distvec[i] = dist;
						pixelLabels[i]  = n;
					}
				}
			}
		}
		//-----------------------------------------------------------------
		// Assign the max color distance for a cluster
		//-----------------------------------------------------------------
		if(0 == numitr)
		{
			maxg.assign(numk,1);
			maxxy.assign(numk,1);
		}
		{for( int i = 0; i < sz; i++ )
		{
			if(maxg[pixelLabels[i]] < distg[i]) maxg[pixelLabels[i]] = distg[i];
			if(maxxy[pixelLabels[i]] < distxy[i]) maxxy[pixelLabels[i]] = distxy[i];
		}}
		//-----------------------------------------------------------------
		// Recalculate the centroid and store in the seed values
		//-----------------------------------------------------------------
		sigmag.assign(numk, 0);
		sigmax.assign(numk, 0);
		sigmay.assign(numk, 0);
		clustersize.assign(numk, 0);

		for( int j = 0; j < sz; j++ )
		{
			int temp = pixelLabels[j];
			_ASSERT(pixelLabels[j] >= 0);
			sigmag[pixelLabels[j]] += m_gvec[j];
			sigmax[pixelLabels[j]] += (j%m_width);
			sigmay[pixelLabels[j]] += (j/m_width);

			clustersize[pixelLabels[j]]++;
		}

		{for( int k = 0; k < numk; k++ )
		{
			//_ASSERT(clustersize[k] > 0);
			if( clustersize[k] <= 0 ) clustersize[k] = 1;
			inv[k] = 1.0/double(clustersize[k]);//computing inverse now to multiply, than divide later
		}}
		
		{for( int k = 0; k < numk; k++ )
		{
			kseedsg[k] = sigmag[k]*inv[k];
			kseedsx[k] = sigmax[k]*inv[k];
			kseedsy[k] = sigmay[k]*inv[k];
		}}
	}
}



void SLICO::EnforceLabelConnectivity(
	const int*					labels,//input labels that need to be corrected to remove stray labels
	const int&					width,
	const int&					height,
	int*						nlabels,//new labels
	int&						kindOfLabels,//the number of labels changes in the end if segments are removed
	const int&					K) //the number of superpixels desired by the user
{
//	const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
//	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};

	const int dx4[4] = {-1,  0,  1,  0};
	const int dy4[4] = { 0, -1,  0,  1};

	const int sz = width*height;
	const int SUPSZ = sz/K;
	//nlabels.resize(sz, -1);
	for( int i = 0; i < sz; i++ ) nlabels[i] = -1;
	int label(0);
	int* xvec = new int[sz];
	int* yvec = new int[sz];
	int oindex(0);
	int adjlabel(0);//adjacent label
	for( int j = 0; j < height; j++ )
	{
		for( int k = 0; k < width; k++ )
		{
			if( 0 > nlabels[oindex] )
			{
				nlabels[oindex] = label;
				//--------------------
				// Start a new segment
				//--------------------
				xvec[0] = k;
				yvec[0] = j;
				//-------------------------------------------------------
				// Quickly find an adjacent label for use later if needed
				//-------------------------------------------------------
				{for( int n = 0; n < 4; n++ )
				{
					int x = xvec[0] + dx4[n];
					int y = yvec[0] + dy4[n];
					if( (x >= 0 && x < width) && (y >= 0 && y < height) )
					{
						int nindex = y*width + x;
						if(nlabels[nindex] >= 0) adjlabel = nlabels[nindex];
					}
				}}

				int count(1);
				for( int c = 0; c < count; c++ )
				{
					for( int n = 0; n < 4; n++ )
					{
						int x = xvec[c] + dx4[n];
						int y = yvec[c] + dy4[n];

						if( (x >= 0 && x < width) && (y >= 0 && y < height) )
						{
							int nindex = y*width + x;

							if( 0 > nlabels[nindex] && labels[oindex] == labels[nindex] )
							{
								xvec[count] = x;
								yvec[count] = y;
								nlabels[nindex] = label;
								count++;
							}
						}

					}
				}
				//-------------------------------------------------------
				// If segment size is less then a limit, assign an
				// adjacent label found before, and decrement label count.
				//-------------------------------------------------------
				if(count <= SUPSZ >> 2)
				{
					for( int c = 0; c < count; c++ )
					{
						int ind = yvec[c]*width+xvec[c];
						nlabels[ind] = adjlabel;
					}
					label--;
				}
				label++;
			}
			oindex++;
		}
	}
	kindOfLabels = label;

	if(xvec) delete [] xvec;
	if(yvec) delete [] yvec;
}


void SLICO::PerformSLICO_ForGivenStepSize(
	const unsigned int*			ubuff,
	const int					width,
	const int					height,
	int*						pixelLabels,
	int&						kindOfLabels,
	const int&					STEP,
	const double&				m)
{
	vector<double> kseedsg(0);
	//vector<double> kseedsl(0);
	//vector<double> kseedsa(0);
	//vector<double> kseedsb(0);
	vector<double> kseedsx(0);
	vector<double> kseedsy(0);

	//--------------------------------------------------
	m_width  = width;
	m_height = height;
	int sz = m_width*m_height;
	//pixelLabels.resize( sz, -1 );
	//--------------------------------------------------
	//pixelLabels = new int[sz];
	for( int s = 0; s < sz; s++ ) pixelLabels[s] = -1;
	//--------------------------------------------------
	DoRGBtoLABConversion(ubuff, m_gvec);
	//--------------------------------------------------

	bool perturbseeds(true);
	vector<double> edgemag(0);
	if(perturbseeds) DetectLabEdges(m_gvec, m_width, m_height, edgemag);
	GetLABXYSeeds_ForGivenStepSize(kseedsg, kseedsx, kseedsy, STEP, perturbseeds, edgemag);

	PerformSuperpixelSegmentation_VariableSandM(kseedsg,kseedsx,kseedsy,pixelLabels,STEP,10);
	kindOfLabels = kseedsg.size();

	int* nlabels = new int[sz];
	EnforceLabelConnectivity(pixelLabels, m_width, m_height, nlabels, kindOfLabels, double(sz)/double(STEP*STEP));
	{for(int i = 0; i < sz; i++ ) pixelLabels[i] = nlabels[i];}
	if(nlabels) delete [] nlabels;
}


void SLICO::PerformSLICO_ForGivenK(
	const unsigned int*			ubuff,
	const int					width,
	const int					height,
	int*						pixelLabels,
	int&						kindOfLabels,
	const int&					K,//required number of superpixels
	const double&				m)//weight given to spatial distance
{
	vector<double> kseedsg(0);
	vector<double> kseedsx(0);
	vector<double> kseedsy(0);

	//--------------------------------------------------
	m_width  = width;
	m_height = height;
	int sz = m_width*m_height;
	//--------------------------------------------------
	//if(0 == pixelLabels) pixelLabels = new int[sz];
	for( int s = 0; s < sz; s++ ) pixelLabels[s] = -1;
	//--------------------------------------------------
	if(1)//LAB
	{
		DoRGBtoLABConversion(ubuff, m_gvec);
	}
	else//RGB
	{
		m_gvec = new double[sz];
		for( int i = 0; i < sz; i++ )
		{
			m_gvec[i] = ubuff[i]       & 0xff;
		}
	}
	//--------------------------------------------------

	bool perturbseeds(true);
	vector<double> edgemag(0);
	if(perturbseeds) DetectLabEdges(m_gvec, m_width, m_height, edgemag);
	GetLABXYSeeds_ForGivenK(kseedsg, kseedsx, kseedsy, K, perturbseeds, edgemag);

	int STEP = sqrt(double(sz)/double(K)) + 2.0;//adding a small value in the even the STEP size is too small.
	//PerformSuperpixelSLIC(kseedsl, kseedsa, kseedsb, kseedsx, kseedsy, pixelLabels, STEP, edgemag, m);
	PerformSuperpixelSegmentation_VariableSandM(kseedsg,kseedsx,kseedsy,pixelLabels,STEP,10);
	kindOfLabels = kseedsg.size();

	int* nlabels = new int[sz];
	EnforceLabelConnectivity(pixelLabels, m_width, m_height, nlabels, kindOfLabels, K);
	{for(int i = 0; i < sz; i++ ) pixelLabels[i] = nlabels[i];}
	if(nlabels) delete [] nlabels;
}



void SLICO::GetArcAndCenterOfSuperpixels( const Mat& img, vector<int>& pixelLabels, int& kindOfLabels, vector<vector<int> >& arcs, vector<Vec3d>& centers, vector<vector<int> >& contains)	
{
	int width = img.cols;
	int height = img.rows;




	//Centers
	vector<int> labelCount(kindOfLabels,0);

	vector<double> sumG(kindOfLabels, 0);
	vector<double> sumX(kindOfLabels, 0);
	vector<double> sumY(kindOfLabels, 0);


	contains.resize(kindOfLabels);

	Point p;
	int index = 0;
    for( p.y = 0; p.y < img.rows; p.y++ )
    {
        for( p.x = 0; p.x < img.cols; p.x++ )
        {
			int index = p.y * width + p.x;
			int label = pixelLabels[index];
			uint gray = img.at<uchar>(p);

			labelCount[label]++;

			sumG[label] += gray;
			sumY[label] += p.y;
			sumX[label] += p.x;

			contains[label].push_back(index);
		}
	}

	Vec3d v3;
	v3[0] = 0, v3[1] = 0, v3[2] = 0;
	if(centers.size() != kindOfLabels) centers.resize(kindOfLabels, v3);
	for(int i = 0; i < kindOfLabels; i++)
	{
		sumG[i] /= labelCount[i];
		sumX[i] /= labelCount[i];
		sumY[i] /= labelCount[i];

		int centerG = ((uint)(sumG[i]+0.5)) & 0xFF;

		centers[i][0] = centerG;
		centers[i][1] = sumX[i];
		centers[i][2] = sumY[i];
	}

	//Arcs

	arcs.resize(kindOfLabels);			//arcs是 kindOfLabels * kindOfLabels 的正方形表格
	for(int i = 0; i < kindOfLabels; i++) 
		arcs[i].resize(kindOfLabels, 0);

	for(int i = 0; i < kindOfLabels; i++)
	{
		for(int j = 0; j < kindOfLabels; j++)
		{
			if(i == j) continue;
			if(abs(centers[j][1]-centers[i][1]) < m_size && abs(centers[j][2]-centers[i][2]) < m_size)
			{
				arcs[i][j] = 1;
				arcs[j][i] = 1;
			}

		}
	}


		////Arcs
	//vector<hash_set<int> > sets(kindOfLabels);
	//for(int r = 0; r < height; r+=2)
	//{
	//	for(int c = 0; c < width; c+=2)
	//	{
	//		int index = r * width + c;

	//		for(int k = 0; k < 4; k++)
	//		{
	//			int newr = r + dy4[k];
	//			int newc = c + dx4[k];
	//			if((height > newr && newr >= 0) && (width > newc && newc >= 0))
	//			{
	//				int newIndex = newr * width + newc;
	//				if(pixelLabels[index] != pixelLabels[newIndex])
	//				{
	//					int lmin = std::min(pixelLabels[index], pixelLabels[newIndex]);
	//					int lmax = std::max(pixelLabels[index], pixelLabels[newIndex]);
	//					sets[lmin].insert(lmax);
	//					sets[lmax].insert(lmin);
	//				}
	//			}
	//		}
	//	}
	//}

	//if(arcs.size() != kindOfLabels) arcs.resize(kindOfLabels);
	//hash_set<int>::iterator iter;
	//for(int i = 0; i < kindOfLabels; i++)
	//{
	//	arcs[i].clear();
	//	for(iter = sets[i].begin(); iter != sets[i].end(); iter++)
	//	{
	//		arcs[i].push_back(*iter);
	//	}
	//}


}

//
//void SLICO::DrawAverageColor(Mat &mat, vector<int>& pixelLabels, vector<Vec6d>& centers)
//{
//		
//	for(int r = 0; r < mat.rows; r++)
//	{
//		for(int c = 0; c < mat.cols; c++)
//		{
//			int pinIndex = r * mat.cols + c;
//			int spIndex = pixelLabels[pinIndex];
//			Vec6d v6 = centers[spIndex];
//			Vec3d v3;
//			v3[0] = v6[0];
//			v3[1] = v6[1];
//			v3[2] = v6[2];
//			mat.at<Vec3b>(r,c) = v3;
//		}
//	}
//}


//=================================================================
//对Mat图像进行分割
//=================================================================
void SLICO::DoSuperpixelSegmentation_ForGivenMat( const Mat& img, vector<int>& pixelLabels, int& kindOfLabels)							
{
	int width = img.cols;
	int height = img.rows;
	int size = width * height;

	uint* ubuff = (uint*)malloc(sizeof(uint)*size);

	Point p;
	int index = 0;
    for( p.y = 0; p.y < img.rows; p.y++ )
    {
        for( p.x = 0; p.x < img.cols; p.x++ )
        {
			uchar gray = img.at<uchar>(p);
			int index = p.y * img.cols + p.x;
			ubuff[index] = gray;
		}
	}


	int* _pixelLabels = (int*)malloc(sizeof(int)*size);
	kindOfLabels = 0;
	int STEP = 10;  // 10 =1536个 // 16=600个
	double compactness = 20.0;

	PerformSLICO_ForGivenStepSize(ubuff, width, height, _pixelLabels, kindOfLabels, STEP, compactness);
	delete ubuff;

	pixelLabels.resize(size);
	for(int i = 0; i < size; i++)
		pixelLabels[i] = _pixelLabels[i];
	delete _pixelLabels;
}


