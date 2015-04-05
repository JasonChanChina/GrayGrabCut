#include "GCApplication.h"

const Scalar RED = Scalar(0,0,255);
const Scalar PINK = Scalar(230,130,255);
const Scalar BLUE = Scalar(255,0,0);
const Scalar LIGHTBLUE = Scalar(255,255,160);
const Scalar GREEN = Scalar(0,255,0);
const Scalar WHITE = Scalar(255,255,255);
const Scalar BLACK = Scalar(0,0,0);

const Scalar COLOR_RECT_FORE = Scalar(255,0,0);	//灰度图像只用到了第一分量，即255=白
const Scalar COLOR_RECT_BACK = Scalar(0,0,0);//灰度图像只用到了第一分量，即0=黑
const Scalar COLOR_FGD = Scalar(128,0,0);
const Scalar COLOR_BGD = Scalar(0,0,0);
const Scalar COLOR_PR_FGD = Scalar(192,0,0);
const Scalar COLOR_PR_BGD = Scalar(0,0,0);

const int KEY_CTRL = CV_EVENT_FLAG_CTRLKEY;
const int KEY_SHIFT = CV_EVENT_FLAG_SHIFTKEY;
const int KEY_ALT = CV_EVENT_FLAG_ALTKEY;

void GCApplication::reset()
{

	mask.create( image.size(), CV_8UC1);
	mask.setTo(Scalar::all(GC_PR_FGD));

	bgdPxls.clear(); fgdPxls.clear();
	prBgdPxls.clear();  prFgdPxls.clear();

	leftState = STATE_UP;
	rightState = STATE_UP;
	iterCount = 0;
	isFirst = true;
}

void GCApplication::clearMarkInImage()
{
	image = originImage.clone();
	superImage = originSuperImage.clone();
}

void GCApplication::setImageAndWinName( const Mat& _image, const string& _winName  )
{
	if( _image.empty() || _winName.empty() )
		return;
	image = _image;
	superImage = image.clone();
	originImage = image.clone();
	originSuperImage = superImage.clone();

	beta = grab.calcBeta(image);

	winName = _winName;
	superName = "Superpixel_Image";

	reset();
}

void GCApplication::show()
{
	showImage();
	showSuperImage();
}

void GCApplication::showImage()
{

	//设置背景颜色
	Mat res;
	res.create(image.size(), image.type());
	res.setTo(255);				//设置空白区域颜色


	Mat binMask;
	getBinMask( mask, binMask );
	image.copyTo( res, binMask );


	imshow( winName, res );
	//string path = "D:\\aa_";
	//path += "_a.png";
	//imwrite(path, res);


}

void GCApplication::showSuperImage()
{
	//设置背景颜色
	Mat res;
	res.create(superImage.size(), superImage.type());
	res.setTo(255);				//设置空白区域颜色

	Mat binMask;
	getBinMask( mask, binMask );
	superImage.copyTo( res, binMask );

	imshow( superName, res );
}


void GCApplication::mouseClick( int event, int x, int y, int flags, void* )
{
	// TODO add bad args check

	//Ctrl+Left : Rect->fgd_pr
	//Shift+Left : fgd_pr . Shift+Right : fgd
	//Alt+Left : bgd_pr , Alt+Right : bgd


	switch( event )
	{
	case CV_EVENT_LBUTTONDOWN: 
		{
			leftState = STATE_DOWN;
			beginPoint.x = x;
			beginPoint.y = y;
		}
		break;
	case CV_EVENT_RBUTTONDOWN:
		{
			rightState = STATE_DOWN;
			beginPoint.x = x;
			beginPoint.y = y;
		}
		break;
	case CV_EVENT_LBUTTONUP:
		{
			leftState = STATE_UP;
			endPoint.x = x;
			endPoint.y = y;
			Point temp(endPoint);
			if((flags & KEY_CTRL) != 0)		//ctrl+left
			{
				Rect rect;
				rect.x = std::min(beginPoint.x, endPoint.x);
				rect.y = std::min(beginPoint.y, endPoint.y);
				rect.width = std::abs(beginPoint.x - endPoint.x);
				rect.height = std::abs(beginPoint.y - endPoint.y);

				Size tempSize = image.size();

				rect.x = std::max(0, rect.x);
				rect.y = std::max(0, rect.y);
				rect.width = std::min(rect.width, tempSize.width-rect.x);
				rect.height = std::min(rect.height, tempSize.height-rect.y);

				rectangle( image, Point( rect.x, rect.y ), Point(rect.x + rect.width, rect.y + rect.height ), COLOR_RECT_FORE, 2);
				rectangle( superImage, Point( rect.x, rect.y ), Point(rect.x + rect.width, rect.y + rect.height ), COLOR_RECT_FORE, 2);

				rects.push_back(rect);

			}else if((flags & KEY_SHIFT) != 0)	//shift+left
			{
				prFgdPxls.push_back(temp);
				//circle( mask, temp, radius, GC_PR_FGD, thickness );
				circle( image, temp, radius, COLOR_PR_FGD, thickness );
				circle( superImage, temp, radius, COLOR_PR_FGD, thickness );
			}else if((flags & KEY_ALT) != 0)	//alt+left
			{
				fgdPxls.push_back(temp);
				//circle( mask, temp, radius, GC_FGD, thickness );
				circle( image, temp, radius, COLOR_FGD, thickness );
				circle( superImage, temp, radius, COLOR_FGD, thickness );
			}else{		//无控制键，纯鼠标
				Rect rect;
				rect.x = std::min(beginPoint.x, endPoint.x);
				rect.y = std::min(beginPoint.y, endPoint.y);
				rect.width = std::abs(beginPoint.x - endPoint.x);
				rect.height = std::abs(beginPoint.y - endPoint.y);

				Size tempSize = image.size();

				rect.x = std::max(0, rect.x);
				rect.y = std::max(0, rect.y);
				rect.width = std::min(rect.width, tempSize.width-rect.x);
				rect.height = std::min(rect.height, tempSize.height-rect.y);
				RotatedRect rRect = RotatedRect(Point(rect.x+rect.width/2, rect.y+rect.height/2), rect.size(), 0);
				ellipse(image, rRect, COLOR_RECT_FORE, 2);
				ellipse(superImage, rRect, COLOR_RECT_FORE, 2);
			}

			show();
		}
		break;
	case CV_EVENT_RBUTTONUP:
		{
			rightState = STATE_UP;
			endPoint.x = x;
			endPoint.y = y;
			Point temp(endPoint);
			if((flags & KEY_CTRL) != 0)		//ctrl+right
			{
				Rect rect;
				rect.x = std::min(beginPoint.x, endPoint.x);
				rect.y = std::min(beginPoint.y, endPoint.y);
				rect.width = std::abs(beginPoint.x - endPoint.x);
				rect.height = std::abs(beginPoint.y - endPoint.y);

				Size tempSize = image.size();

				rect.x = std::max(0, rect.x);
				rect.y = std::max(0, rect.y);
				rect.width = std::min(rect.width, tempSize.width-rect.x);
				rect.height = std::min(rect.height, tempSize.height-rect.y);

				rectangle( image, Point( rect.x, rect.y ), Point(rect.x + rect.width, rect.y + rect.height ), COLOR_RECT_BACK, 2);
				rectangle( superImage, Point( rect.x, rect.y ), Point(rect.x + rect.width, rect.y + rect.height ), COLOR_RECT_BACK, 2);

				rects.push_back(rect);
			}else if((flags & KEY_SHIFT) != 0)	//shift+right
			{
				prBgdPxls.push_back(temp);
				//circle( mask, temp, radius, GC_PR_BGD, thickness );
				circle( image, temp, radius, COLOR_PR_BGD, thickness );
				circle( superImage, temp, radius, COLOR_PR_BGD, thickness );
			}else if((flags & KEY_ALT) != 0)	//alt+right
			{
				bgdPxls.push_back(temp);
				//circle( mask, temp, radius, GC_BGD, thickness );
				circle( image, temp, radius, COLOR_BGD, thickness );
				circle( superImage, temp, radius, COLOR_BGD, thickness );
			}else{	//无控制键，纯鼠标
				Rect rect;
				rect.x = std::min(beginPoint.x, endPoint.x);
				rect.y = std::min(beginPoint.y, endPoint.y);
				rect.width = std::abs(beginPoint.x - endPoint.x);
				rect.height = std::abs(beginPoint.y - endPoint.y);

				Size tempSize = image.size();

				rect.x = std::max(0, rect.x);
				rect.y = std::max(0, rect.y);
				rect.width = std::min(rect.width, tempSize.width-rect.x);
				rect.height = std::min(rect.height, tempSize.height-rect.y);
				RotatedRect rRect = RotatedRect(Point(rect.x+rect.width/2, rect.y+rect.height/2), rect.size(), 0);
				ellipse(image, rRect, COLOR_RECT_BACK, 2);
				ellipse(superImage, rRect, COLOR_RECT_BACK, 2);
			}

			show();
		}
		break;
	case CV_EVENT_MOUSEMOVE:
		{
			if(leftState == STATE_DOWN )
			{
				middlePoint.x = x;
				middlePoint.y = y;

				Point temp(middlePoint);
				if((flags & KEY_CTRL) != 0)		//ctrl+left
				{

				}else if((flags & KEY_SHIFT) != 0)	//shift+left
				{
					prFgdPxls.push_back(temp);
					//circle( mask, temp, radius, GC_PR_FGD, thickness );
					circle( image, temp, radius, COLOR_PR_FGD, thickness );
					circle( superImage, temp, radius, COLOR_PR_FGD, thickness );
				}else if((flags & KEY_ALT) != 0)	//alt+left
				{
					fgdPxls.push_back(temp);
					//circle( mask, temp, radius, GC_FGD, thickness );
					circle( image, temp, radius, COLOR_FGD, thickness );
					circle( superImage, temp, radius, COLOR_FGD, thickness );
				}
			}else if(rightState == STATE_DOWN)
			{
				middlePoint.x = x;
				middlePoint.y = y;
				Point temp(middlePoint);
				if((flags & KEY_CTRL) != 0)		//ctrl+right
				{

				}else if((flags & KEY_SHIFT) != 0)	//shift+right
				{
					prBgdPxls.push_back(temp);
					//circle( mask, temp, radius, GC_PR_BGD, thickness );
					circle( image, temp, radius, COLOR_PR_BGD, thickness );
					circle( superImage, temp, radius, COLOR_PR_BGD, thickness );
				}else if((flags & KEY_ALT) != 0)	//alt+right
				{
					bgdPxls.push_back(temp);
					//circle( mask, temp, radius, GC_BGD, thickness );
					circle( image, temp, radius, COLOR_BGD, thickness );
					circle( superImage, temp, radius, COLOR_BGD, thickness );
				}
			}
			show();
		}
		break;
	}
}


int GCApplication::nextIter()
{
	//double begin = (double)getTickCount();
	updateMaskFromRectsAndPixels();


	grab.graygrabCut( image, mask, bgdModel, fgdModel, 1, beta);

	//double end = (double)getTickCount();
	//double time = (end-begin)/getTickFrequency();
	//cout<<"Time="<<time<<endl;

	iterCount++;

	rects.clear();
	bgdPxls.clear(); fgdPxls.clear();
	prBgdPxls.clear(); prFgdPxls.clear();


	return iterCount;
}


int GCApplication::nextSuperIter()
{
	//double begin = (double)getTickCount();
	updateMaskFromRectsAndPixels();

	mask2SPMask();

	grab.graySupergrabCut(pixelLabels, kindOfLabels, arcs, centers,  contains,  SPMask,  bgdModel, fgdModel,iterCount, beta);


	//double end = (double)getTickCount();
	//double time = (end-begin)/getTickFrequency();
	//cout<<"Time="<<time<<endl;

	iterCount++;

	rects.clear();
	bgdPxls.clear(); fgdPxls.clear();
	prBgdPxls.clear(); prFgdPxls.clear();

	SPMask2Mask();
	return iterCount;
}

int GCApplication::getIterCount() const 
{ 
	return iterCount; 
}


void GCApplication::updateMaskFromRectsAndPixels()
{
	if(isFirst) mask.setTo( GC_BGD );
	isFirst = false;

	Size size = image.size();
	for(int i = 0; i < rects.size(); i++)
	{
		(mask(rects[i])).setTo( GC_PR_FGD);
	}

	vector<Point>::const_iterator it;
	for( it = prBgdPxls.begin(); it != prBgdPxls.end(); ++it )
		circle( mask, *it, radius, GC_PR_BGD, thickness );
	for( it = prFgdPxls.begin(); it != prFgdPxls.end(); ++it )
		circle( mask, *it, radius, GC_PR_FGD, thickness );
	for( it = bgdPxls.begin(); it != bgdPxls.end(); ++it )
		circle( mask, *it, radius, GC_BGD, thickness );
	for( it = fgdPxls.begin(); it != fgdPxls.end(); ++it )
		circle( mask, *it, radius, GC_FGD, thickness );


}


void GCApplication::getBinMask(Mat& comMask, Mat& binMask )
{
	if( comMask.empty() || comMask.type()!=CV_8UC1 )
		CV_Error( CV_StsBadArg, "comMask is empty or has incorrect type (not CV_8UC1)" );
	if( binMask.empty() || binMask.rows!=comMask.rows || binMask.cols!=comMask.cols )
	{
		binMask.create( comMask.size(), CV_8UC1 );
	}
	binMask = comMask & 1;

}


void GCApplication::SPMask2Mask()
{
	Point p;
	for(int i = 0; i < contains.size(); i++)
	{
		uchar label = (uchar)SPMask[i];
		for(int j = 0; j < contains[i].size(); j++)
		{
			p.y = contains[i][j] / mask.cols;
			p.x = contains[i][j] % mask.cols;
			mask.at<uchar>(p) = label;
		}
	}

}
void GCApplication::mask2SPMask()
{
	int count = SPMask.size();

	vector<vector<int> > carts(count);
	for(int i = 0; i < count; i++)
	{
		carts[i].resize(10, 0);
	}
	//统计各个超像素块的各种标记(GC_PR_FGD, GC_PR_BGD, GC_FGD, GC_BGD)的个数
	Point p;
	int cols = mask.cols;
	for(int i = 0; i < count; i++)
	{
		for(int j = 0; j < contains[i].size(); j++)
		{
			p.y = contains[i][j] / cols;
			p.x = contains[i][j] % cols;
			uchar label = mask.at<uchar>(p);
			carts[i][label]++;
		}
	}

	for(int i = 0; i < count; i++)
	{
		int pf = carts[i][GC_PR_FGD];
		int pb = carts[i][GC_PR_BGD];
		int f = carts[i][GC_FGD];
		int b = carts[i][GC_BGD];

		if(f > 0 || b > 0)
		{
			if(f >= b) SPMask[i] = GC_FGD;
			else SPMask[i] = GC_BGD;
		}else
		{
			if(pf >= pb) SPMask[i] = GC_PR_FGD;
			else SPMask[i] = GC_PR_BGD;
		}

	}

}

void GCApplication::dyeInvalidRegion()
{
	int width = image.cols;
	int height = image.rows;
	
	Mat dyeMask;
	dyeMask.create( image.size(), CV_8UC1);
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

	Mat res;
	res.create(image.size(), image.type());
	res.setTo(255);				//设置空白区域颜色
	
	image.copyTo( res, dyeMask );		//清除圆外区域



	//设置一圈黑色圆形边界
	if(width % 2 != 0 && height % 2 != 0)
	{
		Point point((int)(width+1)/2, (int)(height+1)/2);
		int dyeRadius = point.x;
		circle( res, point, dyeRadius, 0, 1 );
	}else
	{

		Point pointLU( (int)(width-1)/2, (int)(height-1)/2  );		//left + up
		Point pointRU( (int)(width+1)/2, (int)(height-1)/2  );		//right + up
		Point pointLD( (int)(width-1)/2, (int)(height+1)/2  );		//left + down
		Point pointRD( (int)(width+1)/2, (int)(height+1)/2  );		//right + down

		int dyeRadius = pointLU.x;
		circle( res, pointLU, dyeRadius, 0, 1 );
		circle( res, pointRU, dyeRadius, 0, 1 );
		circle( res, pointLD, dyeRadius, 0, 1 );
		circle( res, pointRD, dyeRadius, 0, 1 );
	}



	image = res.clone();
	originImage = image.clone();
	superImage = image.clone();
	originSuperImage = image.clone();

}

int GCApplication::countValidSuperpixels()
{
	int width = image.cols;
	int height = image.rows;

	Point points[4];
	int dy[4] = {-1, 1, -1, 1};
	int dx[4] = {1,1,-1,-1};
	for(int i = 0; i < 4; i++)			//left + up		//right + up	//left + down	//right + down
	{
		points[i].x = (int)(width+dx[i])/2;
		points[i].y = (int)(height+dy[i])/2;
	}

	Point p;
	int whiteCount = 0;
	for(int i = 0; i < contains.size(); i++)
	{
		for(int j = 0; j < contains[i].size(); j++)
		{
			p.y = contains[i][j] / mask.cols;
			p.x = contains[i][j] % mask.cols;
			for(int k = 0; k < 4; k++)
			{
				if( (p.x - points[k].x) * (p.x - points[k].x) + (p.y - points[k].y)*(p.y - points[k].y) > (width/2)*(height/2))		//大于圆范围
				{
					whiteCount++;
					goto outRegion;		//为了跳出深层循环
				}
			}

		}
outRegion:
		int f = i+1;
	}
	return whiteCount;
}

int GCApplication::superpixelSegmentation(int _step)
{
	slico.DoSuperpixelSegmentation_ForGivenMat(image, pixelLabels, kindOfLabels, _step);
	slico.GetArcAndCenterOfSuperpixels(superImage, pixelLabels, kindOfLabels, arcs, centers, contains);
	slico.DrawContoursAroundSegments(superImage, pixelLabels, WHITE);
	SPMask.resize(kindOfLabels, GC_PR_FGD);
	originSuperImage = superImage.clone();

	int count = countValidSuperpixels();

	return kindOfLabels-count;
}



