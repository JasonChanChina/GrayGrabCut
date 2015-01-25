
#include "mex.h"

#include "stdafx.h"

#include "GrayGrabCut.h"

//
//void mat2MxArray(Mat matGray)
//{
//	/////////// 方法1
//	int rows = matGray.rows;  
//    int cols = matGray.cols;  
//
//    mxArray* mxGray = mxCreateDoubleMatrix(rows, cols, mxREAL); 
//
//    double *data;  
//    data = mxGetPr(mxGray); 
//
//    for (int i = 0; i < rows; i++)  
//        for (int j = 0; j < cols; j++)  
//            *(data + i + j * rows) = (double)matGray.at<uchar>(i, j);  
//
//
//
//	///////////// 方法2
//	mxArray * mxGray;  
//	if (!matGray.empty())  
//	{  
//		matGray = matGray.t();  
//		int rows = matGray.rows;
//		int cols = matGray.cols;
//
//		mxGray = mxCreateNumericMatrix(cols,rows,mxSINGLE_CLASS, mxREAL);  
//
//		double *data;  
//		data = mxGetPr(mxGray); 
//
//		memcpy(data, matGray.data, mxGetNumberOfElements(mxGray)*sizeof(float));  
//	}  
//
//}
//
//void mxArray2Mat(mxArray mxGray)
//{
//
//
//	///////////// 方法2
//	mat* matGray;  
//  
//
//	memcpy(matGray.data, mxGetPr(mxGray), mxGetNumberOfElements(mxGray)*sizeof(float));  
//	
//}

void mexErrMsgIdAndTxt(std::string msg, std::string info)
{
	std::cout<<msg<<" "<<info<<std::endl;
}



mxClassID getClassId(int depth)
{
	mxClassID targetClassId;
	switch(depth)
	{
	case CV_64F:
		targetClassId = mxDOUBLE_CLASS;
		break;
	case CV_32F:
		targetClassId = mxSINGLE_CLASS;
		break;
	case CV_8S:
		targetClassId = mxINT8_CLASS;
		break;
	case CV_8U:
		targetClassId = mxUINT8_CLASS;
		break;
	case CV_16S:
		targetClassId = mxINT16_CLASS;
		break;
	case CV_16U:
		targetClassId = mxUINT16_CLASS;
		break;
	case CV_32S:
		targetClassId = mxINT32_CLASS;
		break;
	default:
		targetClassId = mxDOUBLE_CLASS;
		break;
	}
	return targetClassId;
}


int getDepth(mxClassID classid)
{
	int targetDepth;
	switch(classid)
	{
	case mxDOUBLE_CLASS:
		targetDepth = CV_64F;
		break;
	case mxSINGLE_CLASS:
		targetDepth = CV_32F;
		break;
	case mxINT8_CLASS:
		targetDepth = CV_8S;
		break;
	case mxUINT8_CLASS:
		targetDepth = CV_8U;
		break;
	case mxINT16_CLASS:
		targetDepth = CV_16S;
		break;
	case mxUINT16_CLASS:
		targetDepth = CV_16U;
		break;
	case mxINT32_CLASS:
		targetDepth = CV_32S;
		break;
	case mxUINT32_CLASS:
		targetDepth = CV_32S;
		break;
	case mxLOGICAL_CLASS:
		targetDepth = CV_8U;
		break;
	default:
		targetDepth = CV_64F;
		break;
	}
	return targetDepth;
}


mxArray* mat2MxArray(cv::Mat& mat, mxClassID classid=mxUNKNOWN_CLASS, bool transpose=true)
{
	mxArray* mx;
	if(mat.empty())
	{
		mx = mxCreateNumericArray(0,0,mxDOUBLE_CLASS, mxREAL);
		if(!mx)
			mexErrMsgIdAndTxt("Error", "Allocation error");
		return mx;
	}

	cv::Mat input = (mat.dims == 2 && transpose) ? mat.t() : mat;
	// Create a new mxArray.
	const int nchannels = input.channels();
	const int* dims_ = input.size;
	std::vector<mwSize> d(dims_, dims_ + input.dims);
	d.push_back(nchannels);

	mxClassID targetClassId = getClassId(input.depth());


	classid = (classid == mxUNKNOWN_CLASS) ? targetClassId : classid;
	std::swap(d[0], d[1]);
	if (classid == mxLOGICAL_CLASS)
	{
		// OpenCV's logical true is any nonzero while matlab's true is 1.
		cv::compare(input, 0, input, cv::CMP_NE);
		input.setTo(1, input);
		mx = mxCreateLogicalArray(d.size(), &d[0]);
	} else {
		mx = mxCreateNumericArray(d.size(), &d[0], classid, mxREAL);
	}
	if (!mx){
		mexErrMsgIdAndTxt("mexopencv:error", "Allocation error");
	}
	// Copy each channel.
	std::vector<cv::Mat> channels;
	cv::split(input, channels);
	std::vector<mwSize> si(d.size(), 0); // subscript index.

	int targetDepth = getDepth(classid);

	int type = CV_MAKETYPE(targetDepth, 1); // destination type.
	for (int i = 0; i < nchannels; ++i)
	{
		si[si.size() - 1] = i; // last dim is a channel index.

		mwIndex subsSi = mxCalcSingleSubscript(mx, si.size(), &si[0]);

		void *ptr = reinterpret_cast<void*>(
			reinterpret_cast<size_t>(mxGetData(mx)) +
			mxGetElementSize(mx) * subsSi);
		cv::Mat m(input.dims, dims_, type, ptr);
		channels[i].convertTo(m, type); // Write to mxArray through m.
	}
	return mx;
}


cv::Mat mxArray2Mat(mxArray* mx, int depth=CV_USRTYPE1, bool transpose=true)
{
	// Create cv::Mat object.
	std::vector<int> d(mxGetDimensions(mx), mxGetDimensions(mx)+mxGetNumberOfDimensions(mx));
	int ndims = (d.size()>2) ? d.size()-1 : d.size();
	int nchannels = (d.size()>2) ? *(d.end()-1) : 1;


	int targetDepth = getDepth(mxGetClassID(mx));


	depth = (depth==CV_USRTYPE1) ? targetDepth : depth;
	std::swap(d[0], d[1]);
	cv::Mat mat(ndims, &d[0], CV_MAKETYPE(depth, nchannels));
	// Copy each channel.
	std::vector<cv::Mat> channels(nchannels);
	std::vector<mwSize> si(d.size(), 0); // subscript index
	int type = CV_MAKETYPE(targetDepth, 1); // Source type
	for (int i = 0; i<nchannels; ++i)
	{
		si[d.size()-1] = i;

		mwIndex subsSi = mxCalcSingleSubscript(mx, si.size(), &si[0]);

		void *pd = reinterpret_cast<void*>(
			reinterpret_cast<size_t>(mxGetData(mx))+
			mxGetElementSize(mx)*subsSi);
		cv::Mat m(ndims, &d[0], type, pd);
		// Read from mxArray through m
		m.convertTo(channels[i], CV_MAKETYPE(depth, 1));
	}
	cv::merge(channels, mat);
	return (mat.dims==2 && transpose) ? cv::Mat(mat.t()) : mat;
}


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])  
{  
	vector<mxArray*> rhs(prhs,prhs+nrhs);

	int iterCount = 10;
	int mode = GC_INIT_WITH_MASK;
	cv::Mat mask;
	mask = mxArray2Mat(rhs[1], CV_8U);
	Rect rect;
	cv::Mat bgdModel, fgdModel;
	cv::Mat img(mxArray2Mat(rhs[0]));
	GrayGrabCut gc;
	gc.graygrabCut(img, mask, rect, bgdModel, fgdModel, iterCount, mode);

	plhs[0] = mat2MxArray(mask);
	nlhs = 1;

}
