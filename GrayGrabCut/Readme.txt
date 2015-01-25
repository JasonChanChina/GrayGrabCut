matlab+CPP混合编译说明
==========================
文件结构说明：
0. common include
	.H    : stdafx.h
1. GrayGrabCut (grabcut算法)
	.CPP : GrayGrabCut.cpp	GMM.cpp		
	.H     : GrayGrabCut.h	GMM.h		
	.HPP : gcgraph.hpp
2. mexopencv   (CPP与Mex数据类型转换框架)
	.CPP : MxArray.cpp
	.H     : 
	.HPP : MxArray.hpp		mexopencv.hpp
3. Cpp2MexInterface (混合编程接口函数)
	.CPP : mexGrayGrabCut.cpp


=========================
=========================
用matlab R2013b x64编译
=========================
环境：Win7 x64 en + matlab R2013b X64  + OpenCV 2.4.8 x64 + mexopencv(2015.01.22) + VS2012 x64
=================
编译步骤：
1. 把matlab的当前目录设置到本工程CPP文件目录
2. mex -setup 设置编译器为VS编译器
3. 运行一下命令编译：
	mex mexGrayGrabCut.cpp MxArray.cpp GrayGrabCut.cpp GMM.cpp  -outdir ./ -O -DNDEBUG -I.\ -IC:\opencv\build\include -largeArrayDims -LC:\opencv\build\x64\vc11\lib -lopencv_core248 -lopencv_highgui248 -lopencv_video248 -lopencv_imgproc248
--------------------------------
命令格式解析(help mex)：
	mex 所有相关CPP文件 .....
	-outdir + 目录	：	输出目录(参数与值之间有空格)
	-O	:	优化编译
	-DNDEBUG	：	不知道
	-largeArrayDims	:	使用X64的数值空间
	-I+目录	：	include目录(参数与值之间无空格)
	-L+目录	：	Lib目录(参数与值之间无空格)
	-l+库	：	链接库(参数与值之间无空格)(这里用到opencv的库)

=========================