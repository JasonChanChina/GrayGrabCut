matlab+CPP��ϱ���˵��
==========================
�ļ��ṹ˵����
0. common include
	.H    : stdafx.h
1. GrayGrabCut (grabcut�㷨)
	.CPP : GrayGrabCut.cpp	GMM.cpp		
	.H     : GrayGrabCut.h	GMM.h		
	.HPP : gcgraph.hpp
2. mexopencv   (CPP��Mex��������ת�����)
	.CPP : MxArray.cpp
	.H     : 
	.HPP : MxArray.hpp		mexopencv.hpp
3. Cpp2MexInterface (��ϱ�̽ӿں���)
	.CPP : mexGrayGrabCut.cpp


=========================
=========================
��matlab R2013b x64����
=========================
������Win7 x64 en + matlab R2013b X64  + OpenCV 2.4.8 x64 + mexopencv(2015.01.22) + VS2012 x64
=================
���벽�裺
1. ��matlab�ĵ�ǰĿ¼���õ�������CPP�ļ�Ŀ¼
2. mex -setup ���ñ�����ΪVS������
3. ����һ��������룺
	mex mexGrayGrabCut.cpp MxArray.cpp GrayGrabCut.cpp GMM.cpp  -outdir ./ -O -DNDEBUG -I.\ -IC:\opencv\build\include -largeArrayDims -LC:\opencv\build\x64\vc11\lib -lopencv_core248 -lopencv_highgui248 -lopencv_video248 -lopencv_imgproc248
--------------------------------
�����ʽ����(help mex)��
	mex �������CPP�ļ� .....
	-outdir + Ŀ¼	��	���Ŀ¼(������ֵ֮���пո�)
	-O	:	�Ż�����
	-DNDEBUG	��	��֪��
	-largeArrayDims	:	ʹ��X64����ֵ�ռ�
	-I+Ŀ¼	��	includeĿ¼(������ֵ֮���޿ո�)
	-L+Ŀ¼	��	LibĿ¼(������ֵ֮���޿ո�)
	-l+��	��	���ӿ�(������ֵ֮���޿ո�)(�����õ�opencv�Ŀ�)

=========================