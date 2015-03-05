/*M///////////////////////////////////////////////////////////////////////////////////////
//
//  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
//
//  By downloading, copying, installing or using the software you agree to this license.
//  If you do not agree to this license, do not download, install,
//  copy or use the software.
//
//
//                        Intel License Agreement
//                For Open Source Computer Vision Library
//
// Copyright (C) 2000, Intel Corporation, all rights reserved.
// Third party copyrights are property of their respective owners.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
//   * Redistribution's of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//
//   * Redistribution's in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
//   * The name of Intel Corporation may not be used to endorse or promote products
//     derived from this software without specific prior written permission.
//
// This software is provided by the copyright holders and contributors "as is" and
// any express or implied warranties, including, but not limited to, the implied
// warranties of merchantability and fitness for a particular purpose are disclaimed.
// In no event shall the Intel Corporation or contributors be liable for any direct,
// indirect, incidental, special, exemplary, or consequential damages
// (including, but not limited to, procurement of substitute goods or services;
// loss of use, data, or profits; or business interruption) however caused
// and on any theory of liability, whether in contract, strict liability,
// or tort (including negligence or otherwise) arising in any way out of
// the use of this software, even if advised of the possibility of such damage.
//
//M*/


// 来自于论文：An Experimental Comparison of Min-Cut/Max-Flow Algorithms for Energy Minimization in Vision
// 其中的第3节

#pragma once

template <class TWeight> class GCGraph
{
public:
    GCGraph();
    GCGraph( unsigned int vtxCount, unsigned int edgeCount );
    ~GCGraph();
    void create( unsigned int vtxCount, unsigned int edgeCount );
    int addVtx();				//添加一个点
    void addEdges( int i, int j, TWeight w/*i-j权重*/, TWeight revw/*j->i权重*/ );	//添加边（权重）
    void addTermWeights( int i, TWeight sourceW, TWeight sinkW );	//添加特殊边（权重）
    TWeight maxFlow();
    bool inSourceSegment( int i );		//判断某点是否属于Source这边
private:
    class Vtx
    {
    public:
        Vtx *next; // initialized and used in maxFlow() only
        int parent;	//存储路径，边
        int first;	//第一个延生出去的边（以此点为起点的边）
        int ts;		//to source，即标志是否直接/间接连接到S/T，0=未连接，1,2,3,4...=代表不同的路径p的标志
        int dist;	//距离root的节点数量distance
        TWeight vWeight;		//  =（Source-Sink） 净权  ： >0 代表是Source   <0 代表是Sink
        uchar t;	//=0 属于source,, !=0 属于sink
    };
    class Edge
    {
    public:
        int dst;		//目的点，边的终点（边的起点未保存，因为程序中是通过点找边，无须保存起点）
        int next;		//下一个兄弟边（起点一样的边）
        TWeight weight;
    };

    std::vector<Vtx> vtcs;		//所有点集合
    std::vector<Edge> edges;	//所有边集合
    TWeight flow;				//保存最大流值
};

template <class TWeight>
GCGraph<TWeight>::GCGraph()
{
    flow = 0;
}
template <class TWeight>
GCGraph<TWeight>::GCGraph( unsigned int vtxCount, unsigned int edgeCount )
{
    create( vtxCount, edgeCount );
}
template <class TWeight>
GCGraph<TWeight>::~GCGraph()
{
}
template <class TWeight>
void GCGraph<TWeight>::create( unsigned int vtxCount, unsigned int edgeCount )
{
    vtcs.reserve( vtxCount );
    edges.reserve( edgeCount + 2 );		//多出某两点到S/T的两条边
    flow = 0;
}

template <class TWeight>
int GCGraph<TWeight>::addVtx()
{
    Vtx v;
    memset( &v, 0, sizeof(Vtx));
    vtcs.push_back(v);
    return (int)vtcs.size() - 1;
}

template <class TWeight>
void GCGraph<TWeight>::addEdges( int i, int j, TWeight w, TWeight revw )
{
    CV_Assert( i>=0 && i<(int)vtcs.size() );
    CV_Assert( j>=0 && j<(int)vtcs.size() );
    CV_Assert( w>=0 && revw>=0 );
    CV_Assert( i != j );

    if( !edges.size() )
        edges.resize( 2 );

    Edge fromI, toI;
    fromI.dst = j;
    fromI.next = vtcs[i].first;
    fromI.weight = w;
    vtcs[i].first = (int)edges.size();
    edges.push_back( fromI );

    toI.dst = i;
    toI.next = vtcs[j].first;
    toI.weight = revw;
    vtcs[j].first = (int)edges.size();
    edges.push_back( toI );
}

template <class TWeight>
void GCGraph<TWeight>::addTermWeights( int i, TWeight sourceW, TWeight sinkW )
{
    CV_Assert( i>=0 && i<(int)vtcs.size() );

    TWeight dw = vtcs[i].vWeight;
    if( dw > 0 )
        sourceW += dw;	//sourceW 是正数
    else
        sinkW -= dw;	//sinkW 是正数
    flow += (sourceW < sinkW) ? sourceW : sinkW;		//累加最小的T-link
    vtcs[i].vWeight = sourceW - sinkW;
}

template <class TWeight>
TWeight GCGraph<TWeight>::maxFlow()
{
    const int TERMINAL = -1, /* 末端 */  ORPHAN = -2;   /* 孤立 */
    Vtx stub, *nilNode = &stub, *first = nilNode, *last = nilNode;	//nilNode指向空点stub
    int curr_ts = 0;
    stub.next = nilNode;
    Vtx *vtxPtr = &vtcs[0];
    Edge *edgePtr = &edges[0];

    std::vector<Vtx*> orphans;

    // initialize the active queue and the graph vertices	//串联节点
    for( int i = 0; i < (int)vtcs.size(); i++ )
    {
        Vtx* v = vtxPtr + i;
        v->ts = 0;	//设置初始值=0
        if( v->vWeight != 0 )
        {
            last = last->next = v;	//串联 上一个节点和当前节点
            v->dist = 1;	//初始距离=1，即一个节点
            v->parent = TERMINAL;
            v->t = v->vWeight < 0;	//标记0/1，source=0 or sink=1
        }
        else
            v->parent = 0;		//未被串联的点
    }
    first = first->next;	//指向第一个有效点
    last->next = nilNode;	//把最后一个点指向空点，即空点在最后，终止标志
    nilNode->next = 0;		//空点不指向任何东西

    // run the search-path -> augment-graph -> restore-trees loop
    for(;;)
    {
        Vtx* v, *u;
        int e0 = -1, ei = 0, ej = 0;
        TWeight minWeight, weight;
        uchar vt;	//即v->t

		// 1. 找出增广路径p，即从图上找出一条跨阵营的流，无需通过所有点，这样的流的起点和终点再分别连接S和T，则成为了一条S->T或者T->S的有效流（还不是最大流）
        // grow S & T search trees, find an edge connecting them
        while( first != nilNode )		//遍历串联上的所有 点及其邻居边(即生长grow)
        {
            v = first;
            if( v->parent )     // v->parent != 0 ，只处理被串联的点
            {
                vt = v->t;
                for( ei = v->first; ei != 0; ei = edgePtr[ei].next )	//遍历v点的所有邻边
                {
                    if( edgePtr[ei^vt].weight == 0 )	//若ei是从偏向source的点(v->t==0)伸出来的边，则就是当前边，若是从偏向sink的点伸出来的边，则求反向边
                        continue;
                    u = vtxPtr+edgePtr[ei].dst;		//目的节点，即v的邻居点，即边i->j 的目的点j
                    if( !u->parent )		//u->parent == 0 ，即此点未被处理过
                    {
                        u->t = vt;			//设置其属于v->t的阵营
                        u->parent = ei ^ 1;	//若k=2n,则k^1=2n+1,  若k=2n+1,则k^=2n, 即找出2n和2n+1的边(即找出对应的反向边)  
						//if(ei % 2 == 0) u->parent = ei+1;
						//else u->parent = ei-1;

                        u->ts = v->ts;
                        u->dist = v->dist + 1;	//路径长度+1
                        if( !u->next )	// u->next == null,即没有Next，表示是最后一个
                        {
                            u->next = nilNode;
                            last = last->next = u;
                        }
                        continue;	//继续查找邻居边
                    }

                    if( u->t != vt )	//不在一个阵营
                    {
                        e0 = ei ^ vt;	//存储到e0  //若ei是从偏向source的点(v->t==0)伸出来的边，则就是当前边，若是从偏向sink的点伸出来的边，则求反向边
                        break;		//找到S->T的路径，结束
                    }

                    if( u->dist > v->dist+1 && u->ts <= v->ts )		//找到一个已被处理过的更差的点，更新路径信息和长度
                    {
                        // reassign the parent
                        u->parent = ei ^ 1;			//若k=2n,则k^1=2n+1,  若k=2n+1,则k^=2n, 即找出2n和2n+1的边(即找出对应的反向边)
                        u->ts = v->ts;
                        u->dist = v->dist + 1;
                    }
                }
                if( e0 > 0 )	//找到一个跨越阵营的边，就结束
                    break;
            }
            // exclude the vertex from the active list
            first = first->next;	
            v->next = 0;
        }

        if( e0 <= 0 )
            break;

		// 2. 找出路径p上的限制容量
        // find the minimum edge weight along the path
        minWeight = edgePtr[e0].weight;			//设置初始最小权重为跨阵营边的权重
        assert( minWeight > 0 );

		//找出从跨阵营边开始，往两边走的所有边权重最小值
        // k = 1: source tree, k = 0: destination tree
        for( int k = 1; k >= 0; k-- )		//e0^k ： 若是source tree，则是从右往左，若是destination tree，则是从左往右
        {
            for( v = vtxPtr+edgePtr[e0^k].dst;; v = vtxPtr+edgePtr[ei].dst )
            {
                if( (ei = v->parent) < 0 )
                    break;
                weight = edgePtr[ei^k].weight;
                minWeight = MIN(minWeight, weight);
                assert( minWeight > 0 );
            }
            weight = fabs(v->vWeight);	//考虑S/T点到v点的边权(t-link)
            minWeight = MIN(minWeight, weight);
            assert( minWeight > 0 );
        }

        // modify weights of the edges along the path and collect orphans
        edgePtr[e0].weight -= minWeight;			
        edgePtr[e0^1].weight += minWeight;				//对应反向边 //若k=2n,则k^1=2n+1,  若k=2n+1,则k^=2n, 即找出2n和2n+1的边(即找出对应的反向边)
        flow += minWeight;

		// 3. 用限制容量修改路径p后得到残留网络，并记录边上流量为0的中立点(vWeight==0,不偏向S/T)
        // k = 1: source tree, k = 0: destination tree
        for( int k = 1; k >= 0; k-- )
        {
            for( v = vtxPtr+edgePtr[e0^k].dst;; v = vtxPtr+edgePtr[ei].dst )
            {
                if( (ei = v->parent) < 0 )
                    break;
                edgePtr[ei^(k^1)].weight += minWeight;
                if( (edgePtr[ei^k].weight -= minWeight) == 0 )	//若边上流量为0，则此点为orphans（中立点(vWeight==0,不偏向S/T)）
                {
                    orphans.push_back(v);
                    v->parent = ORPHAN;
                }
            }

			//考虑S/T点到v点的边权(t-link) 减去 限制容量
            v->vWeight = v->vWeight + minWeight*(1-k*2);
            if( v->vWeight == 0 )		
            {
               orphans.push_back(v);
               v->parent = ORPHAN;
            }
        }

        // restore the search trees by finding new parents for the orphans
        curr_ts++;	//=1,2,3,4,5... 代表不同的路径p的标志
        while( !orphans.empty() )
        {
            Vtx* v2 = orphans.back();
            orphans.pop_back();

            int d, minDist = INT_MAX;
            e0 = 0;
            vt = v2->t;

			//遍历中立点(vWeight==0,不偏向S/T)的所有邻边，找出同阵营的非中立点(vWeight==0,不偏向S/T)邻居
			//（且其起源于S/T，保证不是从中立点(vWeight==0,不偏向S/T)起源过来的），若找到，就以新v->parent连接保留在阵营这边，
			//若没找到，它和其孩子点就设置为中立点(vWeight==0,不偏向S/T)，继续存入中立点(vWeight==0,不偏向S/T)集合
            for( ei = v2->first; ei != 0; ei = edgePtr[ei].next )			
            {
                if( edgePtr[ei^(vt^1)].weight == 0 )		//跳过权重为u->v的权重<=0的点
                    continue;
                u = vtxPtr+edgePtr[ei].dst;					//v2的相邻点
                if( u->t != vt || u->parent == 0 )			//跳过非同阵营点和非路径p上的点(中立点(vWeight==0,不偏向S/T))
                    continue;

				//找出其是否起源于S/T
                // compute the distance to the tree root
                for( d = 0;; )
                {
                    if( u->ts == curr_ts )		//u和当前能连接上S/T的路径标志相同，即确认其可以连接上S/T
                    {
                        d += u->dist;
                        break;
                    }
                    ej = u->parent;
                    d++;
                    if( ej < 0 )		//边缘点
                    {
                        if( ej == ORPHAN )		//若是中立点，即无边
                            d = INT_MAX-1;		//设置为距离无限大
                        else                    //在找路径p时未访问的点
                        {               
                            u->ts = curr_ts;	//设置到当前路径中
                            u->dist = 1;		//所以设置其dist=1,因为他可以直接连接S/T
                        }
                        break;
                    }
                    u = vtxPtr+edgePtr[ej].dst;
                }

                // update the distance
                if( ++d < INT_MAX )			//判断找出的距离是合理的距离
                {
                    if( d < minDist )		//比之前存储的距离还小，则更新
                    {
                        minDist = d;
                        e0 = ei;
                    }
                    for( u = vtxPtr+edgePtr[ei].dst; u->ts != curr_ts; u = vtxPtr+edgePtr[u->parent].dst )
                    {
                        u->ts = curr_ts;
                        u->dist = --d;
                    }
                }
            }

			//设置v2到S/T路径最短的邻居边
            if( (v2->parent = e0) > 0 )
            {
                v2->ts = curr_ts;
                v2->dist = minDist;
                continue;
            }

            /* no parent is found */
			// 若未找到可以连接到S/T，则把它极其同阵营关联点都加入到孤立集合里
            v2->ts = 0;
            for( ei = v2->first; ei != 0; ei = edgePtr[ei].next )
            {
                u = vtxPtr+edgePtr[ei].dst;
                ej = u->parent;
                if( u->t != vt || !ej )		//非同阵营，或边无效(已经是边缘点了)，就跳过
                    continue;
                if( edgePtr[ei^(vt^1)].weight && !u->next )	//若是非中立点，且非末尾点，设置为未处理点，并移入末尾
                {
                    u->next = nilNode;
                    last = last->next = u;
                }
                if( ej > 0 && vtxPtr+edgePtr[ej].dst == v2 )	//
                {
                    orphans.push_back(u);
                    u->parent = ORPHAN;
                }
            }
        }
    }
    return flow;
}

template <class TWeight>
bool GCGraph<TWeight>::inSourceSegment( int i )
{
    CV_Assert( i>=0 && i<(int)vtcs.size() );
    return vtcs[i].t == 0;
};

