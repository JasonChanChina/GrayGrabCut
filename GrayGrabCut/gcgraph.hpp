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
        int parent;	//存储路径(父节点)
        int first;	//第一个延生出去的边（以此点为起点的边）
        int ts;		//此点的流量
        int dist;	//距离root的节点数量distance
        TWeight weight;		//  =（Source-Sink） 净权  ： >0 代表是Source   <0 代表是Sink
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
    edges.reserve( edgeCount + 2 );
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

    TWeight dw = vtcs[i].weight;
    if( dw > 0 )
        sourceW += dw;	//sourceW 是正数
    else
        sinkW -= dw;	//sinkW 是正数
    flow += (sourceW < sinkW) ? sourceW : sinkW;		//累加最小的T-link
    vtcs[i].weight = sourceW - sinkW;
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
        v->ts = 0;	//各个点初始流量设为0
        if( v->weight != 0 )
        {
            last = last->next = v;	//串联 上一个节点和当前节点
            v->dist = 1;	//初始距离=1，即一个节点
            v->parent = TERMINAL;
            v->t = v->weight < 0;	//标记0/1，source=0 or sink=1
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
                    u = vtxPtr+edgePtr[ei].dst;		//目的节点，即边i->j 的目的点j
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
                if( e0 > 0 )	//找到一个不在一个阵营的点，就结束
                    break;
            }
            // exclude the vertex from the active list
            first = first->next;	
            v->next = 0;
        }

        if( e0 <= 0 )
            break;

        // find the minimum edge weight along the path
        minWeight = edgePtr[e0].weight;
        assert( minWeight > 0 );
        // k = 1: source tree, k = 0: destination tree
        for( int k = 1; k >= 0; k-- )
        {
            for( v = vtxPtr+edgePtr[e0^k].dst;; v = vtxPtr+edgePtr[ei].dst )
            {
                if( (ei = v->parent) < 0 )
                    break;
                weight = edgePtr[ei^k].weight;
                minWeight = MIN(minWeight, weight);
                assert( minWeight > 0 );
            }
            weight = fabs(v->weight);
            minWeight = MIN(minWeight, weight);
            assert( minWeight > 0 );
        }

        // modify weights of the edges along the path and collect orphans
        edgePtr[e0].weight -= minWeight;			
        edgePtr[e0^1].weight += minWeight;				//若k=2n,则k^1=2n+1,  若k=2n+1,则k^=2n, 即找出2n和2n+1的边(即找出对应的反向边)
        flow += minWeight;

        // k = 1: source tree, k = 0: destination tree
        for( int k = 1; k >= 0; k-- )
        {
            for( v = vtxPtr+edgePtr[e0^k].dst;; v = vtxPtr+edgePtr[ei].dst )
            {
                if( (ei = v->parent) < 0 )
                    break;
                edgePtr[ei^(k^1)].weight += minWeight;
                if( (edgePtr[ei^k].weight -= minWeight) == 0 )
                {
                    orphans.push_back(v);
                    v->parent = ORPHAN;
                }
            }

            v->weight = v->weight + minWeight*(1-k*2);
            if( v->weight == 0 )
            {
               orphans.push_back(v);
               v->parent = ORPHAN;
            }
        }

        // restore the search trees by finding new parents for the orphans
        curr_ts++;
        while( !orphans.empty() )
        {
            Vtx* v2 = orphans.back();
            orphans.pop_back();

            int d, minDist = INT_MAX;
            e0 = 0;
            vt = v2->t;

            for( ei = v2->first; ei != 0; ei = edgePtr[ei].next )
            {
                if( edgePtr[ei^(vt^1)].weight == 0 )
                    continue;
                u = vtxPtr+edgePtr[ei].dst;
                if( u->t != vt || u->parent == 0 )
                    continue;
                // compute the distance to the tree root
                for( d = 0;; )
                {
                    if( u->ts == curr_ts )
                    {
                        d += u->dist;
                        break;
                    }
                    ej = u->parent;
                    d++;
                    if( ej < 0 )
                    {
                        if( ej == ORPHAN )
                            d = INT_MAX-1;
                        else
                        {
                            u->ts = curr_ts;
                            u->dist = 1;
                        }
                        break;
                    }
                    u = vtxPtr+edgePtr[ej].dst;
                }

                // update the distance
                if( ++d < INT_MAX )
                {
                    if( d < minDist )
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

            if( (v2->parent = e0) > 0 )
            {
                v2->ts = curr_ts;
                v2->dist = minDist;
                continue;
            }

            /* no parent is found */
            v2->ts = 0;
            for( ei = v2->first; ei != 0; ei = edgePtr[ei].next )
            {
                u = vtxPtr+edgePtr[ei].dst;
                ej = u->parent;
                if( u->t != vt || !ej )
                    continue;
                if( edgePtr[ei^(vt^1)].weight && !u->next )
                {
                    u->next = nilNode;
                    last = last->next = u;
                }
                if( ej > 0 && vtxPtr+edgePtr[ej].dst == v2 )
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

