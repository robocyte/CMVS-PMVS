#include <algorithm>
#include <iostream>

#include "graclus.h"

using namespace CMVS;

/*************************************************************************
* This function reads the spd matrix
**************************************************************************/
void Cgraclus::initGraph(GraphType& graph)
{
    graph.gdata  = graph.rdata = nullptr;

    graph.nvtxs  = graph.nedges = -1;
    graph.mincut = graph.minvol = -1;

    graph.xadj      = graph.vwgt = graph.adjncy = graph.adjwgt = nullptr;
    graph.adjwgtsum = nullptr;
    graph.label     = nullptr;
    graph.cmap      = nullptr;

    graph.where  = graph.pwgts = nullptr;
    graph.id     = graph.ed = nullptr;
    graph.bndptr = graph.bndind = nullptr;
    graph.rinfo  = nullptr;
    graph.vrinfo = nullptr;
    graph.nrinfo = nullptr;

    graph.ncon   = -1;
    graph.nvwgt  = nullptr;
    graph.npwgts = nullptr;

    graph.vsize  = nullptr;

    graph.coarser = graph.finer = nullptr;
}

//----------------------------------------------------------------------
// cutType
// 0: NCUT, 1: RASSO
void Cgraclus::runE(std::vector<idxtype>& xadj, std::vector<idxtype>& adjncy,
                    std::vector<idxtype>& adjwgt, const int nparts, const int cutType,
                    std::vector<idxtype>& part)
{
    GraphType graph;
    initGraph(graph);

    graph.ncon = 1;

    graph.xadj   = &xadj[0];
    graph.adjncy = &adjncy[0];
    graph.vwgt   = nullptr;
    graph.adjwgt = &adjwgt[0];

    graph.nvtxs  = (int)xadj.size() - 1;
    graph.nedges = (int)adjncy.size();

    const int wgtflag = 1;
    runSub(graph, nparts, cutType, wgtflag, part);
}

int Cgraclus::mylog2(int a)
{
    int i;
    for (i = 1 ; a > 1; i++, a = a>>1);

    return i-1;
}

void Cgraclus::runSub(GraphType& graph, int nparts, int cutType, int wgtflag, std::vector<idxtype>& part)
{
    const int levels = std::max((graph.nvtxs)/(40*mylog2(nparts)), 20*(nparts));
    part.resize(graph.nvtxs);

    int options[11];    options[0] = 0;
    int numflag = 0;
    int chain_length = 0;  int edgecut;

    MLKKM_PartGraphKway(&graph.nvtxs, graph.xadj, graph.adjncy, graph.vwgt, graph.adjwgt, 
                        &wgtflag, &numflag, &nparts, &chain_length, options, &edgecut, &part[0], levels);

    float lbvec[MAXNCON];
    ComputePartitionBalance(&graph, nparts, &part[0], lbvec);

    float result;

    if (cutType == 0)
    {
        result = ComputeNCut(&graph, &part[0], nparts);
        printf("\nNormalized-Cut... \n   Cut value: %7f, Balance: \n", result);
    } else
    {
        result = ComputeRAsso(&graph, &part[0], nparts);
        printf("\nRatio Association...  \n  Association value: %7f, Balance: \n", result);
    }
}
