/*--------------------------------------------------------
 * xsubgraph.h
 * Interface of xsubgraph.cc
 * Random extraction of a (possibly) connected subgraph
 * See: argraph.h
 *
 * Author: P. Foggia
 --------------------------------------------------------*/


#ifndef XSUBGRAPH_H
#define XSUBGRAPH_H


#include "BALL/STRUCTURE/VFLIB/argraph.h"

namespace BALL {
namespace VFLib {

Graph* ExtractSubgraph(Graph *g, int nodes, bool connected=true);

}
}
#endif
