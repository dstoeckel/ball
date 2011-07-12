/*----------------------------------------------------
 * gene.h
 * Interface of gene.cc
 * Random generation of isomorphic ARGraphs
 * See: argraph.h
 *
 * Author: P. Foggia
 * $Id: gene.h,v 1.2 1998/12/08 13:31:01 foggia Exp $
 ----------------------------------------------------*/

/*-----------------------------------------------------------------
 * REVISION HISTORY
 *   $Log: gene.h,v $
 *   Revision 1.2  1998/12/08 13:31:01  foggia
 *   Minor changes
 *
 *---------------------------------------------------------------*/

#ifndef GENE_H

#include "BALL/STRUCTURE/VFLIB/argraph.h"

namespace BALL {
namespace VFLib {

void Generate(int nodes, int edges, Graph **g1, Graph **g2, 
              bool connected=true);

}
}
#endif
