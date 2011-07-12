/*----------------------------------------------------
 * sortnodes.h
 * Header of sortnodes.cc
 *--------------------------------------------------*/

#ifndef SORTNODES_H
#define SORTNODES_H

#include <BALL/DATATYPE/GRAPH/VFLIB/state.h>

namespace BALL
{
	namespace VFLib
	{

		typedef int (*compare_fn)(const void *, const void *);

		struct NodeInfo
		{
			node_id id;
			node_id in;
			node_id out;
		};

		static int nodeInfoComp1(NodeInfo *a, NodeInfo *b);
		static int nodeInfoComp2(NodeInfo *a, NodeInfo *b);

		/*----------------------------------------------------
		 * Sorts the nodes of a graphs, returning a
		 * heap-allocated vector (using new) with the node ids
		 * in the proper orders.
		 * The sorting criterion takes into account:
		 *    1 - The number of nodes with the same in/out
		 *        degree.
		 *    2 - The valence of the nodes.
		 * The nodes at the beginning of the vector are
		 * the most singular, from which the matching should
		 * start.
		 *--------------------------------------------------*/
		template <class Graph>
		node_id* SortNodesByFrequency(const Graph& g)
		{
			int n = boost::num_vertices(g);
			NodeInfo *vect=new NodeInfo[n];
			int i = 0;
			typename boost::graph_traits<Graph>::vertex_iterator vi, vi_end;
			boost::tie(vi, vi_end) = boost::vertices(g);

			for (; vi != vi_end; ++vi, ++i)
			{
				vect[i].id = i;
				vect[i].in = boost::in_degree(*vi, g);
				vect[i].out= boost::out_degree(*vi, g);
			}

			qsort(vect, n, sizeof(vect[0]), (compare_fn)nodeInfoComp1);

			int run=1;

			for (i=0; i<n; i+=run)
			{
				for (run=1; i+run<n &&
				     vect[i+run].in==vect[i].in &&
				     vect[i+run].out==vect[i].out;
				     run++)
					;

				int j;

				for (j=0; j<run; j++)
				{
					vect[i+j].in += vect[i+j].out;
					vect[i+j].out=run;
				}
			}

			qsort(vect, n, sizeof(vect[0]), (compare_fn)nodeInfoComp2);

			node_id *nodes=new node_id[n];

			for (i=0; i<n; i++)
				nodes[i]=vect[i].id;

			delete[] vect;

			return nodes;
		}

		/**
		 * The ordering by in/out degree
		 */
		static int nodeInfoComp1(NodeInfo *a, NodeInfo *b)
		{
			if (a->out < b->out)
				return -1;
			else if (a->out > b->out)
				return +1;
			else if (a->in < b->in)
				return -1;
			else if (a->in > b->in)
				return +1;
			else
				return 0;
		}

		/**
		 * The ordering by frequency/valence.
		 * The frequency is in the out field, the valence in `in'.
		 */
		static int nodeInfoComp2(NodeInfo *a, NodeInfo *b)
		{
			if (a->in==0 && b->in !=0)
				return +1;
			else if (a->in!=0 && b->in==0)
				return -1;
			else if (a->out < b->out)
				return -1;
			else if (a->out > b->out)
				return +1;
			else if (a->in < b->in)
				return -1;
			else if (a->in > b->in)
				return +1;
			else
				return 0;
		}

	}
}

#endif
