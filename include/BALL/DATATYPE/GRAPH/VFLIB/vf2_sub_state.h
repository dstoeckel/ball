/*------------------------------------------------------------
 * vf2_sub_state.h
 * Interface of vf2_sub_state.cc
 * Definition of a class representing a state of the matching
 * process between two ARGs.
 * See: argraph.h state.h
 *
 * Author: P. Foggia
 *-----------------------------------------------------------------*/

#ifndef VF2_SUB_STATE_H
#define VF2_SUB_STATE_H

#include <BALL/DATATYPE/GRAPH/VFLIB/state.h>
#include <BALL/DATATYPE/GRAPH/VFLIB/sortnodes.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <cassert>

namespace BALL
{
	namespace VFLib
	{

		/*----------------------------------------------------------
		 * class VF2SubState
		 * A representation of the SSR current state
		 * See vf2_sub_state.cc for more details.
		 ---------------------------------------------------------*/
		template<class QueryGraph, class TargetGraph, class VertexComparator, class EdgeComparator>
		class VF2SubState: public State<QueryGraph, TargetGraph>
		{
			private:
				typedef State<QueryGraph, TargetGraph> BaseState;
				typedef typename boost::property_map<QueryGraph, boost::vertex_index_t>::type QueryVertexIndexMap;
				typedef typename boost::graph_traits<QueryGraph>::vertex_descriptor           QueryVertex;
				typedef typename boost::graph_traits<QueryGraph>::edge_descriptor             QueryEdge;
				typedef typename boost::graph_traits<QueryGraph>::vertex_iterator             QueryVertexIterator;
				typedef typename boost::graph_traits<QueryGraph>::edge_iterator               QueryEdgeIterator;
				typedef typename boost::graph_traits<QueryGraph>::out_edge_iterator           QueryOutEdgeIterator;
				typedef typename boost::graph_traits<QueryGraph>::in_edge_iterator            QueryInEdgeIterator;

				typedef typename boost::property_map<TargetGraph, boost::vertex_index_t>::type TargetVertexIndexMap;
				typedef typename boost::graph_traits<TargetGraph>::vertex_descriptor           TargetVertex;
				typedef typename boost::graph_traits<TargetGraph>::edge_descriptor             TargetEdge;
				typedef typename boost::graph_traits<TargetGraph>::edge_iterator               TargetEdgeIterator;
				typedef typename boost::graph_traits<TargetGraph>::out_edge_iterator           TargetOutEdgeIterator;
				typedef typename boost::graph_traits<TargetGraph>::in_edge_iterator            TargetInEdgeIterator;


				int core_len, orig_core_len;
				int added_node1;
				int t1both_len, t2both_len, t1in_len, t1out_len,
				t2in_len, t2out_len; // Core nodes are also counted by these...
				node_id *core_1;
				node_id *core_2;
				node_id *in_1;
				node_id *in_2;
				node_id *out_1;
				node_id *out_2;

				node_id *order;

				const QueryGraph& g1;
				const TargetGraph& g2;
				int n1, n2;

				const VertexComparator* v_comp;
				const EdgeComparator* e_comp;

				long* share_count;

			public:
				VF2SubState(const QueryGraph& g1, const TargetGraph& g2, bool sortNodes=false);
				VF2SubState(const VF2SubState &state);
				~VF2SubState();

				void setVertexComparator(const VertexComparator& vcomp)
				{
					if(*share_count == 1)
					{
						delete v_comp;
						v_comp = new VertexComparator(vcomp);
					}
				}

				void setEdgeComparator(const EdgeComparator& ecomp)
				{
					if(*share_count == 1)
					{
						delete e_comp;
						e_comp = new EdgeComparator(ecomp);
					}
				}

				const QueryGraph *GetGraph1() const
				{
					return &g1;
				}

				const TargetGraph *GetGraph2() const
				{
					return &g2;
				}

				bool NextPair(node_id *pn1, node_id *pn2,
				              node_id prev_n1=NULL_NODE, node_id prev_n2=NULL_NODE);
				bool IsFeasiblePair(node_id n1, node_id n2);
				void AddPair(node_id n1, node_id n2);

				bool IsGoal()
				{
					return core_len==n1 ;
				}

				bool IsDead()
				{
					return n1>n2  ||
					       t1both_len>t2both_len ||
					       t1out_len>t2out_len ||
					       t1in_len>t2in_len;
				}

				int CoreLen()
				{
					return core_len;
				}

				void GetCoreSet(node_id c1[], node_id c2[]);
				BaseState *Clone();

				virtual void BackTrack();
		};

		template<class QueryGraph, class TargetGraph, class VertexComparator, class EdgeComparator>
		VF2SubState<QueryGraph, TargetGraph, VertexComparator, EdgeComparator>::VF2SubState(const QueryGraph& ag1, const TargetGraph& ag2, bool sortNodes)
				: g1(ag1),
				  g2(ag2),
				  v_comp(0),
				  e_comp(0)
		{
			n1=num_vertices(g1);
			n2=num_vertices(g2);

			if (sortNodes)
				order = SortNodesByFrequency(g1);
			else
				order = NULL;

			core_len=orig_core_len=0;
			t1both_len=t1in_len=t1out_len=0;
			t2both_len=t2in_len=t2out_len=0;

			added_node1=NULL_NODE;

			core_1 = new node_id[n1];
			core_2 = new node_id[n2];
			in_1   = new node_id[n1];
			in_2   = new node_id[n2];
			out_1  = new node_id[n1];
			out_2  = new node_id[n2];
			share_count = new long;

			if (!core_1 || !core_2 || !in_1 || !in_2
			    || !out_1 || !out_2 || !share_count)
				throw Exception::OutOfMemory(__FILE__, __LINE__);

			int i;

			for (i=0; i<n1; i++)
			{
				core_1[i]=NULL_NODE;
				in_1[i]=0;
				out_1[i]=0;
			}

			for (i=0; i<n2; i++)
			{
				core_2[i]=NULL_NODE;
				in_2[i]=0;
				out_2[i]=0;
			}

			*share_count = 1;
		}


		/*----------------------------------------------------------
		 * VF2SubState::VF2SubState(state)
		 * Copy constructor.
		 ---------------------------------------------------------*/
		template<class QueryGraph, class TargetGraph, class VertexComparator, class EdgeComparator>
		VF2SubState<QueryGraph,TargetGraph,VertexComparator,EdgeComparator>::VF2SubState(const VF2SubState &state)
				: g1(state.g1),
				  g2(state.g2),
				  n1(state.n1),
				  n2(state.n2),
				  order(state.order),
				  core_len(state.core_len),
				  orig_core_len(state.core_len),
				  t1in_len(state.t1in_len),
				  t1out_len(state.t1out_len),
				  t1both_len(state.t1both_len),
				  t2in_len(state.t2in_len),
				  t2out_len(state.t2out_len),
				  t2both_len(state.t2both_len),
				  added_node1(NULL_NODE),
				  core_1(state.core_1),
				  core_2(state.core_2),
				  in_1(state.in_1),
				  in_2(state.in_2),
				  out_1(state.out_1),
				  out_2(state.out_2),
				  share_count(state.share_count),
				  v_comp(state.v_comp),
				  e_comp(state.e_comp)
		{
			++*share_count;
		}

		/*---------------------------------------------------------------
		 * VF2SubState::~VF2SubState()
		 * Destructor.
		 --------------------------------------------------------------*/
		template<class QueryGraph, class TargetGraph, class VertexComparator, class EdgeComparator>
		VF2SubState<QueryGraph, TargetGraph, VertexComparator, EdgeComparator>::~VF2SubState()
		{
			if (-- *share_count == 0)
			{
				delete [] core_1;
				delete [] core_2;
				delete [] in_1;
				delete [] out_1;
				delete [] in_2;
				delete [] out_2;
				delete share_count;
				delete [] order;
				delete v_comp;
				delete e_comp;
			}
		}

		/*--------------------------------------------------------------------------
		 * bool VF2SubState::NextPair(pn1, pn2, prev_n1, prev_n2)
		 * Puts in *pn1, *pn2 the next pair of nodes to be tried.
		 * prev_n1 and prev_n2 must be the last nodes, or NULL_NODE (default)
		 * to start from the first pair.
		 * Returns false if no more pairs are available.
		 -------------------------------------------------------------------------*/
		template<class QueryGraph, class TargetGraph, class VertexComparator, class EdgeComparator>
		bool VF2SubState<QueryGraph, TargetGraph, VertexComparator, EdgeComparator>::NextPair(node_id *pn1, node_id *pn2,
		                                                                                      node_id prev_n1, node_id prev_n2)
		{
			if (prev_n1==NULL_NODE)
				prev_n1=0;

			if (prev_n2==NULL_NODE)
				prev_n2=0;
			else
				prev_n2++;

			if (t1both_len>core_len && t2both_len>core_len)
			{
				while (prev_n1<n1 &&
				       (core_1[prev_n1]!=NULL_NODE || out_1[prev_n1]==0
				        || in_1[prev_n1]==0) )
				{
					prev_n1++;
					prev_n2=0;
				}
			}
			else if (t1out_len>core_len && t2out_len>core_len)
			{
				while (prev_n1<n1 &&
				       (core_1[prev_n1]!=NULL_NODE || out_1[prev_n1]==0) )
				{
					prev_n1++;
					prev_n2=0;
				}
			}
			else if (t1in_len>core_len && t2in_len>core_len)
			{
				while (prev_n1<n1 &&
				       (core_1[prev_n1]!=NULL_NODE || in_1[prev_n1]==0) )
				{
					prev_n1++;
					prev_n2=0;
				}
			}
			else if (prev_n1==0 && order!=NULL)
			{
				int i=0;

				while (i<n1 && core_1[prev_n1=order[i]] != NULL_NODE)
					i++;

				if (i==n1)
					prev_n1=n1;
			}
			else
			{
				while (prev_n1<n1 && core_1[prev_n1]!=NULL_NODE )
				{
					prev_n1++;
					prev_n2=0;
				}
			}



			if (t1both_len>core_len && t2both_len>core_len)
			{
				while (prev_n2<n2 &&
				       (core_2[prev_n2]!=NULL_NODE || out_2[prev_n2]==0
				        || in_2[prev_n2]==0) )
				{
					prev_n2++;
				}
			}
			else if (t1out_len>core_len && t2out_len>core_len)
			{
				while (prev_n2<n2 &&
				       (core_2[prev_n2]!=NULL_NODE || out_2[prev_n2]==0) )
				{
					prev_n2++;
				}
			}
			else if (t1in_len>core_len && t2in_len>core_len)
			{
				while (prev_n2<n2 &&
				       (core_2[prev_n2]!=NULL_NODE || in_2[prev_n2]==0) )
				{
					prev_n2++;
				}
			}
			else
			{
				while (prev_n2<n2 && core_2[prev_n2]!=NULL_NODE )
				{
					prev_n2++;
				}
			}


			if (prev_n1<n1 && prev_n2<n2)
			{
				*pn1=prev_n1;
				*pn2=prev_n2;
				return true;
			}

			return false;
		}



		/*---------------------------------------------------------------
		 * bool VF2SubState::IsFeasiblePair(node1, node2)
		 * Returns true if (node1, node2) can be added to the state
		 * NOTE:
		 *   The attribute compatibility check (methods CompatibleNode
		 *   and CompatibleEdge of ARGraph) is always performed
		 *   applying the method to g1, and passing the attribute of
		 *   g1 as first argument, and the attribute of g2 as second
		 *   argument. This may be important if the compatibility
		 *   criterion is not symmetric.
		 --------------------------------------------------------------*/
		template<class QueryGraph, class TargetGraph, class VertexComparator, class EdgeComparator>
		bool VF2SubState<QueryGraph, TargetGraph, VertexComparator, EdgeComparator>::IsFeasiblePair(node_id node1, node_id node2)
		{
			assert(node1<n1);
			assert(node2<n2);
			assert(core_1[node1]==NULL_NODE);
			assert(core_2[node2]==NULL_NODE);

			if (v_comp && !v_comp->compare(boost::vertex(node1, g1), boost::vertex(node2, g2)))
				return false;

			int i, other1, other2;
			void *attr1;
			int termout1=0, termout2=0, termin1=0, termin2=0, new1=0, new2=0;

			// Check the 'out' edges of node1
			QueryOutEdgeIterator oe, oe_end;
			boost::tie(oe, oe_end) = boost::out_edges(boost::vertex(node1, g1), g1);
			QueryVertexIndexMap g1_indices = boost::get(boost::vertex_index, g1);
			TargetVertex node2_vert = boost::vertex(node2, g2);

			for (; oe != oe_end; ++oe)
			{
				other1=boost::get(g1_indices, boost::target(*oe, g1));

				if (core_1[other1] != NULL_NODE)
				{
					other2=core_1[other1];
					TargetEdge tmp_edge;
					bool found;
					boost::tie(tmp_edge, found) = boost::edge(node2_vert, boost::vertex(other2, g2), g2);

					if (!found || (e_comp && !e_comp->compare(*oe, tmp_edge)))
						return false;
				}
				else
				{
					if (in_1[other1])
						termin1++;

					if (out_1[other1])
						termout1++;

					if (!in_1[other1] && !out_1[other1])
						new1++;
				}
			}

			QueryInEdgeIterator ie, ie_end;
			boost::tie(ie, ie_end) = boost::in_edges(boost::vertex(node1, g1), g1);

			// Check the 'in' edges of node1
			for (; ie != ie_end; ++ie)
			{
				other1=boost::get(g1_indices, boost::source(*ie,g1));

				if (core_1[other1]!=NULL_NODE)
				{
					other2=core_1[other1];
					TargetEdge tmp_edge;
					bool found;
					boost::tie(tmp_edge, found) = boost::edge(boost::vertex(other2, g2), node2_vert, g2); // note: other way round than in the 'out' edge test!

					if (!found || (e_comp && !e_comp->compare(*ie, tmp_edge)))
						return false;
				}
				else
				{
					if (in_1[other1])
						termin1++;

					if (out_1[other1])
						termout1++;

					if (!in_1[other1] && !out_1[other1])
						new1++;
				}
			}

			TargetOutEdgeIterator te, te_end;
			boost::tie(te, te_end) = boost::out_edges(node2_vert,g2);
			TargetVertexIndexMap g2_indices = boost::get(boost::vertex_index, g2);
			QueryVertex node1_vert = boost::vertex(node1, g1);

			// Check the 'out' edges of node2
			for (; te != te_end; ++te)
			{
				other2=boost::get(g2_indices, boost::target(*te, g2));

				if (core_2[other2]!=NULL_NODE)
				{
					other1=core_2[other2];

					if (!boost::edge(node1_vert, boost::vertex(other1, g1), g1).second)
						return false;
				}
				else
				{
					if (in_2[other2])
						termin2++;

					if (out_2[other2])
						termout2++;

					if (!in_2[other2] && !out_2[other2])
						new2++;
				}
			}

			TargetInEdgeIterator ite, ite_end;
			boost::tie(ite, ite_end) = boost::in_edges(boost::vertex(node2, g2),g2);

			// Check the 'in' edges of node2
			for (; ite != ite_end; ++ite)
			{
				other2=boost::get(g2_indices, boost::source(*ite,g2));

				if (core_2[other2] != NULL_NODE)
				{
					other1=core_2[other2];

					if (!boost::edge(node1_vert, boost::vertex(other1, g1), g1).second)
						return false;
				}
				else
				{
					if (in_2[other2])
						termin2++;

					if (out_2[other2])
						termout2++;

					if (!in_2[other2] && !out_2[other2])
						new2++;
				}
			}

			return termin1<=termin2 && termout1<=termout2 && new1<=new2;
		}



		/*--------------------------------------------------------------
		 * void VF2SubState::AddPair(node1, node2)
		 * Adds a pair to the Core set of the state.
		 * Precondition: the pair must be feasible
		 -------------------------------------------------------------*/
		template<class QueryGraph, class TargetGraph, class VertexComparator, class EdgeComparator>
		void VF2SubState<QueryGraph, TargetGraph, VertexComparator, EdgeComparator>::AddPair(node_id node1, node_id node2)
		{
			assert(node1<n1);
			assert(node2<n2);
			assert(core_len<n1);
			assert(core_len<n2);

			core_len++;
			added_node1=node1;

			if (!in_1[node1])
			{
				in_1[node1]=core_len;
				t1in_len++;

				if (out_1[node1])
					t1both_len++;
			}

			if (!out_1[node1])
			{
				out_1[node1]=core_len;
				t1out_len++;

				if (in_1[node1])
					t1both_len++;
			}

			if (!in_2[node2])
			{
				in_2[node2]=core_len;
				t2in_len++;

				if (out_2[node2])
					t2both_len++;
			}

			if (!out_2[node2])
			{
				out_2[node2]=core_len;
				t2out_len++;

				if (in_2[node2])
					t2both_len++;
			}

			core_1[node1]=node2;
			core_2[node2]=node1;


			int i, other;
			QueryInEdgeIterator qe, qe_end;
			boost::tie(qe, qe_end) = boost::in_edges(boost::vertex(node1, g1),g1);
			QueryVertexIndexMap g1_indices = boost::get(boost::vertex_index,g1);

			for (; qe != qe_end; ++qe)
			{
				other=boost::get(g1_indices,boost::source(*qe,g1));

				if (!in_1[other])
				{
					in_1[other]=core_len;
					t1in_len++;

					if (out_1[other])
						t1both_len++;
				}
			}

			QueryOutEdgeIterator ie, ie_end;
			boost::tie(ie, ie_end) = boost::out_edges(boost::vertex(node1, g1),g1);

			for (; ie != ie_end; ++ie)
			{
				other=boost::get(g1_indices,boost::target(*ie,g1));

				if (!out_1[other])
				{
					out_1[other]=core_len;
					t1out_len++;

					if (in_1[other])
						t1both_len++;
				}
			}

			TargetInEdgeIterator ite, ite_end;
			boost::tie(ite, ite_end) = boost::in_edges(boost::vertex(node2,g2),g2);
			TargetVertexIndexMap g2_indices = boost::get(boost::vertex_index, g2);

			for (; ite != ite_end; ++ite)
			{
				other=boost::get(g2_indices,boost::source(*ite,g2));

				if (!in_2[other])
				{
					in_2[other]=core_len;
					t2in_len++;

					if (out_2[other])
						t2both_len++;
				}
			}

			TargetOutEdgeIterator te, te_end;
			boost::tie(te, te_end) = boost::out_edges(boost::vertex(node2,g2),g2);

			for (; te != te_end; ++te)
			{
				other=boost::get(g2_indices,boost::target(*te,g2));

				if (!out_2[other])
				{
					out_2[other]=core_len;
					t2out_len++;

					if (in_2[other])
						t2both_len++;
				}
			}

		}



		/*--------------------------------------------------------------
		 * void VF2SubState::GetCoreSet(c1, c2)
		 * Reads the core set of the state into the arrays c1 and c2.
		 * The i-th pair of the mapping is (c1[i], c2[i])
		 --------------------------------------------------------------*/
		template<class QueryGraph, class TargetGraph, class VertexComparator, class EdgeComparator>
		void VF2SubState<QueryGraph, TargetGraph, VertexComparator, EdgeComparator>::GetCoreSet(node_id c1[], node_id c2[])
		{
			int i,j;

			for (i=0,j=0; i<n1; i++)
				if (core_1[i] != NULL_NODE)
				{
					c1[j]=i;
					c2[j]=core_1[i];
					j++;
				}
		}


		/*----------------------------------------------------------------
		 * Clones a VF2SubState, allocating with new the clone.
		 --------------------------------------------------------------*/
		template<class QueryGraph, class TargetGraph, class VertexComparator, class EdgeComparator>
		State<QueryGraph, TargetGraph>* VF2SubState<QueryGraph, TargetGraph, VertexComparator, EdgeComparator>::Clone()
		{
			return new VF2SubState(*this);
		}

		/*----------------------------------------------------------------
		 * Undoes the changes to the shared vectors made by the
		 * current state. Assumes that at most one AddPair has been
		 * performed.
		 ----------------------------------------------------------------*/
		template<class QueryGraph, class TargetGraph, class VertexComparator, class EdgeComparator>
		void VF2SubState<QueryGraph, TargetGraph, VertexComparator, EdgeComparator>::BackTrack()
		{
			assert(core_len - orig_core_len <= 1);
			assert(added_node1 != NULL_NODE);

			if (orig_core_len < core_len)
			{
				int node2;

				QueryInEdgeIterator qe, qe_end;
				boost::tie(qe, qe_end) = boost::in_edges(boost::vertex(added_node1, g1),g1);
				QueryVertexIndexMap g1_indices = boost::get(boost::vertex_index,g1);

				if (in_1[added_node1] == core_len)
					in_1[added_node1] = 0;

				for (; qe != qe_end; ++qe)
				{
					int other=boost::get(g1_indices, boost::source(*qe, g1));

					if (in_1[other]==core_len)
						in_1[other]=0;
				}

				QueryOutEdgeIterator oqe, oqe_end;
				boost::tie(oqe, oqe_end) = boost::out_edges(boost::vertex(added_node1, g1),g1);

				if (out_1[added_node1] == core_len)
					out_1[added_node1] = 0;

				for (; oqe != oqe_end; ++oqe)
				{
					int other=boost::get(g1_indices, boost::target(*oqe, g1));

					if (out_1[other]==core_len)
						out_1[other]=0;
				}

				node2 = core_1[added_node1];

				TargetInEdgeIterator te, te_end;
				boost::tie(te, te_end) = boost::in_edges(boost::vertex(node2, g2),g2);
				TargetVertexIndexMap g2_indices = boost::get(boost::vertex_index,g2);

				if (in_2[node2] == core_len)
					in_2[node2] = 0;

				for (; te != te_end; ++te)
				{
					int other=boost::get(g2_indices, boost::source(*te, g2));

					if (in_2[other]==core_len)
						in_2[other]=0;
				}

				TargetOutEdgeIterator ote, ote_end;
				boost::tie(ote, ote_end) = boost::out_edges(boost::vertex(node2, g2),g2);

				if (out_2[node2] == core_len)
					out_2[node2] = 0;

				for (; ote != ote_end; ++ote)
				{
					int other=boost::get(g2_indices, boost::target(*ote, g2));

					if (out_2[other]==core_len)
						out_2[other]=0;
				}

				core_1[added_node1] = NULL_NODE;
				core_2[node2] = NULL_NODE;

				core_len=orig_core_len;
				added_node1 = NULL_NODE;
			}

		}


	} // namespace VFLib
} // namespace BALL

#endif //VF2_SUB_STATE_H
