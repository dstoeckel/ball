/*------------------------------------------------------------
 * state.h
 * Definition of an abstract class representing a state of the
 * matching process between two ARGs.
 * See: argraph.h
 *
 * Author: P. Foggia
 * $Id: state.h,v 1.3 1998/09/29 09:50:16 foggia Exp $
 *-----------------------------------------------------------------*/


/*-----------------------------------------------------------------
 * REVISION HISTORY
 *   $Log: state.h,v $
 *   Revision 1.3  1998/09/29 09:50:16  foggia
 *   Ora definisce una classe astratta State, da cui discendono
 *   VFState e UllState
 *
 *   Revision 1.2  1998/09/26 09:02:32  foggia
 *   minor changes
 *
 *   Revision 1.1  1998/09/19 14:40:35  foggia
 *   Initial revision
 *
 *----------------------------------------------------------------*/


#ifndef BALL_DATATYPE_GRAPH_VFLIB_STATE_H
#define BALL_DATATYPE_GRAPH_VFLIB_STATE_H

namespace BALL
{
	namespace VFLib
	{

		typedef unsigned short node_id;
		const node_id NULL_NODE=0xFFFF;

		/*----------------------------------------------------------
		 * class State
		 * An abstract representation of the SSR current state.
		 * NOTE: Respect to pre-2.0 version of the library, class
		 *   State assumes explicitly a depth-first search. The
		 *   BackTrack method has been added to allow a state
		 *   to clean up things before reverting to its parent. This
		 *   can be used, for instance, for sharing resources between
		 *   the parent and the child. The BackTrack implementation
		 *   can safely assume that at most one AddPair has been
		 *   performed on the state.
		 ---------------------------------------------------------*/
		template<class QueryGraph_t, class TargetGraph_t>
		class State
		{
			public:
				typedef QueryGraph_t QueryGraph;
				typedef TargetGraph_t TargetGraph;

				virtual ~State() {}
				virtual const QueryGraph *GetGraph1() const =0;
				virtual const TargetGraph *GetGraph2() const =0;
				virtual bool NextPair(node_id *pn1, node_id *pn2,
				                      node_id prev_n1=NULL_NODE, node_id prev_n2=NULL_NODE)=0;
				virtual bool IsFeasiblePair(node_id n1, node_id n2)=0;
				virtual void AddPair(node_id n1, node_id n2)=0;
				virtual bool IsGoal() =0;
				virtual bool IsDead() =0;
				virtual int CoreLen() const =0;
				virtual void GetCoreSet(node_id c1[], node_id c2[]) const =0;
				virtual State *Clone() =0;  // Changed clone to Clone for uniformity

				virtual void BackTrack() { }
		};

	}
}
#endif //BALL_DATATYPE_GRAPH_VFLIB_STATE_H
