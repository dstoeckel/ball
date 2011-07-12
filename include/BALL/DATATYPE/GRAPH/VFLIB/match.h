/*------------------------------------------------------------------
 * match.h
 * Header of match.cc
 * Declaration of the match function
 *
 * Author: P. Foggia
 * $Id: match.h,v 1.1 1998/09/29 09:49:48 foggia Exp $
 *-----------------------------------------------------------------*/


/*-----------------------------------------------------------------
 * REVISION HISTORY
 *   $Log: match.h,v $
 *   Revision 1.1  1998/09/29 09:49:48  foggia
 *   Initial revision
 *
 *----------------------------------------------------------------*/


#ifndef MATCH_H
#define MATCH_H

#include <BALL/DATATYPE/GRAPH/VFLIB/state.h>

namespace BALL
{
	namespace VFLib
	{

		typedef bool (*match_visitor)(int n, node_id c1[], node_id c2[],
		                              void *usr_data);

		/*-------------------------------------------------------------
		 * static bool match(pn, c1, c2, s)
		 * Finds a matching between two graphs, if it exists, starting
		 * from state s.
		 * Returns true a match has been found.
		 * *pn is assigned the numbero of matched nodes, and
		 * c1 and c2 will contain the ids of the corresponding nodes
		 * in the two graphs.
		 ------------------------------------------------------------*/
		template<class StateTemp>
		bool match(int *pn, node_id c1[], node_id c2[], StateTemp *s)
		{
			if (s->IsGoal())
			{
				*pn=s->CoreLen();
				s->GetCoreSet(c1, c2);
				return true;
			}

			if (s->IsDead())
				return false;

			node_id n1=NULL_NODE, n2=NULL_NODE;
			bool found=false;

			while (!found && s->NextPair(&n1, &n2, n1, n2))
			{
				if (s->IsFeasiblePair(n1, n2))
				{
					StateTemp *s1=s->Clone();
					s1->AddPair(n1, n2);
					found=match(pn, c1, c2, s1);
					s1->BackTrack();
					delete s1;
				}
			}

			return found;
		}


		template<class StateTemp, class MatchVisitor>
		bool match(node_id c1[], node_id c2[], MatchVisitor& vis,
		           StateTemp *s, int *pcount);

		template<class StateTemp>
		bool maximize(StateTemp *s0, int *pn, node_id c1[], node_id c2[]);

		/**
		 * Visits all the matches between two graphs, given the
		 * initial state of the match.
		 * Stops when there are no more matches, or the visitor vis
		 * returns true.
		 * @param vis A match visitor is a function object that is invoked for
		 *            each match that has been found.
		 *            If it returns false, then the next match is
		 *            searched; else the seach process terminates.
		 * @return the number of visited matches.
		 * @throw Exception::OutOfMemory if not enough memory for the matches could be allocated
		 */
		template<class StateTemp, class MatchVisitor>
		int match(StateTemp *s0, MatchVisitor& vis)
		{
			const typename StateTemp::QueryGraph*  g1 = s0->GetGraph1();
			const typename StateTemp::TargetGraph* g2 = s0->GetGraph2();

			/* Choose a conservative dimension for the arrays */
			int n;

			if (boost::num_vertices(*g1) < boost::num_vertices(*g2))
				n=boost::num_vertices(*g2);
			else
				n=boost::num_vertices(*g1);

			node_id *c1 = new node_id[n];
			node_id *c2 = new node_id[n];

			if (!c1 || !c2)
				throw Exception::OutOfMemory(__FILE__, __LINE__);

			int count=0;
			match(c1, c2, vis, s0, &count);

			delete[] c1;
			delete[] c2;
			return count;
		}

		/*-------------------------------------------------------------
		 * bool match(s0, pn, c1, c2)
		 * Finds a matching between two graph, if it exists, given the
		 * initial state of the matching process.
		 * Returns true a match has been found.
		 * *pn is assigned the number of matched nodes, and
		 * c1 and c2 will contain the ids of the corresponding nodes
		 * in the two graphs
		 ------------------------------------------------------------*/
		template<class StateTemp>
		bool match(StateTemp *s0, int *pn, node_id c1[], node_id c2[])
		{
			return match(pn,c1,c2,s0);
		}

		/*-------------------------------------------------------------
		 * static bool match(c1, c2, vis, usr_data, pcount)
		 * Visits all the matchings between two graphs,  starting
		 * from state s.
		 * Returns true if the caller must stop the visit.
		 * Stops when there are no more matches, or the visitor vis
		 * returns true.
		 ------------------------------------------------------------*/
		template<class StateTemp, class MatchVisitor>
		static bool match(node_id c1[], node_id c2[],
		                  MatchVisitor& vis, StateTemp *s, int *pcount)
		{
			if (s->IsGoal())
			{
				++*pcount;
				int n=s->CoreLen();
				s->GetCoreSet(c1, c2);
				return vis(n, c1, c2);
			}

			if (s->IsDead())
				return false;

			node_id n1=NULL_NODE, n2=NULL_NODE;

			while (s->NextPair(&n1, &n2, n1, n2))
			{
				if (s->IsFeasiblePair(n1, n2))
				{
					State<typename StateTemp::QueryGraph, typename StateTemp::TargetGraph>* s1 = s->Clone();
					s1->AddPair(n1, n2);

					if (match(c1, c2, vis, s1, pcount))
					{
						s1->BackTrack();
						delete s1;
						return true;
					}
					else
					{
						s1->BackTrack();
						delete s1;
					}
				}
			}

			return false;
		}


		/*-------------------------------------------------------------
		 * static bool maximize(s0, pn, c1, c2)
		 * Visits all the matchings between two graphs,  starting
		 * from state s.
		 * Returns true if the caller must stop the visit.
		 * Stops when there are no more matches, or the visitor vis
		 * returns true.
		 ------------------------------------------------------------*/
		template<class StateTemp>
		static bool maximize(StateTemp* s, node_id c1[], node_id c2[], int *pn)
		{
			if (s->CoreLen() > *pn) /* if we've found a core better than the current */
			{
				*pn=s->CoreLen();
				s->GetCoreSet(c1, c2);
				// return false; /* keep going, we might find a bigger one yet */
			}

			if (s->IsGoal())
			{
				// we found a complete overlap, it doesn't get much bigger than that.
				return true;
			}

			if (s->IsDead())
				return false;

			node_id n1=NULL_NODE, n2=NULL_NODE;
			bool found=false;

			while (!found && s->NextPair(&n1, &n2, n1, n2))
			{
				if (s->IsFeasiblePair(n1, n2))
				{
					State<typename StateTemp::QueryGraph, typename StateTemp::TargetGraph>* s1 = s->Clone();
					s1->AddPair(n1, n2);
					found = maximize(s1, c1, c2, pn);
					s1->BackTrack();
					delete s1;
				}
			}

			return found;
		}

		template<class StateTemp>
		bool maximize(StateTemp* s0, int *pn, node_id c1[], node_id c2[])
		{
			// here, pn serves as a save point for the size of the current maximal core.
			// likewise, c1 and c2 store the current maximum state.
			return maximize(s0, c1, c2, pn);
		}

	}
}
#endif //MATCH_H
