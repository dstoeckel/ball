/*------------------------------------------------------------------
 * vf2_mcs_state.cc
 * Implementation of the class VF2MCSState
 *
 * A State Space Search State for the Maximum Common Subgraph (MCS)
 * problem. Adapted from VF2SubState; using information from
 *   Conte, Donatello and Foggia, Pasquale and Vento, Mario;
 *   Challenging Complexity of Maximum Common Subgraph Detection 
 *   Algorithms [...];
 *   Journal of Graph Algorithms and Applications, Vol 11, No. 1,
 *   pp 99-143
 *
 * Author: W. Herget <wolfgang@r007.de>
 *-----------------------------------------------------------------*/



/*-----------------------------------------------------------------
 * NOTE:
 *   The attribute compatibility check (methods CompatibleNode
 *   and CompatibleEdge of ARGraph) is always performed
 *   applying the method to g1, and passing the attribute of
 *   g1 as first argument, and the attribute of g2 as second
 *   argument. This may be important if the compatibility
 *   criterion is not symmetric.
 -----------------------------------------------------------------*/


/*---------------------------------------------------------
 *   IMPLEMENTATION NOTES:
 * The vectors core_1, core_2, and exhausted_1
 * are shared among the instances of this class; they are
 * owned by the instance with core_len==0 (the root of the
 * SSR).
 * In the vector exhausted_1 there is a value indicating
 * the level at which the corresponding node exhausted all
 * possible pairing choices in the second graph.
 * largest_known_mcs holds the size (number of paired nodes)
 * of the current (global) maximum common subgraph.
 * This information is used for backtracking.
 ---------------------------------------------------------*/


#include <stddef.h>

#include "BALL/STRUCTURE/VFLIB/vf2_mcs_state.h"
#include "BALL/STRUCTURE/VFLIB/sortnodes.h"

#include "BALL/STRUCTURE/VFLIB/error.h"

namespace BALL {
namespace VFLib {

/*----------------------------------------------------------
 * Methods of the class VF2SubState
 ---------------------------------------------------------*/

/*----------------------------------------------------------
 * VF2MCSState::VF2MCSState(g1, g2, sortNodes)
 * Constructor. Makes an empty state.
 * If sortNodes is true, computes an initial ordering
 * for the nodes based on the frequency of their valence.
 ---------------------------------------------------------*/
VF2MCSState::VF2MCSState(Graph *ag1, Graph *ag2, bool sortNodes)
  { g1=ag1;
    g2=ag2;
    n1=g1->NodeCount();
    n2=g2->NodeCount();

    if (sortNodes)
      order = SortNodesByFrequency(ag1);
    else
      order = NULL;

    core_len=orig_core_len=0;
    exhaust_len1=0;

    added_node1=NULL_NODE;

    core_1=new node_id[n1];
    core_2=new node_id[n2];
    exhausted_1=new node_id[n1];
    share_count = new long;
    longest_known_mcs = new int;
    if (!core_1 || !core_2 || !exhausted_1  || !share_count || !longest_known_mcs)
      error("Out of memory");

    int i;
    for(i=0; i<n1; i++)
      {
        core_1[i]=NULL_NODE;
        exhausted_1[i]=0;
      }
    for(i=0; i<n2; i++)
      {
        core_2[i]=NULL_NODE;
      }

    *share_count = 1;
    *longest_known_mcs = 0;
  }


/*----------------------------------------------------------
 * VF2MCSState::VF2MCSState(state)
 * Copy constructor.
 ---------------------------------------------------------*/
VF2MCSState::VF2MCSState(const VF2MCSState &state)
  { g1=state.g1;
    g2=state.g2;
    n1=state.n1; // #node in g1
    n2=state.n2; // #node in g2

    order=state.order; // optional sorting of nodes: g1 -> g1

    // length of mapping (# of nodes in the current state)
    core_len=orig_core_len=state.core_len;

    /* terminal sets:
      (predecessor of node: is source of path to node)
      (successor of node: is destination of path from node)

       t_x^in: set of nodes not mapped, but predecessors of a node that's mapped
       t_x^out: set of nodes not mapped, but successor of a node that's mapped
    */
    exhaust_len1=state.exhaust_len1; // # incoming edges into the mapped portion of g1

    added_node1=NULL_NODE; // the node added in this state, from g1
      // its partner is stored in core_1[added_node1]

    // carry over the pointers (for memory saving)
    core_1=state.core_1; // g1[node_id] -> g2[node_id]
    core_2=state.core_2;
    exhausted_1=state.exhausted_1; // earliest search space tree depth at which an edge FROM
    share_count=state.share_count;
    longest_known_mcs=state.longest_known_mcs;

    ++ *share_count;

  }


/*---------------------------------------------------------------
 * VF2SubState::~VF2SubState()
 * Destructor.
 --------------------------------------------------------------*/
VF2MCSState::~VF2MCSState()
  { if (-- *share_count == 0) // destruction of the initial state cleans up.
    { delete [] core_1;
      delete [] core_2;
      delete [] exhausted_1;
      delete share_count;
      delete [] order;
    }
  }

/*--------------------------------------------------------------------------
 * bool VF2SubState::IsDead()
 * determines if it is possibly fruitful to continue exploring
 * the search space from this state onward.
 -------------------------------------------------------------------------*/
bool VF2MCSState::IsDead()
{
  // cannot possibly make mcs larger
  if (n1 - exhaust_len1 <= *longest_known_mcs) return true;
  // or used up all nodes
  if (core_len+exhaust_len1 >= n1) return true;

  return false;
}

/*--------------------------------------------------------------------------
 * bool VF2SubState::NextPair(pn1, pn2, prev_n1, prev_n2)
 * Puts in *pn1, *pn2 the next pair of nodes to be tried.
 * prev_n1 and prev_n2 must be the last nodes, or NULL_NODE (default)
 * to start from the first pair.
 * Returns false if no more pairs are available.
 -------------------------------------------------------------------------*/
bool VF2MCSState::NextPair(node_id *pn1, node_id *pn2,
              node_id prev_n1, node_id prev_n2)
  {

  node_id next_n1, next_n2;

  // if next_nX given (!= NULL_NODE), continue from there. else restart
  if (prev_n1 == NULL_NODE) {
    next_n1 = 0;
  } else {
    next_n1 = prev_n1;
  }

  if (prev_n2 == NULL_NODE) {
    next_n2 = 0;
  } else {
    next_n2 = prev_n2 + 1;
  }

  // this next_n1 has exhausted its possibilities at this level.
  if (next_n2 >= n2) {
    exhausted_1[next_n1] = core_len;
    exhaust_len1++;
  }

  // n1, the number of nodes in g1 is always bigger then the last node_id of g1
  while ((next_n1 < n1) && (exhausted_1[next_n1] || (core_1[next_n1] != NULL_NODE))) {
    next_n1++; next_n2 = 0;
  }
  while ((next_n2 < n2) && (core_2[next_n2] != NULL_NODE)) next_n2++;

  if ((next_n1 < n1) && (next_n2 < n2)) {
    *pn1 = next_n1;
    *pn2 = next_n2;
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
bool VF2MCSState::IsFeasiblePair(node_id node1, node_id node2)
  {

  if (!g1->CompatibleNode(g1->GetNodeAttr(node1),g2->GetNodeAttr(node2)))
    return false;

  // assure mutual connections:
  // if node1 to be added is connected to something
  // already in core1 and its image is also in core2
  // make sure a connection also exists in g2.

  void* attr;
  node_id partner1, partner2;

  // edges incoming to node1 in g1
  for (int e1 = 0; e1 < g1->InEdgeCount(node1); e1++)
  {
    partner1 = g1->GetInEdge(node1,e1,&attr);
    if ((partner2 = core_1[partner1]) != NULL_NODE) {
      if (!g2->HasEdge(node2,partner2) ||
          !g2->CompatibleEdge(attr,g2->GetEdgeAttr(node2,partner2)))
      {
        return false;
      }
    }
  }

  // edges outgoing from node1 in g2
  for (int e1 = 0; e1 < g1->OutEdgeCount(node1); e1++)
  {
    partner1 = g1->GetOutEdge(node1,e1,&attr);
    if ((partner2 = core_1[partner1]) != NULL_NODE) {
      if (!g2->HasEdge(node2,partner2) ||
          !g2->CompatibleEdge(attr,g2->GetEdgeAttr(node2,partner2)))
      {
        return false;
      }
    }
  }

  // edges incoming to node2 in g2
  for (int e2 = 0; e2 < g2->InEdgeCount(node2); e2++)
  {
    partner2 = g2->GetInEdge(node2,e2,&attr);
    if ((partner1 = core_2[partner2]) != NULL_NODE) {
      if (!g1->HasEdge(node1,partner1) ||
          !g1->CompatibleEdge(attr,g1->GetEdgeAttr(node1,partner1)))
      {
        return false;
      }
    }
  }

  // edges outgoing from node2 in g2
  for (int e2 = 0; e2 < g2->OutEdgeCount(node2); e2++)
  {
    partner2 = g2->GetOutEdge(node2,e2,&attr);
    if ((partner1 = core_2[partner2]) != NULL_NODE) {
      if (!g1->HasEdge(node1,partner1) ||
          !g1->CompatibleEdge(attr,g1->GetEdgeAttr(node1,partner1)))
      {
        return false;
      }
    }
  }

  // only completely disjoint nodes are allowed in sans
  // connectivity check.
  return true;

  }



/*--------------------------------------------------------------
 * void VF2SubState::AddPair(node1, node2)
 * Adds a pair to the Core set of the state.
 * Precondition: the pair must be feasible
 -------------------------------------------------------------*/
void VF2MCSState::AddPair(node_id node1, node_id node2)
  {

  core_1[node1] = node2;
  core_2[node2] = node1;

  added_node1 = node1;

  core_len++;

  if (core_len > *longest_known_mcs)
    *longest_known_mcs = core_len;
  }



/*--------------------------------------------------------------
 * void VF2SubState::GetCoreSet(c1, c2)
 * Reads the core set of the state into the arrays c1 and c2.
 * The i-th pair of the mapping is (c1[i], c2[i])
 --------------------------------------------------------------*/
void VF2MCSState::GetCoreSet(node_id c1[], node_id c2[])
  { int i,j;
    for (i=0,j=0; i<n1; i++)
      if (core_1[i] != NULL_NODE)
        { c1[j]=i;
          c2[j]=core_1[i];
          j++;
        }
  }


/*----------------------------------------------------------------
 * Clones a VF2SubState, allocating with new the clone.
 --------------------------------------------------------------*/
State* VF2MCSState::Clone()
  { return new VF2MCSState(*this);
  }

/*----------------------------------------------------------------
 * Undoes the changes to the shared vectors made by the
 * current state. Assumes that at most one AddPair has been
 * performed.
 ----------------------------------------------------------------*/
void VF2MCSState::BackTrack()
  {
  node_id added_node2 = core_1[added_node1];
  core_1[added_node1] = NULL_NODE;
  core_2[added_node2] = NULL_NODE;

  // clear exhausted flags for nodes added deeper in the SST
  for (int n = 0; n<n1; n++) if (exhausted_1[n] >= core_len) {
    exhausted_1[n] = 0;
    exhaust_len1--;
  }

  core_len--;
  }

}
}
