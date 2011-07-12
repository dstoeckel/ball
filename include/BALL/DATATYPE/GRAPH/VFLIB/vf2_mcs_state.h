#ifndef VF2_MCS_STATE_H
#define VF2_MCS_STATE_H

#include "BALL/STRUCTURE/VFLIB/state.h"
#include "BALL/STRUCTURE/VFLIB/argraph.h"

namespace BALL {
namespace VFLib {

class VF2MCSState: public State
  { typedef ARGraph_impl Graph;

    private:
      int core_len, orig_core_len;
      int added_node1;
      int exhaust_len1;
      node_id *core_1;
      node_id *core_2;
      node_id *exhausted_1;

      node_id *order;

      Graph *g1, *g2;
      int n1, n2;

    int *longest_known_mcs;
    long *share_count;

    public:
      VF2MCSState(Graph *g1, Graph *g2, bool sortNodes=false);
      VF2MCSState(const VF2MCSState &state);
      ~VF2MCSState();
      Graph *GetGraph1() { return g1; }
      Graph *GetGraph2() { return g2; }
      bool NextPair(node_id *pn1, node_id *pn2,
                    node_id prev_n1=NULL_NODE, node_id prev_n2=NULL_NODE);
      bool IsFeasiblePair(node_id n1, node_id n2);
      void AddPair(node_id n1, node_id n2);
      bool IsGoal() { return core_len==n1 ; }
      bool IsDead();
      int CoreLen() { return core_len; }
      void GetCoreSet(node_id c1[], node_id c2[]);
      State *Clone();

      virtual void BackTrack();
  };

}
}
#endif // VF2_MCS_STATE_H
