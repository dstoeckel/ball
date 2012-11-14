#ifndef BALL_DATATYPE_GRAPH_MOLECULARGRAPH_H
#define BALL_DATATYPE_GRAPH_MOLECULARGRAPH_H

#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>

#include <list>

#ifndef BALL_COMMON_GLOBAL_H
#	include <BALL/COMMON/global.h>
#endif

#ifndef BALL_DATATYPE_HASHSET_H
#	include <BALL/DATATYPE/hashSet.h>
#endif

#ifndef BALL_DATATYPE_GRAPH_GRAPHALGORITHMS_H
# include <BALL/DATATYPE/GRAPH/graphAlgorithms.h>
#endif

namespace BALL
{
	//Forward declarations
	class Atom;
	class Bond;
	class AtomContainer;

	//Define a properly adjusted boost graph type
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
		boost::property<boost::vertex_atom_ptr_t, Atom*>,
		boost::property<boost::edge_bond_ptr_t, Bond*>
	> MolecularGraphBase;

	/**
	 * This class provides a boost graph library representation of
	 * the molecular graph. It takes an AtomContainer as an input
	 * and creates an undirected graph modelling the atom connectivity.
	 * Bonds pointing into another atom container will be omited.
	 * The vertices of the graph are labeled with pointers to the original
	 * atoms; The edges are labeled with pointers to the original bonds.
	 *
	 * These properties can be obtained by using the following code:
	 * @code
	 *   AtomPtrMap atom_ptrs = get(boost::vertex_atom_ptr, graph);
	 *   Atom* atom = boost::get(atom_ptrs, vertex_handle);
	 * @endcode
	 */
	class BALL_EXPORT MolecularGraph
		: public MolecularGraphBase
	{
		public:
			///@TODO do something useful with this
			enum ExportOptions { INCLUDE_ALL = -1 };

			typedef boost::graph_traits<MolecularGraph>::vertex_descriptor Vertex;
			typedef boost::graph_traits<MolecularGraph>::edge_descriptor   Edge;

			typedef boost::graph_traits<MolecularGraph>::vertex_iterator   VertexIterator;
			typedef boost::graph_traits<MolecularGraph>::edge_iterator     EdgeIterator;
			typedef boost::graph_traits<MolecularGraph>::out_edge_iterator OutEdgeIterator;
			typedef boost::graph_traits<MolecularGraph>::in_edge_iterator  InEdgeIterator;

			typedef boost::property_map<MolecularGraphBase, boost::vertex_atom_ptr_t>::type AtomPtrMap;
			typedef boost::property_map<MolecularGraphBase, boost::edge_bond_ptr_t>::type   BondPtrMap;
			typedef boost::property_map<MolecularGraphBase, boost::vertex_atom_ptr_t>::const_type ConstAtomPtrMap;
			typedef boost::property_map<MolecularGraphBase, boost::edge_bond_ptr_t>::const_type   ConstBondPtrMap;

			typedef GRAPH::GraphTraits<MolecularGraph>::EditableGraph EditableGraph;

			MolecularGraph(AtomContainer& ac, ExportOptions opt = INCLUDE_ALL);
			MolecularGraph(HashSet<Atom*>&, HashSet<Bond*>&);

			const Edge&   getEdge  (Bond* bond) const;
			const Vertex& getVertex(Atom* atom) const;

			void editableCopy(EditableGraph& eg);

		private:
			std::map<Bond*, Edge>   bond_to_edge_;
			std::map<Atom*, Vertex> atom_to_vertex_;
	};

	typedef GRAPH::GraphTraits<MolecularGraph> MolecularGraphTraits;
	typedef MolecularGraph::EditableGraph EditableMolecularGraph;
}

namespace boost {
 // FIXME: not sure if this works as intended.
 size_t hash_value(const detail::edge_desc_impl<undirected_tag, long unsigned int>& edge);
}

#endif //BALL_DATATYPE_GRAPH_MOLECULARGRAPH_H
