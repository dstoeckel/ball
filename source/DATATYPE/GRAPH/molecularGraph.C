#include <BALL/DATATYPE/GRAPH/molecularGraph.h>

#include <BALL/KERNEL/atom.h>
#include <BALL/KERNEL/bond.h>
#include <BALL/KERNEL/atomContainer.h>

#include <BALL/COMMON/exception.h>

namespace BALL
{
	class BondVisitor : public UnaryProcessor<Bond>
	{
		public:
			BondVisitor(MolecularGraph& molgraph, std::map<Bond*, MolecularGraph::Edge>& bte)
				: molgraph_(molgraph),
				  bte_(bte),
				  bond_ptrs_(boost::get(boost::edge_bond_ptr, molgraph))
			{
			}

			Processor::Result operator()(Bond& bond)
			{
				MolecularGraph::Vertex a = molgraph_.getVertex(bond.getFirstAtom());
				MolecularGraph::Vertex b = molgraph_.getVertex(bond.getSecondAtom());

				MolecularGraph::Edge edge = boost::add_edge(a, b, molgraph_).first;
				boost::put(bond_ptrs_, edge, &bond);
				bte_.insert(std::make_pair(&bond, edge));

				return Processor::CONTINUE;
			}

		private:
			MolecularGraph& molgraph_;
			std::map<Bond*, MolecularGraph::Edge>& bte_;
			MolecularGraph::BondPtrMap bond_ptrs_;
	};

	MolecularGraph::MolecularGraph(AtomContainer& ac, MolecularGraph::ExportOptions /*opts*/)
		: MolecularGraphBase(ac.countAtoms())
	{
		// Create a bijective mapping between Vertices and Atom*
		AtomPtrMap atom_ptrs = get(boost::vertex_atom_ptr, *this);

		AtomIterator ait = ac.beginAtom();
		VertexIterator vi, vi_end;
		boost::tie(vi, vi_end) = boost::vertices(*this);

		for(; +ait; ++vi, ++ait)
		{
			atom_to_vertex_.insert(std::make_pair(&*ait,*vi));
			boost::put(atom_ptrs, *vi, &*ait);
		}

		//Create all bonds
		BondVisitor visitor(*this, bond_to_edge_);
		ac.applyIntraBond(visitor);
	}

	MolecularGraph::MolecularGraph(HashSet<Atom*> &atoms, HashSet<Bond*> &bonds)
		: MolecularGraphBase(atoms.size())
	{
		AtomPtrMap vert2atom = boost::get(boost::vertex_atom_ptr, *this);
		BondPtrMap edge2bond = boost::get(boost::edge_bond_ptr, *this);

		VertexIterator vi, vi_end;
		boost::tie(vi, vi_end) = boost::vertices(*this);

		HashSet<Atom*>::Iterator ait = atoms.begin();
		for(; ait != atoms.end(); ++ait, ++vi) {
			atom_to_vertex_.insert(std::make_pair(*ait, *vi));
			boost::put(vert2atom, *vi, *ait);
		}

		HashSet<Bond*>::Iterator bond = bonds.begin();
		for (; bond != bonds.end(); ++bond) {
			Bond* theBond = *bond;
			std::map<Atom*, Vertex>::const_iterator v_first, v_second;
			v_first = atom_to_vertex_.find(theBond->getFirstAtom());
			v_second = atom_to_vertex_.find(theBond->getSecondAtom());
			if ((v_first != atom_to_vertex_.end()) &&
			    (v_second != atom_to_vertex_.end()))
			{
				Edge e = boost::add_edge(v_first->second, v_second->second, *this).first;
				boost::put(edge2bond, e, theBond);
				bond_to_edge_.insert(std::make_pair(theBond, e));
			}
		}
	}

	const MolecularGraph::Edge& MolecularGraph::getEdge(Bond* bond) const
	{
		std::map<Bond*, Edge>::const_iterator it = bond_to_edge_.find(bond);

		if(it == bond_to_edge_.end())
		{
			throw Exception::InvalidArgument(__FILE__, __LINE__, "An unknown Bond pointer was specified");
		}

		return it->second;
	}

	const MolecularGraph::Vertex& MolecularGraph::getVertex(Atom* atom) const
	{
		std::map<Atom*, Vertex>::const_iterator it = atom_to_vertex_.find(atom);

		if(it == atom_to_vertex_.end())
		{
			throw Exception::InvalidArgument(__FILE__, __LINE__, "An unknown Atom pointer was specified");
		}

		return it->second;
	}

	void MolecularGraph::editableCopy(EditableGraph& eg)
	{
		boost::copy_graph(*dynamic_cast<MolecularGraphBase*>(this), eg, 
				vertex_copy(GRAPH::makeEditableVertexCopier(*this, eg)).edge_copy(GRAPH::makeEditableEdgeCopier(*this, eg)));
	}
}

namespace boost {
 // FIXME: not sure if this works as intended.
  size_t hash_value(const detail::edge_desc_impl<undirected_tag, long unsigned int>& edge) {
	 return hash_value(edge.m_source)+hash_value(edge.m_target)+hash_value(edge.get_property());
 }
}
