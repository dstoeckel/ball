
#include <BALL/STRUCTURE/atomBiMap.h>
#include <BALL/MATHS/vector3.h>
#include <BALL/KERNEL/atom.h>
#include <BALL/KERNEL/atomContainer.h>
#include <BALL/KERNEL/extractors.h>
#include <BALL/KERNEL/residue.h>
#include <BALL/DATATYPE/hashGrid.h>

namespace BALL {

	typedef boost::property_map<MolecularGraph, boost::vertex_index_t>::type VertexIndexMap;
	
	Size AtomBiMap::assignFromMorphism(const VFLib::State<MolecularGraph,MolecularGraph>& state)
	{
		clear();

		const MolecularGraph* queryGraph = state.GetGraph1();
		const MolecularGraph* targetGraph = state.GetGraph2();
	
		MolecularGraph::ConstAtomPtrMap queryVertexLabels = boost::get(boost::vertex_atom_ptr, *queryGraph);
		MolecularGraph::ConstAtomPtrMap targetVertexLabels = boost::get(boost::vertex_atom_ptr, *targetGraph);
	
		VertexIndexMap queryVertexIndices = boost::get(boost::vertex_index, *queryGraph);
		VertexIndexMap targetVertexIndices = boost::get(boost::vertex_index, *targetGraph);
	
		VFLib::node_id* queryCore = new VFLib::node_id[state.CoreLen()];
		VFLib::node_id* targetCore = new VFLib::node_id[state.CoreLen()];
	
		state.GetCoreSet(queryCore, targetCore);
		
		for (int pair_no = 0; pair_no < state.CoreLen(); ++pair_no) {
			insert(AtomPair(boost::get(queryVertexLabels,boost::get(queryVertexIndices,queryCore[pair_no])),
											boost::get(targetVertexLabels,boost::get(targetVertexIndices,targetCore[pair_no]))));
		}

		// TODO: retrieve unassigned Atoms.

		delete[] queryCore;
		delete[] targetCore;

		return size();
	}

	/* Calculate the root mean squared deviation */
	double AtomBiMap::calculateRMSD() const
	{
		double sum_of_squares = 0.0;
		
		if (!empty())
		{
			left_const_iterator it = left.begin();

			for(; it != left.end(); ++it) {
				Vector3& r(it->first->getPosition());
				sum_of_squares += r.getSquareDistance(it->second->getPosition());
			}

			// calculate mean square deviation
			sum_of_squares = sqrt(sum_of_squares / (double)size());
		}

		// return RMSD
		return sum_of_squares;
	}

	Size AtomBiMap::assignByName(AtomContainer& A, AtomContainer& B)
	{
		// Clear old bijection.
		clear();
		unmapped_A_.clear();
		unmapped_B_.clear();

		// Remember the names of A and their atom pointers.
		StringHashMap<Atom*> A_names;
		for (AtomIterator ai = A.beginAtom(); +ai; ++ai)
		{
			A_names.insert(std::pair<String, Atom*>(ai->getFullName(Atom::ADD_VARIANT_EXTENSIONS_AND_ID), &*ai));
		}
		
		// Iterate over all atoms of B and try to find an 
		// atom in A identical names.
		for (AtomIterator ai = B.beginAtom(); +ai; ++ai)
		{
			if (A_names.has(ai->getFullName(Atom::ADD_VARIANT_EXTENSIONS_AND_ID)))
			{
				// We found two matching atoms. Remember them.
				insert(AtomPair(A_names[ai->getFullName(Atom::ADD_VARIANT_EXTENSIONS_AND_ID)], &*ai));

				// Throw away the hash map entry in order to avoid
				// 1:n mappings.
				A_names.erase(ai->getFullName(Atom::ADD_VARIANT_EXTENSIONS_AND_ID));
			}
			else
			{
				// There's no matching atom in A for this one.
				unmapped_B_.push_back(&*ai);
			}
		}

		// Check whether we could map anything. 
		// If not, try to map by atom name alone.
		if (size() == 0)	
		{
			// Next stage: try to map by atom name only.
			A_names.clear();
			for (AtomIterator ai = A.beginAtom(); +ai; ++ai)
			{
				A_names.insert(std::pair<String, Atom*>(ai->getName(), &*ai));
			}
			clear();
			unmapped_B_.clear();

			for (AtomIterator ai = B.beginAtom(); +ai; ++ai)
			{
				if (A_names.has(ai->getName()))
				{
					// We found two matching atoms. Remember them.
					insert(AtomPair(A_names[ai->getName()], &*ai));
					// Throw away the hash map entry in order to avoid
					// 1:n mappings.
					A_names.erase(ai->getName());
				}
				else
				{
					// There's no matching atom in A for this one.
					unmapped_B_.push_back(&*ai);
				}
			}
		}

		// Collect remaining unmapped Atoms from A.
		if (A_names.size() > 0) {
			StringHashMap<Atom*>::Iterator leftover = A_names.begin();
			for (; leftover != A_names.end(); ++leftover)
			{
				unmapped_A_.push_back(leftover->second);
			}
		}

		return size();
	}

	Size AtomBiMap::assignTrivial(AtomContainer& A, AtomContainer& B)	
	{
		// Delete old bijection.
		clear();
		unmapped_A_.clear();
		unmapped_B_.clear();

		// Map in order -- first atom of A onto
		// first atom of B and so on.
		AtomIterator ai(A.beginAtom());
		AtomIterator bi(B.beginAtom());
		for (; +ai && +bi; ++ai, ++bi)
		{
			insert(AtomPair(&*ai, &*bi));
		}
		if (+ai) {
			for (; +ai; ++ai) {
				unmapped_A_.push_back(&*ai);
			}
		}
		else if (+bi)
		{
			for (; +bi; ++bi) {
				unmapped_B_.push_back(&*bi);
			}
		}

		return size();
	}

	Size AtomBiMap::assignCAlphaAtoms(AtomContainer& A, AtomContainer& B)
	{
		// Delete old bijection.
		clear();
		
		// Extract all residues in A and B
		ResidueList rla(residues(A));
		ResidueList rlb(residues(B));

		// Walk over the residues in parallel
		ResidueList::iterator ita(rla.begin());
		ResidueList::iterator itb(rlb.begin());
		for (; ita != rla.end() && itb != rlb.end(); ++ita, ++itb)
		{
			// If the two residues do have an atom named CA, push back the pair.
			Atom* caa = (*ita)->getAtom("CA");
			Atom* cab = (*itb)->getAtom("CA");
			if (caa != 0 && cab != 0)
			{
				insert(AtomPair(caa, cab));
			}
		}
		//
		return size();
	}

	Size AtomBiMap::assignBackboneAtoms(AtomContainer& A, AtomContainer& B)
	{
		// Delete old bijection.
		clear();
		
		// Extract all residues in A and B
		ResidueList rla(residues(A));
		ResidueList rlb(residues(B));

		// Walk over the residues in parallel
		ResidueList::iterator ita(rla.begin());
		ResidueList::iterator itb(rlb.begin());
		for (; ita != rla.end() && itb != rlb.end(); ++ita, ++itb)
		{
			// If the two residues do have backbone atoms (CA, C, N, O, H)
			// map then onto each other.
			Atom* a = (*ita)->getAtom("CA");
			Atom* b = (*itb)->getAtom("CA");
			if (a != 0 && b != 0)
			{
				insert(AtomPair(a, b));
			}
			a = (*ita)->getAtom("C");
			b = (*itb)->getAtom("C");
			if (a != 0 && b != 0)
			{
				insert(AtomPair(a, b));
			}
			a = (*ita)->getAtom("N");
			b = (*itb)->getAtom("N");
			if (a != 0 && b != 0)
			{
				insert(AtomPair(a, b));
			}
			a = (*ita)->getAtom("O");
			b = (*itb)->getAtom("O");
			if (a != 0 && b != 0)
			{
				insert(AtomPair(a, b));
			}
			a = (*ita)->getAtom("H");
			b = (*itb)->getAtom("H");
			if (a != 0 && b != 0)
			{
				insert(AtomPair(a, b));
			}
		}
			
		//
		return size();
	}

	const std::vector<Atom*>& AtomBiMap::unmappedAtomsFromLeft() const
	{
		return unmapped_A_;
	}

	const std::vector<Atom*>& AtomBiMap::unmappedAtomsFromRight() const
	{
		return unmapped_B_;
	}

	Atom* AtomBiMap::leftFromRight(Atom* atom) {
		right_const_iterator it = right.find(atom);
		if (it == right.end())
		{
			return 0;
		}
		else
		{
			return it->second;
		}
	}

	Atom* AtomBiMap::rightFromLeft(Atom* atom) {
		left_const_iterator it = left.find(atom);
		if (it == left.end())
		{
			return 0;
		}
		else
		{
			return it->second;
		}
	}

}
