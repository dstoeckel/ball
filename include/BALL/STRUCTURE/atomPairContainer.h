// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//

#ifndef BALL_STRUCTURE_ATOMPAIRCONTAINER_H
#define BALL_STRUCTURE_ATOMPAIRCONTAINER_H

#ifndef BALL_COMMON_H
#	include <BALL/common.h>
#endif

#ifndef BALL_KERNEL_ATOMCONTAINER_H
#	include <BALL/KERNEL/atomContainer.h>
#endif

#ifndef BALL_KERNEL_RESIDUE_H
#	include <BALL/KERNEL/residue.h>
#endif

#ifndef BALL_KERNEL_EXTRACTORS_H
# include <BALL/KERNEL/extractors.h>
#endif

#ifdef BALL_USE_NEW_MOLECULAR_GRAPH
#ifndef BALL_DATATYPE_MOLECULARGRAPH_H
# include <BALL/DATATYPE/GRAPH/molecularGraph.h>
#endif
#endif


#ifndef BALL_DATATYPE_GRAPH_VFLIB_STATE_H
# include <BALL/DATATYPE/GRAPH/VFLIB/state.h>
#endif

namespace BALL {
	/** A generic interface for filling containers of pairs of atoms.
	    This class centralizes the algorithms that are being used in the
	    construction of AtomBijections and AtomBiMaps, to avoid maintaining
	    the code in several different locations.
	    @p
	    The Template parameter, AtomContainerPair is the type of pair that
	    goes inside the container that inherits from this, and has to have
	    a constructor AtomContainerPair(Atom*,Atom*).
	    @p
	    The protected pure virtual functions, addPair, clearContainer and
	    containerSize have to be reimplemented inside the inheriting class,
	    to abstract away from the actual access to the container.
	*/
	template <typename AtomContainerPair>
	class AtomPairContainer {

		public:
		/** @name Type Definitions */
		//@{
		typedef std::vector<Atom*> RemainderList;
#ifdef BALL_USE_NEW_MOLECULAR_GRAPH
		typedef boost::property_map<MolecularGraph, boost::vertex_index_t>::type VertexIndexMap;
#endif
		//@}

		/**	@name Bijection construction */
		//@{
		/**	Assign the atom pairs through a name matching.
				This method creates a mapping based on the atom names.
				If the atom is contained in a Residue/Protein, the name consists
				of the fully qualified name (<chain>:<residue name>:<residue id>:<atom name>).
				If no pair of atoms could by matched this way, it will try to match by
				atom names only (not considering residues, chains or the like).
				\p
				The method constructs a hash map for all atom names, so run time is linear
				in the number of atoms.
				\p
				The number of atoms mapped is returned.
		*/
		Size assignByName(AtomContainer& A, AtomContainer& B, bool limit_to_selection = false)
		{
			// Clear old bijection.
			clearContainer();
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
					if (   !limit_to_selection
							|| (ai->isSelected() || A_names[ai->getFullName(Atom::ADD_VARIANT_EXTENSIONS_AND_ID)]->isSelected()))
					{
						// We found two matching atoms. Remember them.
						addPair(AtomContainerPair(A_names[ai->getFullName(Atom::ADD_VARIANT_EXTENSIONS_AND_ID)], &*ai));
						// Throw away the hash map entry in order to avoid
						// 1:n mappings.
						A_names.erase(ai->getFullName(Atom::ADD_VARIANT_EXTENSIONS_AND_ID));
					}
					else
					{
						unmapped_B_.push_back(&*ai);
					}
				}
				else
				{
					// There's no matching atom in A for this one.
					unmapped_B_.push_back(&*ai);
				}
			}

			// Check whether we could map anything. 
			// If not, try to map by atom name alone.
			if (containerSize() == 0)
			{
				// Next stage: try to map by atom name only.
				A_names.clear();
				for (AtomIterator ai = A.beginAtom(); +ai; ++ai)
				{
					A_names.insert(std::pair<String, Atom*>(ai->getName(), &*ai));
				}
				clearContainer();
				unmapped_B_.clear();

				for (AtomIterator ai = B.beginAtom(); +ai; ++ai)
				{
					if (A_names.has(ai->getName()))
					{
					if (   !limit_to_selection
							|| (ai->isSelected() || A_names[ai->getName()]->isSelected()))
					{
						// We found two matching atoms. Remember them.
						addPair(AtomContainerPair(A_names[ai->getName()], &*ai));
						// Throw away the hash map entry in order to avoid
						// 1:n mappings.
						A_names.erase(ai->getName());
						}
						else
						{
							unmapped_B_.push_back(&*ai);
						}
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

			return containerSize();
		}

		/** Assign all atoms in the two atom containers in order.
				Construct a simple bijection mapping the atoms of the two
				atom containers onto each other. The mapping iterates
				over the atoms and stops assigning pairs of atoms as soon
				as the smalles of the two atom sets is fully assigned. 
				The larger of the two atom container can thus contain
				unassigned atoms. No checking with respect to atom names,
				elements or the like are being made.
				\p
				This trivial bijection is useful, if the two atom containers
				correspond to exactly the same structure (i.e. they
				just differ in their conformations).
				Care must be taken that the order of atoms is correct.
				Beware of adding hydrogens, which might mess up atom
				order in some cases.
				\p
				The number of atoms mapped is returned.
		*/
		Size assignTrivial(AtomContainer& A, AtomContainer& B, bool limit_to_selection = false)
		{
			// Delete old bijection.
			clearContainer();
			unmapped_A_.clear();
			unmapped_B_.clear();

			// Map in order -- first atom of A onto
			// first atom of B and so on.
			AtomIterator ai(A.beginAtom());
			AtomIterator bi(B.beginAtom());
			for (; +ai && +bi; ++ai, ++bi)
			{
				if (   !limit_to_selection
						|| (ai->isSelected() || bi->isSelected()))
				{
					addPair(AtomContainerPair(&*ai, &*bi));
				}
				else
				{
					unmapped_A_.push_back(&*ai);
					unmapped_B_.push_back(&*bi);
				}
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

			return containerSize();
		}

		/** Assign the C-alpha atoms ordered by sequence.
				This method iterated over all residues and assigns
				the C-alpha atoms (i.e. all atoms named "CA" in a
				residue with the property AMINO_ACID) of the two proteins
				in the order they are traversed. The size of the mapping
				corresponds to the minimum of the number of C-alpha atoms
				of both atom containers.
				\note For efficiency reasons, an AtomBiMap assigned this way
				does not allow for retrieval of #unmappedAtomsFromA() or
				#unmappedAtomsFromB().
				\p
				\return The number of atom pairs mapped
		*/
		Size assignCAlphaAtoms(AtomContainer& A, AtomContainer& B, bool limit_to_selection = false)
		{
			// Delete old bijection.
			clearContainer();

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
				if (caa != 0 && cab != 0
					&&  (   !limit_to_selection
					     || (caa->isSelected() || cab->isSelected())))
				{
					addPair(AtomContainerPair(caa, cab));
				}
			}
			//
			return containerSize();
		}

		/** Assign the backbone atoms ordered by sequence.
				This method iterated over all residues and assigns
				the backbone atoms (i.e. the first(!) atoms named "CA", "C",
				"N", "H", or "O" in every residue with the property AMINO_ACID)
				of the two proteins
				in the order they are traversed. The mapping terminates,
				if the traversal of the residues in one of the two atom containers
				terminates.
				\note For efficiency reasons, an AtomBiMap assigned this way
				does not allow for retrieval of #unmappedAtomsFromA() or
				#unmappedAtomsFromB().
				\p
				\return The number of atom pairs mapped
		*/
		Size assignBackboneAtoms(AtomContainer& A, AtomContainer& B, bool limit_to_selection = false)
		{
			// Delete old bijection.
			clearContainer();
			
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
				if (a != 0 && b != 0
				  &&  (   !limit_to_selection
					     || (a->isSelected() || b->isSelected())))
				{
					addPair(AtomContainerPair(a, b));
				}
				a = (*ita)->getAtom("C");
				b = (*itb)->getAtom("C");
				if (a != 0 && b != 0
				  &&  (   !limit_to_selection
					     || (a->isSelected() || b->isSelected())))
				{
					addPair(AtomContainerPair(a, b));
				}
				a = (*ita)->getAtom("N");
				b = (*itb)->getAtom("N");
				if (a != 0 && b != 0
				  &&  (   !limit_to_selection
					     || (a->isSelected() || b->isSelected())))
				{
					addPair(AtomContainerPair(a, b));
				}
				a = (*ita)->getAtom("O");
				b = (*itb)->getAtom("O");
				if (a != 0 && b != 0
				  &&  (   !limit_to_selection
					     || (a->isSelected() || b->isSelected())))
				{
					addPair(AtomContainerPair(a, b));
				}
				a = (*ita)->getAtom("H");
				b = (*itb)->getAtom("H");
				if (a != 0 && b != 0
				  &&  (   !limit_to_selection
					     || (a->isSelected() || b->isSelected())))
				{
					addPair(AtomContainerPair(a, b));
				}
			}

			return containerSize();
		}

		/** Assign the atom pairs through a name matching and based on the 
		    property "ATOMBIJECTION_RMSD_SELECTION".
				@see assignTrivial()
				This is restriction is useful, if the focus of investigation
				is limited to e.g. a binding pocket.
		*/
		Size assignAtomsByProperty(AtomContainer& A, AtomContainer& B, bool limit_to_selection = false)
		{
			// Delete old bijection.
			clearContainer();
			// Map in order -- first atom of A onto
			// first atom of B and so on.
			AtomIterator ai(A.beginAtom());
			AtomIterator bi(B.beginAtom());
			for (; +ai && +bi; ++ai, ++bi)
			{
				if ( ai->hasProperty("ATOMBIJECTION_RMSD_SELECTION") || bi->hasProperty("ATOMBIJECTION_RMSD_SELECTION"))
				{
					if(!limit_to_selection || ai->isSelected() || bi->isSelected())
					{
						addPair(AtomContainerPair(&*ai, &*bi));
					}
				}
			}
			return containerSize();
		}

#ifdef BALL_USE_NEW_MOLECULAR_GRAPH
		/** Assign the mapping from a previously computed Morphism.
				\p
				\return The number of atom pairs mapped
		*/
		Size assignFromMorphism(const VFLib::State<MolecularGraph,MolecularGraph>& state)
		{
			clearContainer();

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
				addPair(AtomContainerPair(boost::get(queryVertexLabels,boost::get(queryVertexIndices,queryCore[pair_no])),
												boost::get(targetVertexLabels,boost::get(targetVertexIndices,targetCore[pair_no]))));
			}

			// TODO: retrieve unassigned Atoms.

			delete[] queryCore;
			delete[] targetCore;

			return containerSize();
		}
#endif
		//@}

		const RemainderList& unmappedAtomsFromLeft() const
		{
			return unmapped_A_;
		}

		const RemainderList& unmappedAtomsFromRight() const
		{
			return unmapped_B_;
		}

		protected:
		/**	@name Interface requirements */
		//@{
		/** Add a Pair to the Container */
		virtual void addPair(const AtomContainerPair&) = 0;
		/** Clear the container */
		virtual void clearContainer() = 0;
		/** Get the size of the Container */
		virtual Size containerSize() = 0;
		//@}
		
		RemainderList unmapped_A_;
		RemainderList unmapped_B_;
	};
}

#endif // ATOMPAIRCONTAINER_H
