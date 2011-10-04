// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//

#ifndef BALL_STRUCTURE_FRAGMENTDB_ATOMBIJECTION_H
#define BALL_STRUCTURE_FRAGMENTDB_ATOMBIJECTION_H

#include<boost/bimap.hpp>
#include<boost/bimap/unordered_set_of.hpp>

#include<BALL/DATATYPE/GRAPH/VFLIB/state.h>
#include<BALL/DATATYPE/GRAPH/molecularGraph.h>

namespace BALL
{
	class Atom;

	typedef boost::bimaps::bimap<
				boost::bimaps::unordered_set_of<Atom*>,
				boost::bimaps::unordered_set_of<Atom*>
			> AtomBiMapBase;

	/**	Atom bijection.
			This class implements a mapping of two sets of atoms onto each other.
			Since it is a BiMap, it allows for efficient retrieval of mapping 
			partners when only one partner is known.
			If you don't need that functionality, you can use an
			\link AtomBijection AtomBijection
			which is less complex. Behaviour similar to AtomBijection can be had
			when accessing the "left" projection of this bimap, i.e.
			\code
				AtomBiMap map;
				AtomBijection bijection;

				map.left.begin(); // behaves like 
				bijection.begin();
			\endcode
			in that it returns an iterator over something that behaves like an
			std::pair<Atom*,Atom*>.
			\p
			Adding Atoms to the relation works similar to AtomBijection:
			\code
				Atom* atom1 = ...;
				Atom* atom2 = ...;
	
				// Create an empty bijection
				AtomBiMap map;

				// Map atom1 onto atom2.
				bijection.insert(AtomBijection::AtomPair(atom1, atom2));
			\endcode
			\p
			The class offers the full boost::bimaps::bimap interface.
	\ingroup StructureMapping
	*/
	class AtomBiMap : 
		public AtomBiMapBase
	{
		public:
		/** @name Type definitions
		 */
		//@{
		typedef AtomBiMap::value_type AtomPair;
		typedef std::vector<Atom*> RemainderList;
		//@}

		/**	@name	Constructors and Destructors
		*/
		//@{
		
		/**	Default constructor
		*/
		AtomBiMap() {}
		
		/**	Construct a trivial bijection between to atom containers.
		    \p
		    This corresponds to calling assignTrivial after default
		    construction
		*/
		AtomBiMap(AtomContainer& A, AtomContainer& B);

		/**	@name Bijection construction */
		//@{
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
		Size assignTrivial(AtomContainer& A, AtomContainer& B);

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
		Size assignByName(AtomContainer& A, AtomContainer& B);

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
		Size assignCAlphaAtoms(AtomContainer& A, AtomContainer& B);


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
		Size assignBackboneAtoms(AtomContainer& A, AtomContainer& B);

		/** Assign the mapping from a previously computed Morphism.
				\p
				\return The number of atom pairs mapped
		*/
		Size assignFromMorphism(const VFLib::State<MolecularGraph,MolecularGraph>&);
		//@}

		/**	@name Accessors */
		//@{
		///	Calculate the root mean squared deviation of the mapped atoms.
		double calculateRMSD() const;

		const RemainderList& unmappedAtomsFromLeft() const;
		const RemainderList& unmappedAtomsFromRight() const;
		//@}

		Atom* leftFromRight(Atom*);
		Atom* rightFromLeft(Atom*);

		private:
		RemainderList unmapped_A_;
		RemainderList unmapped_B_;

	};

}

#endif // ATOMBIJECTION_H
