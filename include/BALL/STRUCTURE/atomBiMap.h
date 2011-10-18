// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//

#ifndef BALL_STRUCTURE_FRAGMENTDB_ATOMBIJECTION_H
#define BALL_STRUCTURE_FRAGMENTDB_ATOMBIJECTION_H

#include<boost/bimap.hpp>
#include<boost/bimap/unordered_set_of.hpp>

#include<BALL/DATATYPE/GRAPH/VFLIB/state.h>
#include<BALL/DATATYPE/GRAPH/molecularGraph.h>
#include<BALL/STRUCTURE/atomPairContainer.h>

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
		public AtomBiMapBase, public AtomPairContainer<AtomBiMapBase::value_type>
	{
		public:
		/** @name Type definitions
		 */
		//@{
		typedef AtomBiMap::value_type AtomPair;
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


		/**	@name Accessors */
		//@{
		///	Calculate the root mean squared deviation of the mapped atoms.
		double calculateRMSD() const;

		Atom* leftFromRight(Atom*);
		Atom* rightFromLeft(Atom*);
		//@}

		/** @name AtomPairContainer interface */
		//@{
		virtual void addPair(const AtomPair& pair) { insert(pair); }
		/** Clear the container */
		virtual void clearContainer() { clear(); }
		/** Get the size of the Container */
		virtual Size containerSize() { return size(); }
		//@}

	};

}

#endif // ATOMBIJECTION_H
