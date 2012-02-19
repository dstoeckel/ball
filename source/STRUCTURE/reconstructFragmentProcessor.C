// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//

#include <BALL/STRUCTURE/reconstructFragmentProcessor.h>

#include <list>
#include <vector>

#include <BALL/KERNEL/PTE.h>
#include <BALL/SYSTEM/path.h>
#include <BALL/KERNEL/bond.h>
#include <BALL/KERNEL/forEach.h>
#include <BALL/MATHS/matrix44.h>
#include <BALL/DATATYPE/stringHashMap.h>
#include <BALL/STRUCTURE/fragmentDB.h>
#include <BALL/STRUCTURE/structureMapper.h>
#include <BALL/STRUCTURE/atomBiMap.h>
	
using namespace std;

//#define BALL_DEBUG_RECONSTRUCTFRAGMENTPROCESSOR

#ifdef BALL_DEBUG_RECONSTRUCTFRAGMENTPROCESSOR
#	define DEBUG(a) Log.info() << "ReconstructFragmentProcessor: " << a << std::endl;
#else
#	define DEBUG(a)
#endif

namespace BALL 
{

	/////////////////////////////////////////////////////////////////
	//		ReconstructFragmentProcessor												     //
	/////////////////////////////////////////////////////////////////	

	void ReconstructFragmentProcessor::setFragmentDB(const FragmentDB& db)
	{
		fragment_db_ = &const_cast<FragmentDB&>(db);
	}

	const FragmentDB* ReconstructFragmentProcessor::getFragmentDB() const
	{
		return fragment_db_;
	}

	/**	Identify two reference atoms.
			Performs a breadth-first search for two additional heavy atoms
			starting from the center atom. These atoms are used as 
			anchor points for attaching the next atom.
	*/
	Triple<bool, const Atom*, const Atom*> 
	ReconstructFragmentProcessor::getTwoReferenceAtoms
		(const Atom& ref_center_atom, const HashSet<Atom*>& allowed)
		
	{
		Triple<bool, const Atom*, const Atom*> result(false, 0, 0);

		// a hash set to remember all those atoms we have already visited
		list<const Atom*> atom_list;
		atom_list.push_back(&ref_center_atom);

		// abort if we found the three first atoms (beyond the center atom)
		// or we are running out of fresh atoms
		list<const Atom*>::iterator current(atom_list.begin());
		while ((atom_list.size() < 3) && (current != atom_list.end()))
		{
			Atom::BondConstIterator bond((*current)->beginBond());
			for (; +bond; ++bond)
			{
				Atom* next_atom = bond->getPartner(**current);
				if (allowed.has(next_atom) 
						&& (find(atom_list.begin(), atom_list.end(), next_atom) == atom_list.end()))
				{
					atom_list.push_back(next_atom);
					if (atom_list.size() > 2)
					{
						bond = (*current)->endBond();
						break;
					}
				}
			}

			// try the bonds of the next atom in the list
			current++;
		}

		// copy the two  reference atoms to the result 
		// (omit the first atom, which is the center atom!)
		result.first = (atom_list.size() == 3);
		current = atom_list.begin();
		current++;
		if (current != atom_list.end())
		{
			result.second = *current;
			current++;
		}
		if (current != atom_list.end())
		{
			result.third  = *current;
		}

		return result;
	}

	// start function of ReconstructFragmentProcessor
	// nothing important is done here
	bool ReconstructFragmentProcessor::start()
	{
		inserted_atoms_.clear();

		if (fragment_db_ == 0)
		{
			Log.error() << "ReconstructFragmentProcessor: no FragmentDB defined. "
									<< "Use setFragmentDB() to associate a fragment database." << std::endl;
			return false;
		}

		return true;
	}
	
	// Processor finish method
	bool ReconstructFragmentProcessor::finish()
	{
		return true;
	}

	// Processor application method
	Processor::Result ReconstructFragmentProcessor::operator () (Fragment& object)
	{
		// abort if the object is not a residue (only residues are 
		// contained in the fragment DB)
		if (!RTTI::isKindOf<Residue>(object))
		{
			return Processor::CONTINUE;
		}

		// cast the object to a residue
		Residue& residue = dynamic_cast<Residue&>(object);

		// get the reference fragment from the fragment DB
		boost::shared_ptr<Residue> reference_fragment = fragment_db_->findReferenceFragment(residue);

		// complain if no reference fragment could be found
		if (!reference_fragment)
		{
			Log.warn() << "ReconstructFragmentProcessor: no reference fragment found for " 
							   << residue.getName() << ":" << residue.getID() << std::endl;
			return Processor::CONTINUE;
		}

		// Reconstruct the atoms and count the number of new atoms.
		// number_of_inserted_atoms_ += reconstructFragment(residue, *reference_fragment);
		list<Atom*> inserted_atoms;
		list<Atom*>::iterator it;

		inserted_atoms = reconstructFragment(residue, *reference_fragment);

		for (it = inserted_atoms.begin(); it != inserted_atoms.end(); ++it)
		{
			inserted_atoms_.push_back(*it);
		}

		// Next residue.
		return Processor::CONTINUE;
	}

	ReconstructFragmentProcessor::ReconstructFragmentProcessor(const FragmentDB& db)
		:	fragment_db_(&db),
			inserted_atoms_()
	{
	}

	// copy constructor	
	ReconstructFragmentProcessor::ReconstructFragmentProcessor(const ReconstructFragmentProcessor& rfp)
		:	UnaryProcessor<Fragment>(rfp),
			fragment_db_(rfp.fragment_db_),
			inserted_atoms_(rfp.inserted_atoms_)
	{
	}

	// default constructor	
	ReconstructFragmentProcessor::ReconstructFragmentProcessor()
		:	fragment_db_(0),
			inserted_atoms_()
	{
	}

	// destructor	
	ReconstructFragmentProcessor::~ReconstructFragmentProcessor()
	{
		fragment_db_ = 0;
	}

	list<Atom*>& ReconstructFragmentProcessor::getInsertedAtoms()
	{
		return inserted_atoms_;
	}

	// return the number of inserted atoms
	Size ReconstructFragmentProcessor::getNumberOfInsertedAtoms() const
	{
		return inserted_atoms_.size();
	}

	list<Atom*> ReconstructFragmentProcessor::reconstructFragment
		(Fragment& fragment, const Fragment& tplate)
	{
		// We count the number of atoms created.
		list<Atom*> inserted_atoms;
		HashSet<Atom*> transformed;

		AtomBiMap mapping;
		// have to use const_cast here, as AtomBiMap is more general,
		// but we won't do anything bad in here, will we?
		mapping.assignByName(fragment, const_cast<Fragment&>(tplate));

		AtomBiMap::iterator template_atom = mapping.begin();
		for (; template_atom != mapping.end(); ++template_atom)
		{
			transformed.insert(template_atom->right);
		}

		AtomBiMap::RemainderList::const_iterator unmapped_atom_from_template
			= mapping.unmappedAtomsFromRight().begin();
		for (; unmapped_atom_from_template
		       != mapping.unmappedAtomsFromRight().end()
		     ; ++unmapped_atom_from_template)
		{
			// FIXME: add option for only adding hydrogens (bug#17)
			Atom* copy_of_unmapped = reinterpret_cast<Atom*>(
				(*unmapped_atom_from_template)->create(false)
			);
			fragment.insert(*copy_of_unmapped);
			mapping.insert(
				AtomBiMap::AtomPair(copy_of_unmapped, *unmapped_atom_from_template)
			);
			inserted_atoms.push_back(copy_of_unmapped);
		}

		DEBUG("Found "<<mapping.size()<<" same atoms, "
		      << mapping.unmappedAtomsFromLeft().size() <<" remain from fragment, "
		      << mapping.unmappedAtomsFromRight().size() <<" remain from template."
		     )

		// TODO: IDEA: if the remainders are of equal length, maybe try harder.

		// We've now made sure that all atoms of the tplate exist in the 
		// reconstructed residue as well (careful, not the other way round!)
		// we can now start to adjust the atom coordinates.

		// If no atoms were in common, there's not much we can do...
		// Trivial solution: no atoms are actually matched to each 
		// other, so we just leave the coordinates the way they
		// are (copy of the tpl coordinates) and return.
		if (!transformed.isEmpty())
		{
			// Otherwise, we start adjusting coordinates
			// We use a hash set for BFS
			HashSet<Atom*> visited;

			list<Atom*> stack;
			stack.push_back(*transformed.begin());
			while (!stack.empty())
			{
				Atom* template_atom = stack.front();
				stack.pop_front();
				visited.insert(template_atom);

				DEBUG("center is " << (void*)template_atom << " (" << template_atom->getFullName() << ") visited = "
				      << (visited.has(template_atom)) << " transformed = " << transformed.has(template_atom)
				      << " @ " << template_atom->getPosition())
				DEBUG("residue atom is @ " << mapping.leftFromRight(template_atom)->getPosition()  << " (dist = "
				      << mapping.leftFromRight(template_atom)->getPosition().getDistance(template_atom->getPosition()) << ")")

				Atom::BondConstIterator bond;
				for (bond = template_atom->beginBond(); +bond; ++ bond)
				{
					Atom* translation_template = bond->getPartner(*template_atom);
					DEBUG("examining "
					      << (void*)translation_template
					      << " (" << translation_template->getFullName() << ")"
					      << " visited = " 	<< (visited.has(translation_template))
					      << " transformed = " << transformed.has(translation_template)
					     )

					if (!visited.has(translation_template))
					{
						stack.push_back(translation_template);
						visited.insert(translation_template);
						if (!transformed.has(translation_template))
						{
							DEBUG("searching reference atoms for "<< translation_template->getFullName())
							Triple<bool, const Atom*, const Atom*> correctly_positioned_partners;
							correctly_positioned_partners = getTwoReferenceAtoms(*template_atom, transformed);
							Atom* partner_a = const_cast<Atom*>(correctly_positioned_partners.second);
							Atom* partner_b = const_cast<Atom*>(correctly_positioned_partners.third);
							DEBUG("reference atoms: "
							      << (partner_a == 0 ? String("-") : partner_a->getFullName()) << " / "
							      << (correctly_positioned_partners.third == 0 ? String("-") : correctly_positioned_partners.third->getFullName())
							     )

							Matrix4x4 transform_to_match_neighbors;
							if (correctly_positioned_partners.first == true)
							{
								// we can map all three atoms, great!
								DEBUG("mapping three atoms: " << mapping.leftFromRight(template_atom)->getFullName() << "/"
								      << mapping.leftFromRight(partner_a)->getFullName() << "/"
								      << mapping.leftFromRight(partner_b)->getFullName())
								DEBUG("onto:                " << template_atom->getFullName() << "/"
								      << partner_a->getFullName() << "/"
								      << partner_b->getFullName())
								DEBUG("from: " << mapping.leftFromRight(template_atom)->getPosition() << "/"
								      << mapping.leftFromRight(partner_a)->getPosition() << "/"
								      << mapping.leftFromRight(partner_b)->getPosition())
								DEBUG("to:   " << template_atom->getPosition() << "/"
								      << partner_a->getPosition() << "/"
								      << partner_b->getPosition())

								transform_to_match_neighbors = StructureMapper::matchPoints(
								       // map these three known positions:
								       template_atom->getPosition(),
								       partner_a->getPosition(),
								       partner_b->getPosition(),
								       // onto these three known positions:
								       mapping.leftFromRight(template_atom)->getPosition(),
								       mapping.leftFromRight(partner_a)->getPosition(),
								       mapping.leftFromRight(partner_b)->getPosition());

								DEBUG("found two reference atoms, mapped with T =\n" << transform_to_match_neighbors)
							}
							else 
							{
								// We could map the two center atoms only, which corresponds to 
								// a simple translation by the difference of the two atom positions.
								transform_to_match_neighbors.setIdentity();
								transform_to_match_neighbors.setTranslation(mapping.leftFromRight(template_atom)->getPosition() - template_atom->getPosition());
								DEBUG("translating by " << transform_to_match_neighbors)
							}

							Atom* fragment_translatee = mapping.leftFromRight(translation_template);
							// Apply the transformation to the coordinates of the atom we're interested in
							fragment_translatee->setPosition(transform_to_match_neighbors * fragment_translatee->getPosition());

							// Remember that we already took care of that guy.
							transformed.insert(translation_template);
							DEBUG(translation_template->getFullName() << " is transformed: " << fragment_translatee->getPosition()
							      << "/" << translation_template->getPosition())
							DEBUG("distance = " << fragment_translatee->getPosition().getDistance(translation_template->getPosition()))
						}
					}
				}
			}
		}
		
		// Return the number of atoms created.
		// return number_of_inserted_atoms;
		return inserted_atoms;
	}

} // namespace BALL

