// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//

#include <BALL/STRUCTURE/fragmentDB.h>

#include <BALL/KERNEL/PTE.h>
#include <BALL/COMMON/limits.h>
#include <BALL/SYSTEM/path.h>
#include <BALL/KERNEL/forEach.h>
#include <BALL/KERNEL/nucleotide.h>

#include <BALL/STRUCTURE/FRAGMENTDB/resourceFileFragmentStorage.h>
#include <BALL/STRUCTURE/FRAGMENTDB/nameMapQuery.h>
	
/*			Things still missing (among others)
				===================================
				- check for unique atom names
				- dynamic import of databases
*/

//#define BALL_DEBUG_FRAGMENTDB

#ifdef BALL_DEBUG_FRAGMENTDB
# define DEBUG(a) Log.info() << a << std::endl;
#else
# define DEBUG(a)
#endif

using namespace std;

namespace BALL 
{

	FragmentDB::NoFragmentNode::NoFragmentNode(const char* file, int line, const string& filename)
		: Exception::GeneralException(file, line, "NoFragmentNode", 
											 string("the resource database does not contain a valid Fragment entry: ") + filename),
			filename_(filename)
	{
	}

	// default constructor
	FragmentDB::FragmentDB()
		:	valid_(false)
	{
		normalize_names.setFragmentDB(*this);
		add_hydrogens.setFragmentDB(*this);
		build_bonds.setFragmentDB(*this);
	}


	FragmentDB::FragmentDB(const String& filename)
		:	valid_(false)
	{
		if (filename == "") {
			Log.error() << "FragmentDB: using an empty filename is deprecated. Simulating previous behaviour by only using ResourceFile storage." << std::endl;
			Log.error() << "FragmentDB: Note that this may change in a near future release." << std::endl;
			/* RFFS's default constructor simulates the old behaviour, "fragments/Fragments.db" */
			stores_.insert(boost::shared_ptr<FragmentStorage>(new ResourceFileFragmentStorage()));
		}
		else
		{
			/* construct it with a custom filename. */
			stores_.insert(boost::shared_ptr<FragmentStorage>(new ResourceFileFragmentStorage(filename)));
		}
		normalize_names.setFragmentDB(*this);
		add_hydrogens.setFragmentDB(*this);
		build_bonds.setFragmentDB(*this);
	}

	void FragmentDB::addStore(boost::shared_ptr<FragmentStorage>& store)
	{
		stores_.insert(store);
	}

	FragmentDB::FragmentDB(const FragmentDB& db, bool /* deep */)
		:	valid_(false)
	{
		valid_ = db.isValid();
		stores_.insert(db.stores_.begin(),db.stores_.end());
	}

	FragmentDB::~FragmentDB()
	{
		/* the destruction of stores_ and shared_ptr take care of everything... */
	}

	bool FragmentDB::isValid() const 
	{
		return valid_;
	}

	bool FragmentDB::has(const String& fragment_name) const 
	{
		if (!isValid()) return false;

		NameFragmentQuery q(fragment_name);
		return query(q);
	}

	bool FragmentDB::query(FragmentQuery &the_query) const
	{
		FragmentStoreSet::iterator store;
		bool found = false;
		for (store = stores_.begin(); store != stores_.end(); ++store)
		{
			found |= (*store)->query(the_query);
		}
		return found;
	}

	FragmentDB::Type BALL_DEPRECATED FragmentDB::getFragmentType(const String& fragment_name) const
	{
		// did anyone ever use this? The original backend does not have that data...
		// I'll deprecate it for now. -- wolfgang, 2011
		NameFragmentQuery theFragment(fragment_name);
		if (!query(theFragment))
		{
			return FragmentDB::TYPE__UNKNOWN;
		}
		if (theFragment.getResults().size() > 1)
		{
			Log.info() << "FragmentDB: More than one result for query \"" << fragment_name << "\". Returning first match!" << std::endl;
		}
		boost::shared_ptr<Residue> fragment(*(theFragment.getResults().begin()));
		if (fragment->hasProperty("TYPE"))
		{
			return fragment->getProperty("TYPE").getInt();
		}
		else
		{
			return FragmentDB::TYPE__UNKNOWN;
		}
	}

	FragmentDB& FragmentDB::operator = (const FragmentDB& db)
	{
		/* save the other set of pointers.
		   in case db == *this, we'd throw away the store pointers otherwise. */
		FragmentStoreSet otherSet = db.stores_;
		stores_.clear();
		stores_.insert(otherSet.begin(),otherSet.end());
		return *this;
	}

	static String fragmentdb_pdb_dummy_ = "PDB";

	const String& FragmentDB::getDefaultNamingStandard() const 
	{
		/* FIXME: compute this from the stores */
		return fragmentdb_pdb_dummy_;
	}


	const Fragment* FragmentDB::getFragment(const String& fragment_name) const
	{
		Fragment* p_copy = getFragmentCopy(fragment_name);
		if (p_copy != NULL)
		{
			Log.warn() << "FragmentDB::getFragment() returning a copy Fragment* where non-copy requested. This _WILL_ leak memory! Use getFragmentCopy() or query() instead." << std::endl;
		}
		return p_copy;
	}

	Fragment* FragmentDB::getFragmentCopy(const String& fragment_name) const
	{
		/* this is OK, since Residue is-a Fragment */
		return getResidueCopy(fragment_name);
	}

	Residue* FragmentDB::getResidueCopy(const String& fragment_name) const
	{
		NameFragmentQuery q(fragment_name);
		if (!query(q))
		{
			return NULL;
		}
		if (q.getResults().size() > 1)
		{
			Log.info() << "FragmentDB: More than one result for query \"" << fragment_name << "\". Returning first match!" << std::endl;
			Log.info() << "FragmentDB: You should consider using query() instead of getFragment() to obtain all matching fragments." << std::endl;
		}
		return new Residue(*q.getResults().begin()->get());
	}

	Molecule* FragmentDB::getMoleculeCopy(const String& fragment_name) const
	{
		Fragment* ref_fragment = getFragmentCopy(fragment_name);
		Molecule* copy = 0;

		// copy the reference fragment if we found a reference fragment
		// (pointer != 0) and insert it into a new molecule.
		// Otherwise, return the NULL pointer.
		if (ref_fragment !=	0)
		{
			copy = new Molecule;
			copy->insert(*ref_fragment);
		}

		return copy;
	}

	const Fragment* FragmentDB::getReferenceFragment(const Fragment& fragment) const
	{
		boost::shared_ptr<Residue> reference = findReferenceFragment(fragment);
		if (reference)
		{
			Log.warn() << "FragmentDB: returning a copy Fragment* where non-copy requested. This _WILL_ leak memory! Use findReferenceFragment() or query() instead." << std::endl;
			return new Fragment(*reference.get());
		}
		else
		{
			return NULL;
		}
	}

	boost::shared_ptr<Residue> FragmentDB::findReferenceFragment(const Fragment &fragment) const
	{
	
		/* FIXME: this should be propagated downward into a specific FragmentQuery that can
			make use of more properties of &fragment than is used here. (we only use the name here).
		*/
	
		NameFragmentQuery q(fragment.getName(),"PDB",0);
		if (!query(q))
		{
			return boost::shared_ptr<Residue>();
		}

		if (q.getResults().size() == 1)
		{
			return *q.getResults().begin();
		}

		// now find the variant that best matches the fragment
		// This returns N/C terminal variants for fragments
		// that have the corresponding properties set or 
		// cystein variants without thiol hydrogen if the
		// disulphide bond property is set
		
		// First, check for two special properties of amino acids:
		// C_TERMINAL and N_TERMINAL 
		// They are usually not set, so set them here
		// As the fragment should be const, we store the properties
		// in a bit vector and OR them later with the fragment's properties
		BitVector	additional_properties;
		const Residue* residue = dynamic_cast<const Residue*>(&fragment);
		if (residue != 0)
		{
			if (residue->isCTerminal())
			{
				additional_properties.setBit(Residue::PROPERTY__C_TERMINAL);
			}
			if (residue->isNTerminal())
			{
				additional_properties.setBit(Residue::PROPERTY__N_TERMINAL);
			}
		}
		else
		{
			const Nucleotide* nucleotide = dynamic_cast<const Nucleotide*>(&fragment);
			if (nucleotide != 0)
			{
				if (nucleotide->is3Prime())
				{
					additional_properties.setBit(Nucleotide::PROPERTY__3_PRIME);
				}
				if (nucleotide->is5Prime())
				{
					additional_properties.setBit(Nucleotide::PROPERTY__5_PRIME);
				}
			}
			else
			{
				DEBUG(" neither residue nor nucleotide!")
			}
		}
	
		boost::shared_ptr<Residue> variant;

		// the number of properties that matched.
		// the fragment with the largest number of matched
		// properties is returned
		Index number_of_properties = -1;
		Index property_difference = -1;
		Index best_number_of_properties = -1;
		Index best_property_difference = 10000;

		// Iterate over all variants of the fragment and compare the properties.
		FragmentQuery::ResultSet::iterator it = q.getResults().begin();
		for (; it != q.getResults().end(); ++it)
		{
			// determine how many properties both have in common
			// by ANDing both bitvectors and counting ones
			BitVector props = fragment.getBitVector();
			props |= additional_properties;
			property_difference = (int)abs((int)props.countValue(true) - (int)(*it)->getBitVector().countValue(true));
			DEBUG(" props = " << props << "  bv = " << (*it)->getBitVector() << "   add = " << additional_properties)

			props &= (*it)->getBitVector();
			number_of_properties = (int)props.countValue(true);
			DEBUG(" considering variant " << (*it)->getName() << ". # properties: " << number_of_properties)

			if ((number_of_properties > best_number_of_properties)
					|| ((number_of_properties == best_number_of_properties) 
							&& (property_difference < best_property_difference)))
			{
				variant = *it;
				best_number_of_properties = number_of_properties;
				best_property_difference = property_difference;
			}
		}

		return variant;
	}

	const Residue* FragmentDB::getResidue(const String& fragment_name) const 
	{
		Residue* p_copy = getResidueCopy(fragment_name);
		if (p_copy != NULL)
		{
			Log.warn() << "FragmentDB::getResidue(): returning a copy Residue* where non-copy requested. This _WILL_ leak memory! Use getFragmentCopy() or query() instead." << std::endl;
		}
		return p_copy;
	}

	list<String> FragmentDB::getVariantNames(const String& name) const
	{
		NameFragmentQuery q(name, "PDB", /* unlimited */ 0);
		list<String> names;

		if (!query(q))
		{
			return names;
		}

		FragmentQuery::ResultSet::iterator result;
		for (result = q.getResults().begin(); result != q.getResults().end(); ++result)
		{
			names.push_back((*result)->getName());
		}

		return names;
	}


		const std::vector<Residue*>& FragmentDB::getFragments() const
		{
			// TODO
		}


} // namespace BALL
