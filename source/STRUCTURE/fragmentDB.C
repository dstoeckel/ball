// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//

#include <BALL/STRUCTURE/fragmentDB.h>

#include <BALL/KERNEL/PTE.h>
#include <BALL/KERNEL/nucleotide.h>
#include <BALL/KERNEL/nucleicAcid.h>
#include <BALL/KERNEL/chain.h>
#include <BALL/SYSTEM/path.h>
#include <BALL/KERNEL/forEach.h>
#include <BALL/MATHS/matrix44.h>
#include <BALL/FORMAT/resourceFile.h>

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

	FragmentDB::Type FragmentDB::getFragmentType(const String& fragment_name) const 
	{
/*
		FIXME: This probably needs to go in an extra FragmentQuery::QuerySelector
		*******
		if (!isValid() || 
				!tree->isValid() ||
				!has(fragment_name))
		{
			return FragmentDB::TYPE__UNKNOWN;
		}

		String path = (*name_to_path_.find(fragment_name)).second;
		path += "/Type";
		ResourceEntry* entry = tree->findChild(path);
		entry = tree->findChild("");
			
		if (entry!= 0)
		{
			if (entry->getValue() == "residue")
			{
				return FragmentDB::TYPE__RESIDUE;
			}

			if (entry->getValue() == "molecule")
			{
				return FragmentDB::TYPE__MOLECULE;
			}

			if (entry->getValue() == "fragment")
			{
				return FragmentDB::TYPE__MOLECULE;
			}
		}
*/
		return FragmentDB::TYPE__UNKNOWN;
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
 
	/////////////////////////////////////////////////////////////////////
	// NormalizeNamesProcessor	
	/////////////////////////////////////////////////////////////////////
	FragmentDB::NormalizeNamesProcessor::NormalizeNamesProcessor() 
		: UnaryProcessor<Fragment>() 
	{
		fragment_db_ = 0;
		naming_standard_ = "";
	}

	FragmentDB::NormalizeNamesProcessor::NormalizeNamesProcessor(FragmentDB& db)
	{
		setFragmentDB(db);
	}

	FragmentDB::NormalizeNamesProcessor::~NormalizeNamesProcessor()
	{
	}

	void FragmentDB::NormalizeNamesProcessor::setFragmentDB(FragmentDB& db)
	{
		fragment_db_ = &db;
		setNamingStandard(db.getDefaultNamingStandard());
	}

	void FragmentDB::NormalizeNamesProcessor::setNamingStandard(const String& naming_standard) 
	{
		naming_standard_ = naming_standard;
	}

	const String& FragmentDB::NormalizeNamesProcessor::getNamingStandard() 
	{
		return naming_standard_;
	}

	StringHashMap<StringHashMap<String> >& FragmentDB::getNamingStandards() 
	{
		throw std::exception();
		// FIXME return standards_;
	}

	const StringHashMap<String>& FragmentDB::getNamingStandard(const String& standard) const
	{
		throw std::exception();
		// FIXME return standards_[std];
	}

	std::vector<String> FragmentDB::getAvailableNamingStandards() const
	{
		throw std::exception();
		// FIXME return result;
	}

	Processor::Result FragmentDB::NormalizeNamesProcessor::operator () (Fragment& fragment) 
	{
		fragments_.push_back(&fragment);

		return Processor::CONTINUE;
	}

	bool FragmentDB::NormalizeNamesProcessor::start() 
	{
		fragments_.clear();
		return true;
	}


	// match an RES/ATOM pair in a map
	bool FragmentDB::NormalizeNamesProcessor::matchName(String& res_name, String&	atom_name, const FragmentDB::NameMap&	map) const
	{
		// residue name (non const)
		res_name.trim();
		String	s[2];

		NameMap::ConstIterator it;
		it = map.find(res_name + ":*");
		if (it != map.end())
		{
			it->second.split(s, 2, ":");
			res_name = s[0];
		}
		
		// atom name (non const)
		atom_name.trim();

		bool hit = false;

		// first, try to match exactly
		it = map.find(res_name + ":" + atom_name);
		if (it != map.end())
		{
			it->second.split(s, 2, ":");
			atom_name = s[1];
			res_name = s[0];
			hit = true;
		} 
		else 
		{
			// second, try wildcard match for residue names
			it = map.find("*:" + atom_name);
			if (it != map.end())
			{
				it->second.split(s, 2, ":");
				atom_name = s[1];
				hit = true;
			}
		}

		return hit;
	}

	bool FragmentDB::NormalizeNamesProcessor::finish() 
	{
		/*
			The _good_ way to do it would be to operate on a by-fragment basis, maybe similar to the current code,
			using the non-chained fragments to determine the map first, by getting "all translations of the fragment",
			graph isomorphism wrt. names, and choosing the standard based on that. After that, just get the Fragment
			readily translated from the backend with NameFragmentQuery.
			
			For now, I'll just work my necromancy and conjure the translation tables back up from the backend.
			-- wolfgang, Jul 2011
		*/
		if (fragment_db_ == 0)
		{
			return false;
		}

		const char* error_msg = "FragmentDB: cannot locate an appropriate name conversion table!";
		const String map_name = "-" + naming_standard_;

		NameMapQuery query(map_name);
		if (!fragment_db_->query(query)) {
			Log.error() << error_msg << endl;
			return false;
		}
		const StringHashMap<NameMapQuery::NameMap*>& table = query.getMaps();

		HashMap<NameMapQuery::NameMap*, Index> usable_maps;

		Log.info() << "Got Name Maps:" << std::endl;
		for (StringHashMap<NameMapQuery::NameMap*>::ConstIterator it = table.begin(); it != table.end(); ++it)
		{
			usable_maps[it->second] = 0;
			Log.info() << it->first << std::endl;
		}


		//We now sort the fragments into parent containers if available. The rational
		//is that we should get a more stable estimate of the applicable naming scheme
		//if we iterate over a set of fragments than applying naming schemes to a
		//single Fragment, as there might by errors which could lead to the selection
		//of a wrong naming scheme.
		map<AtomContainer*, list<Fragment*> > parent_containers;

		list<Fragment*>::iterator frag_it = fragments_.begin();
		while(frag_it != fragments_.end()) {
			Residue* residue = 0;
			Nucleotide* nacid = 0;

			if((residue = RTTI::castTo<Residue>(**frag_it)) && residue->getChain())
			{
				parent_containers[static_cast<AtomContainer*>(residue->getChain())].push_back(*frag_it);
				frag_it = fragments_.erase(frag_it);
			}
			else if ((nacid = RTTI::castTo<Nucleotide>(**frag_it)) && nacid->getNucleicAcid())
			{
				parent_containers[static_cast<AtomContainer*>(nacid->getNucleicAcid())].push_back(*frag_it);
				frag_it = fragments_.erase(frag_it);
			}
			else
			{
				++frag_it;
			}
		}

		const NameMap* map;

		//First deal with the remaining fragments
		for(list<Fragment*>::iterator frag_it = fragments_.begin(); frag_it != fragments_.end(); ++frag_it) {
			countHits_(usable_maps, *frag_it, OVERWRITE);
			// we found an appropriate map, so use it
			if ((map = getBestMap_(usable_maps)) != 0)
			{
				normalizeFragment_(map, *frag_it);
			}
			else
			{
				// if we couldn't find an appropriate table, complain about it!
				Log.error() << error_msg << endl;
			}
		}

		//Now look at all fragments that could be assigned to a chain
		std::map<AtomContainer*, list<Fragment*> >::iterator chain_it;
		for (chain_it = parent_containers.begin(); chain_it != parent_containers.end(); ++chain_it)
		{
			countHits_(usable_maps, chain_it->second);
			// we found an appropriate map, so use it
			if ((map = getBestMap_(usable_maps)) != 0)
			{
				normalizeFragments_(map, chain_it->second);
			}
			else
			{
				// if we couldn't find an appropriate table, complain about it!
				Log.error() << error_msg << endl;
			}

		}

		return true;
	}

	String FragmentDB::NormalizeNamesProcessor::getSuffix_(const Fragment* frag) const
	{
		// determine whether the fragment is an amino acid
		// if it is: determine the correct name for N-,C-terminal AA
		const Residue* const residue = RTTI::castTo<Residue>(*frag);
		if (residue != 0)
		{
			if (residue->isCTerminal()) return String("-C");
			if (residue->isNTerminal()) return String("-N");
		}

		return String();
	}

	bool FragmentDB::NormalizeNamesProcessor::doMatch_(String& res_name, const String& res_name_suffix, String& atom_name, const NameMap& name_map) const
	{
		// first, try to match exactly
		String match_name = res_name + res_name_suffix;
		if (res_name_suffix.isEmpty() || !matchName(match_name, atom_name, name_map))
		{
			// try to match non-terminal residues
			if (!matchName(res_name, atom_name, name_map))
			{
				match_name = "*" + res_name_suffix;
				if (res_name_suffix.isEmpty() || !matchName(match_name, atom_name, name_map))
				{
					match_name = "*";
					if (!matchName(match_name, atom_name, name_map))
					{
						return false;
					}
				}
			}
		}

		return true;
	}

	void FragmentDB::NormalizeNamesProcessor::countHits_(HashMap<NameMap*, Index>& maps, const std::list<Fragment*>& frags)
	{
		if(frags.size() == 0) {
			return;
		}

		list<Fragment*>::const_iterator it = frags.begin();
		countHits_(maps, *it, OVERWRITE);

		for (++it; it != frags.end(); ++it)
		{
			countHits_(maps, *it, ADD);
		}
	}

	void FragmentDB::NormalizeNamesProcessor::countHits_(HashMap<NameMap*, Index>& maps, const Fragment* frag, CountingMode mode)
	{
		String atom_name;
		String res_name = frag->getName();
		AtomConstIterator atom_it;
		HashMap<NameMap*, Index>::Iterator map_iterator;

		for (map_iterator = maps.begin(); map_iterator != maps.end(); ++map_iterator)
		{
			Size hit_counter = 0;
			String res_name_suffix = getSuffix_(frag);

			for (atom_it = frag->beginAtom(); +atom_it; ++atom_it)
			{
				atom_name = atom_it->getName();

				if(doMatch_(res_name, res_name_suffix, atom_name, *map_iterator->first)) {
					hit_counter++;
				}
			}

			// update hit_count for each map
			if(mode == ADD) {
				map_iterator->second += hit_counter;
			} else {
				map_iterator->second = hit_counter;
			}
		}
	}

	const FragmentDB::NameMap* FragmentDB::NormalizeNamesProcessor::getBestMap_(const HashMap<NameMap*, Index>& maps) const
	{
		// these two variables are needed to store the best map
		Index max_hits = 0;
		NameMap* result = 0;

		// look for the best map
		HashMap<NameMap*, Index>::ConstIterator map_iterator;
		for (map_iterator = maps.begin(); map_iterator != maps.end(); ++map_iterator)
		{
			if (map_iterator->second > max_hits)
			{
				max_hits = (*map_iterator).second;
				result   = (*map_iterator).first;
			}
		}

		return result;
	}

	void FragmentDB::NormalizeNamesProcessor::normalizeFragments_(const NameMap* name_map, const std::list<Fragment*>& frags)
	{
		for (list<Fragment*>::const_iterator frag_it = frags.begin(); frag_it != frags.end(); ++frag_it)
		{
			normalizeFragment_(name_map, *frag_it);
		}
	}

	void FragmentDB::NormalizeNamesProcessor::normalizeFragment_(const NameMap* name_map, Fragment* frag)
	{
		String atom_name;
		// extract the residue name
		String res_name = frag->getName();
		String res_name_suffix = getSuffix_(frag);

		// now, iterate over the fragment`s atoms
		for (AtomIterator atom_it = frag->beginAtom(); +atom_it; ++atom_it)
		{
			// get the atom name
			atom_name = atom_it->getName();

			if(doMatch_(res_name, res_name_suffix, atom_name, *name_map)) {
				atom_it->setName(atom_name);
				atom_it->getFragment()->setName(res_name);
			}
		}
	}
	/////////////////////////////////////////////////////////////////////
	//	BuildBondsProcessor
	/////////////////////////////////////////////////////////////////////

	FragmentDB::BuildBondsProcessor::BuildBondsProcessor()
		:	fragment_db_(0),
			bonds_built_(0)
	{
	}

	FragmentDB::BuildBondsProcessor::BuildBondsProcessor(const FragmentDB& db)
		: fragment_db_(const_cast<FragmentDB*>(&db)),
			bonds_built_(0)
	{
	}

	FragmentDB::BuildBondsProcessor::~BuildBondsProcessor()
	{
		fragment_db_ = 0;
	}

	void FragmentDB::BuildBondsProcessor::setFragmentDB(const FragmentDB& db)
	{
		fragment_db_ = &const_cast<FragmentDB&>(db);
	}

	bool FragmentDB::BuildBondsProcessor::start()
	{
		bonds_built_ = 0;
		return true;
	}

	Size FragmentDB::BuildBondsProcessor::getNumberOfBondsBuilt()
	{
		return bonds_built_;
	}

	bool FragmentDB::BuildBondsProcessor::finish()
	{
		bool ok = true;

		try
		{
			// if there are no inter-fragment bonds, return
			if (connections_.size() >= 2)
			{
				ConnectionList::iterator it1(connections_.begin());
				ConnectionList::iterator it2;
				for (; it1 != connections_.end(); ++it1)
				{
					for (it2 = it1, ++it2; it2 != connections_.end(); ++it2)
					{
						if ((it1->atom != 0) && (it2->atom != 0))
						{
							if (buildConnection_(*it1, *it2))
							{
								// Remember we built a bond
								bonds_built_++;
								// Remove the connection we made from the list of connections
								it1->atom = 0;	
								it2->atom = 0;
							}
						}
					}
				}
			}
		}
		catch(Exception::TooManyBonds&)
		{
			ok = false;
		}

		// Clear the connection list
		connections_.clear();
		return ok;
	}


	bool FragmentDB::BuildBondsProcessor::buildConnection_(Connection& con1, 
																												 Connection& con2)
	{
		if (con1.type_name != con2.connect_to || 
				con1.connect_to != con2.type_name)
		{
			return false;
		}

		// if the two connection types match,
		// check for distance condition and the two atoms
		const float distance = con1.atom->getPosition().getDistance(con2.atom->getPosition());
		if (fabs(con1.dist - distance) > con1.delta || 
			  fabs(con2.dist - distance) > con2.delta)
		{
			return false;
		}

		// create the bond only if it does not exist
		if (con1.atom->isBoundTo(*con2.atom)) return false;
		
		// create the bond
		Bond* const bond = con1.atom->createBond(*con2.atom);

		if (bond != 0)
		{
			bond->setOrder(con1.order);
			if (con1.order != con2.order)
			{
				Log.warn() << "FragmentDB::BuildBondsProcessor: inconsistent bond orders for connection between " 
									 << con1.atom->getFullName() << " and " << con2.atom->getFullName() << std::endl;
			}
			return true;
		}
		
		return false;
	}


	Processor::Result FragmentDB::BuildBondsProcessor::operator () (Fragment& fragment)
	{
		try
		{
			// build all bonds in the fragment
			bonds_built_ += buildFragmentBonds(fragment);
		}
		catch(Exception::TooManyBonds&)
		{
			return Processor::ABORT;
		}

		return Processor::CONTINUE;
	}


	Size FragmentDB::BuildBondsProcessor::buildFragmentBonds(Fragment& fragment, const Fragment& tplate)
	{
		// abort immediately if no fragment DB is known
		if (fragment_db_ == 0) return 0;

		DEBUG("FragmentDB::BuildBondsProcessor: building bonds for " 
							 << fragment.getName() << " from template " << tplate.getName())

		StringHashMap<const Atom*> template_names;

		for (AtomConstIterator catom_it = tplate.beginAtom(); +catom_it; ++catom_it)
		{
			const String atom_name = catom_it->getName().trim();
#ifdef BALL_DEBUG_FRAGMENTDB
			if (template_names.has(atom_name))
			{
				DEBUG("FragmentDB::BuildBondsProcessor: duplicate atom name in template " << tplate.getName())
			}
#endif
			template_names.insert(atom_name, &*catom_it);
		}
		
		// count the counds we build...
		Size bonds_built = 0;
		
		// iterate over all atoms in the tplate
		Atom::BondConstIterator	tplate_bond_it;

		// iterate over all fragment atoms 
		for (AtomIterator frag_atom_it = fragment.beginAtom(); +frag_atom_it; ++frag_atom_it) 
		{
			const String atom_name = frag_atom_it->getName().trim();
			const StringHashMap<const Atom*>::Iterator to_find = template_names.find(atom_name);
			if (to_find == template_names.end()) 
			{
				continue;	
			}
			
			const Atom* tplate_atom = to_find->second;
			
			// if tplate_atom has a connection, store it for &*frag_atom_it.
			// so inter-fragment connections can be built in the finish() step.
			if (tplate_atom->hasProperty("CONNECTION")) {
				DEBUG(tplate_atom->getName() << " has a Connection entry (inter-fragment bond).")
				boost::shared_ptr<Connection> con;
				con = boost::dynamic_pointer_cast<Connection>(tplate_atom->getProperty("CONNECTION").getSmartObject());
				if (con) {
					Connection conn = *con;
					conn.atom = &*frag_atom_it;
					DEBUG("Storing connection entry (" << conn.type_name << ") for later.")

					connections_.push_back(conn);
				}
			}

			// we found two matching atoms. Great! Now check their bonds..
			// iterate over all bonds of the template
			for (tplate_bond_it = tplate_atom->beginBond(); +tplate_bond_it; ++tplate_bond_it) 
			{
				const Atom* partner = tplate_bond_it->getPartner(*tplate_atom);
				// if we found the bond partner, create the new bond
				if (partner == 0) continue;
			
				const String name = partner->getName();

				// look in the fragment for the correct partner of the current bond in the template
				AtomIterator second_frag_it(fragment.beginAtom());
				for (; +second_frag_it; ++second_frag_it) 
				{
					if (second_frag_it->getName().trim() != name) continue;
				
					// ok, we found the correct partner atom
					// does the bond already exists?
					Bond* bond = second_frag_it->getBond(*frag_atom_it);
					if (bond == 0)
					{
						// no, so create it											
						bond = frag_atom_it->createBond(*second_frag_it);
					}

					// assign the correct bond order, name, and type
					// (even if the bond exists -- to correct PDB CONECT entries)
					if (bond != 0)
					{
						// assign the bond type and order
						bond->setOrder(tplate_bond_it->getOrder());
						bond->setType(tplate_bond_it->getType());
						bond->setName(tplate_bond_it->getName());

						// count this bond 
						bonds_built++;
					}
					break;
				}
			}
		}

		return bonds_built;
	}


	Size FragmentDB::BuildBondsProcessor::buildFragmentBonds(Fragment& fragment)
	{
		// abort immediately if no fragment DB is known
		if (fragment_db_ == 0) return 0;

		// check whether our DB knows the fragment and retrieve the template
		boost::shared_ptr<Residue> tplate = fragment_db_->findReferenceFragment(fragment);
		if (!tplate) return 0;
		
		return buildFragmentBonds(fragment, *tplate);
	}

	Size FragmentDB::BuildBondsProcessor::buildInterFragmentBonds(Fragment& first, Fragment& second) const
	{
		if (fragment_db_ == 0) return 0;

		// apply the whole processor to just the two arguments.
		// retained for compatibility reasons
		// this inner instance is needed for const-correctness, obviously.

		BuildBondsProcessor temp_pocessor(*fragment_db_);
		temp_pocessor.start();
		temp_pocessor(first);
		temp_pocessor(second);
		temp_pocessor.finish();

		return temp_pocessor.getNumberOfBondsBuilt();
	}

		const std::vector<Residue*>& FragmentDB::getFragments() const
		{
			// TODO
		}


} // namespace BALL
