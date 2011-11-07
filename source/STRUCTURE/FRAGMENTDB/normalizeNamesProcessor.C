#include <BALL/STRUCTURE/fragmentDB.h>
#include <BALL/STRUCTURE/FRAGMENTDB/nameMapQuery.h>
#include <BALL/MATHS/matrix44.h>
#include <BALL/KERNEL/nucleotide.h>
#include <BALL/KERNEL/nucleicAcid.h>
#include <BALL/KERNEL/chain.h>
#include <BALL/KERNEL/atomContainer.h>
#include <map>

namespace BALL {

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

		const char* error_msg = "FragmentDB: no name conversion found, assuming PDB nomenclature.";
		const String map_name = "-" + naming_standard_;

		NameMapQuery query(map_name);
		if (!fragment_db_->query(query)) {
			Log.error() << error_msg << "(Tried to find: "<<map_name <<")" << endl;
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
		std::map<AtomContainer*, list<Fragment*> > parent_containers;
	
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
				StringHashMap<NameMapQuery::NameMap*>::ConstIterator it = table.begin();
				for (; it != table.end(); ++it) {
					Log.error() << "Score:" << it->first << " = " << usable_maps[it->second] << std::endl;
				}
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

}
