#include <BALL/STRUCTURE/fragmentDB.h>
#include <BALL/KERNEL/nucleotide.h>
#include <BALL/KERNEL/nucleicAcid.h>
#include <BALL/KERNEL/chain.h>

//#define BALL_DEBUG_BUILDBONDSPROCESSOR

#ifdef BALL_DEBUG_BUILDBONDSPROCESSOR
# define DEBUG(a) Log.info() << a << std::endl;
#else
# define DEBUG(a)
#endif


namespace BALL {
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

}
