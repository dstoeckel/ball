#include <BALL/STRUCTURE/FRAGMENTDB/resourceFileFragmentStorage.h>
#include <BALL/STRUCTURE/FRAGMENTDB/fragmentQuery.h>
#include <BALL/STRUCTURE/FRAGMENTDB/nameMapQuery.h>
#include <BALL/STRUCTURE/FRAGMENTDB/propertyFragmentQuery.h>
#include <BALL/STRUCTURE/fragmentDB.h>

#include <BALL/KERNEL/atom.h>
#include <BALL/KERNEL/bond.h>
#include <BALL/KERNEL/PTE.h>
#include <BALL/KERNEL/residue.h>
#include <BALL/KERNEL/nucleotide.h>
#include <BALL/DATATYPE/stringHashMap.h>
#include <BALL/SYSTEM/path.h>
#include <BALL/COMMON/limits.h>

#include <boost/smart_ptr/shared_ptr.hpp>

namespace BALL
{

	static const char* FRAGMENT_DB_INCLUDE_TAG = "#include";

	ResourceFileFragmentStorage::ResourceFileFragmentStorage(const String& path_to_root_resourcefile)
		: tree_(0),
		  filename_(path_to_root_resourcefile)
	{
		init();
	}

	ResourceFileFragmentStorage::~ResourceFileFragmentStorage() {
		delete tree_;
		for (std::vector<Residue*>::iterator it = fragments_.begin();
				 it != fragments_.end(); ++it)
		{
			delete (*it);
		}
	}

	String ResourceFileFragmentStorage::getVersionTag()
	{
		return tree_->getRoot().findEntry("Version")->getValue();
	}
	
	bool ResourceFileFragmentStorage::query(FragmentQuery &query)
	{
		bool found;
		if (query.selectsOn(FragmentQuery::QueryFragmentName)) {
			NameFragmentQuery* q;
			try {
				q = boost::any_cast<NameFragmentQuery*>(query.getSelectorDetail(FragmentQuery::QueryFragmentName));
			} catch (boost::bad_any_cast& error) {
				Log.error() << "ResourceFileFragmentStorage: Unable to process query: " << query.toString() << std::endl;
				return false;
			}

			// TODO: Naming Standard!
			// TODO: Delayed-Load is entirely possible here.
			
			if (query.getMaxResults() == 1)
			{
				if (name_to_frag_index_.has(q->getFragmentName())) {
					boost::shared_ptr<Residue> frag_copy_p(
						new Residue(*fragments_[name_to_frag_index_[q->getFragmentName()]])
					);
					query.addResult(frag_copy_p);
					found = true;
				}
			}
			else
			{
				// want several results, get all variants to that name.
				// as per FragmentQuery docs, max_results_ == 0 means unlimited.
				unsigned int remaining_results = query.getMaxResults();
				bool fetch_all(remaining_results == 0);
				
				StringHashMap<list<Position> >::ConstIterator to_find = name_to_variants_.find(q->getFragmentName());
				if (to_find == name_to_variants_.end())
				{	
					// TODO: maybe we want to try exact matching here then...
					return false;	
				}

				list<Position>::const_iterator it = (*to_find).second.begin();
				const list<Position>::const_iterator end_it = (*to_find).second.end();

				while (fetch_all || (remaining_results > 0))
				{
					if (it == end_it)
					{
					 break;
					}
					boost::shared_ptr<Residue> frag_copy_p(
						new Residue(*fragments_[*it])
					);
					query.addResult(frag_copy_p);
					found = true;
					++it;
					--remaining_results;
				}
			}
		}
		else if (query.selectsOn(FragmentQuery::QueryNameMap))
		{
			NameMapQuery* q;
			try {
				q = boost::any_cast<NameMapQuery*>(query.getSelectorDetail(FragmentQuery::QueryNameMap));
			} catch (boost::bad_any_cast& error) {
				Log.error() << "ResourceFileFragmentStorage: Unable to process query: " << query.toString() << std::endl;
				return false;
			}
			for (StringHashMap<NameMap>::Iterator it = standards_.begin(); it != standards_.end(); ++it)
			{
				if (it->first.hasSubstring(q->getMapName()))
				{
					found = true;
					q->addMap(it->first,& it->second);
				}
			}
		}
		else if (query.selectsOn(FragmentQuery::QueryFragmentProperties))
		{
			PropertyFragmentQuery* q;
			try {
				q = boost::any_cast<PropertyFragmentQuery*>(query.getSelectorDetail(FragmentQuery::QueryFragmentProperties));
			} catch (boost::bad_any_cast& error) {
				Log.error() << "ResourceFileFragmentStorage: Unable to process query: " << query.toString() << std::endl;
				return false;
			}
			BitVector tplate = q->getPropertyManager().getBitVector();

			// TODO: yes, this is as inefficient as it looks.
			std::vector<Residue*>::iterator fragment = fragments_.begin();
			for (; fragment != fragments_.end(); ++fragment)
			{
				if ((tplate & (*fragment)->getBitVector()) == tplate)
				{
					// all bits in the template are also present in the fragment
					// now check the named properties. (TODO: test equality here?)
					bool namedPropsOK = true;
					NamedPropertyConstIterator namedProp = q->getPropertyManager().beginNamedProperty();
					NamedPropertyConstIterator end = q->getPropertyManager().endNamedProperty();
					for (; (namedProp != end) && namedPropsOK; ++namedProp)
					{
						if (!(*fragment)->hasProperty(namedProp->getName()))
						{
							namedPropsOK = false;
						}
					}
					if (namedPropsOK)
					{
						boost::shared_ptr<Residue> frag_copy_p(
							new Residue(**fragment)
						);
						query.addResult(frag_copy_p);
					}
				}
			}
		}

		return found;
	}
	
	// ------------------------- from ex FragmentDB.C --------------------
	
	void ResourceFileFragmentStorage::expandTree_(ResourceEntry& root_entry)
	{
		bool expanded_one = true;
		while (expanded_one)
		{
			expanded_one = false;
			ResourceEntry::Iterator entry_it;
			for (entry_it = ++root_entry.begin(); +entry_it && !expanded_one; ++entry_it)
			{
				if (entry_it->getKey().hasPrefix(FRAGMENT_DB_INCLUDE_TAG))
				{
					expandFirst_(*entry_it);
					expanded_one = true;
					break;
				}
			}	
		}
	}

	bool ResourceFileFragmentStorage::expandFirst_(ResourceEntry& root_entry)
	{
		String key = root_entry.getKey();
		vector<String> key_fields;

		if (key.split(key_fields, ":") != 2)
		{
			// if the include directive is invalid, remove the entry
			Log.error() << "FragmentDB: illegal #include directive: " << key << std::endl;
			root_entry.getParent()->removeChild(key, 0);
			return false;
		} 
		else 
		{
			String value_fields[2];
			String value = root_entry.getValue();
			value.split(value_fields, 2, ":");
				
			ResourceEntry*	parent = root_entry.getParent();
			parent->removeChild(key, 0);

			// search in the standard fragment DB file
			Path path;
			String filename = path.find(value_fields[0]);
			if (filename == "")
			{
				throw Exception::FileNotFound(__FILE__, __LINE__, value_fields[0]);
			}

			ResourceFile file(filename);
			if (!file.isValid())
			{
				Log.error() << "FragmentDB: cannot open include file " << value_fields[0] << std::endl;
				return false;
			}
				
			ResourceEntry* tree_entry = file.getRoot().getEntry(value_fields[1]);
			if (tree_entry == 0)
			{
				Log.error() << "FragmentDB: cannot find node " << value_fields[1] << " in file " << value_fields[0] << std::endl;
			} 
			else 
			{
				ResourceEntry* entry = parent->insertChild(key_fields[1], tree_entry->getValue());
				entry->mergeChildrenOf(*tree_entry);
			}
		}

		return true;
	}

	Position ResourceFileFragmentStorage::addNewFragment_(Residue* fragment)
	{
		fragments_.push_back(fragment);
		return fragments_.size() - 1;
	}

	void ResourceFileFragmentStorage::parseAtoms_(ResourceEntry& entry, Fragment& fragment)
	{
		ResourceEntry::Iterator	entry_it;

		for (entry_it = ++entry.begin(); +entry_it; ++entry_it)
		{
			if (entry_it->getDepth() == entry.getDepth() + 1)
			{
				if (entry_it->getValue().countFields(" ") != 4)
				{
					Log.error() << "FragmentDB: wrong entry for atom " << entry_it->getKey() 
							 << ": " << entry_it->getValue() << std::endl;
				} 
				else 
				{
					// create a new atom...
					Atom*	atom = new Atom;
								
					// ...set its name and insert it into the fragment.
					atom->setName(entry_it->getKey());
					fragment.insert(*atom);
		
					// Now extract element and position (x, y, z) from the string.
					String		s[4];
					entry_it->getValue().split(s, BALL_SIZEOF_ARRAY(s), " ");
					Vector3	r(s[1].toFloat(), s[2].toFloat(), s[3].toFloat());
		
					// and assign its values to the atom
					atom->setPosition(r);
					atom->setElement(PTE.getElement(s[0]));
				}
			}
		}
	}

	void ResourceFileFragmentStorage::parseBonds_(ResourceEntry& entry, Fragment& fragment)
	{

		ResourceEntry::Iterator	entry_it;

		for (entry_it = ++entry.begin(); +entry_it; ++entry_it)
		{
			if (entry_it->getDepth() == entry.getDepth() + 1)
			{
				// check whether the fragment contains both bonds
				Atom*	atom1 = 0;
				Atom*	atom2 = 0;
				AtomIterator	atom_it;
				
				// first field contains a serial number, field 2 the first atom, field 3 the second atom
				// and the third field (optional) the bond type (s[ingle], d[ouble], a[romatic])
				String fields[3];
				entry_it->getValue().split(fields, 3);
		
				for (atom_it = fragment.beginAtom(); +atom_it; ++atom_it)
				{
					if (atom_it->getName() == fields[0])
					{
						atom1 = &*atom_it;
					}
					if (atom_it->getName() == fields[1])
					{
						atom2 = &*atom_it;
					}
				}
		
				if ((atom1 == 0) || (atom2 == 0))
				{
					// if at least on of the atoms doesn`t exist: complain about it
					Log.error() << "FragmentDB: Bond to a non-existing atom: " 
											<< fields[0] << "-" << fields[1] 
											<< " (in " << entry_it->getPath() << ")" << std::endl;
				} 
				else	
				{
					// otherwise create the bond, if valences free
					if ((atom1->countBonds() > Atom::MAX_NUMBER_OF_BONDS) || (atom2->countBonds() > Atom::MAX_NUMBER_OF_BONDS))
					{
						Log.error() << "FragmentDB: too many bonds - cannot create bond: " 
												<< atom1->getName() << "-" << atom2->getName()
												<< " in fragment " << fragment.getName() 
												<< " (in " << entry_it->getPath() << ")" << std::endl;
					} 
					else 
					{
						// create the bond
						Bond* bond = atom1->createBond(*atom2);

						if (bond != 0)
						{
							// by default, we create single bonds
							bond->setOrder(Bond::ORDER__SINGLE);

							// if the bond order is specified, set it
							// s == single, a == aromatic, d = double, t = triple
							if (fields[2] != "")
							{
								switch (fields[2][0])
								{
									case 'a':
										bond->setOrder(Bond::ORDER__AROMATIC); break;
									case 'd':
										bond->setOrder(Bond::ORDER__DOUBLE); break;
									case 't':
										bond->setOrder(Bond::ORDER__TRIPLE); break;
									case 's':
										bond->setOrder(Bond::ORDER__SINGLE); break;
									default:
										Log.error() << "FragmentDB::parseBonds_: unknown bond type " 
																<< fields[2] << " (in " << entry_it->getPath() << ")" << std::endl;
								}
							}
						}
					}
				}
			}
		}
	}

	void ResourceFileFragmentStorage::parseDelete_(ResourceEntry& entry, Fragment& fragment)
	{

		ResourceEntry::Iterator	entry_it;

		for (entry_it = ++entry.begin(); +entry_it; ++entry_it)
		{
			if (entry_it->getDepth() == entry.getDepth() + 1)
			{
				// check whether the fragment contains both bonds
				AtomIterator	atom_it;
				Atom*				atom = 0;
		
				for (atom_it = fragment.beginAtom(); +atom_it; ++atom_it)
				{
					if (atom_it->getName() == entry_it->getKey())
					{
						atom = &*atom_it;
					}
				}
		
				if (atom == 0)
				{
					// if the atom to be deleted doesn`t exist - complain about it!
					Log.error() << "FragmentDB: cannot delete non-existing atom: "
																			<< entry_it->getKey() << std::endl;
				} 
				else 
				{
					// otherwise delete the atom
					fragment.remove(*atom);
					delete atom;
				}
			}
		}
	}

	void ResourceFileFragmentStorage::parseRename_(ResourceEntry& entry, Fragment& fragment)
	{

		ResourceEntry::Iterator	entry_it;

		for (entry_it = ++entry.begin(); +entry_it; ++entry_it)
		{
			if (entry_it->getDepth() == entry.getDepth() + 1)
			{
				// check whether the fragment contains both bonds
				AtomIterator	atom_it;
				Atom*				atom = 0;
		
				for (atom_it = fragment.beginAtom(); +atom_it; ++atom_it)
				{
					if (atom_it->getName() == entry_it->getKey())
					{
						atom = &*atom_it;
					}
				}
		
				if (atom == 0)
				{
					// if the atom to be renamed doesn`t exist - complain about it!
					Log.error() << "FragmentDB: cannot rename non-existing atom: "
																			<< entry_it->getKey() << std::endl;
				} 
				else 
				{
					// otherwise rename the atom
					atom->setName(entry_it->getValue());
				}
			}
		}
	}
	
	void ResourceFileFragmentStorage::parseConnections_(ResourceEntry&  entry, Fragment& fragment)
	{
		const ResourceEntry* first_entry = &entry;

		if (first_entry == 0) return;

		ResourceEntry::ConstIterator	it1 = first_entry->begin();
		for (++it1; +it1; ++it1)
		{
			// split the fields of the "Connections" entry.
			// It should have the following format:
			//   (<name> <atom_name> <match_name> <distance> <tolerance>)
			//	<name>:				Name of the connection type (eg C-term)
			//	<atom_name>:	Name of the atom that might create the connection
			//  <bond_order>: s/d/t/a (single/double/triple/aromatic)
			//	<match_name>:	Name of a matching connection type: this connection is 
			//								created if the two names match
			//	<distance>:		Distance of the connection in Angstrom
			//	<tolerance>:	Tolerance: connection will be built only if the distance
			//								of the two atoms within <tolerance> of <distance>
			//	Example entry:
			//		(C-term C s N-term 1.33 0.5):
			//			This will build a connection to a fragment with a N-term connection
			//			if the two atoms are 1.33+/-0.5 Angstrom apart. The bond is a single bond.
			
			String	s[6];
			it1->getValue().split(s, 6);
			Connection conn;
			conn.atom = 0;
			for (AtomIterator ai = fragment.beginAtom(); +ai; ++ai)
			{
				if (ai->getName() == s[0])
				{
					conn.atom = &*ai;
					break;
				}
			}
			// If there is a matching atom, store the connection.
			if (conn.atom != 0)
			{
				conn.type_name = it1->getKey();
				conn.connect_to = s[1];
				
				conn.dist = s[3].toFloat();
				conn.delta = s[4].toFloat();
				// set the bond order
				switch (s[2][0])
				{
					case 's': conn.order = Bond::ORDER__SINGLE; break;
					case 'd': conn.order = Bond::ORDER__DOUBLE; break;
					case 't': conn.order = Bond::ORDER__TRIPLE; break;
					case 'a': conn.order = Bond::ORDER__AROMATIC; break;
					default:
						Log.warn() << "FragmentDB::BuildBondsProcessor: unknown bond order " 
											 << s[2] << " (in " << first_entry->getPath() << ")" << std::endl;
				}
				
				boost::shared_ptr<PersistentObject> connptr(new Connection(conn));
				conn.atom->setProperty(NamedProperty("CONNECTION", connptr));
			}
		}

	}

	void ResourceFileFragmentStorage::parseProperties_(ResourceEntry&  entry, PropertyManager& property_man)
	{
		ResourceEntry::Iterator	entry_it;
		String property;
		bool invert = false;

		for (entry_it = ++entry.begin(); +entry_it; ++entry_it)
		{
			if (entry_it->getDepth() == entry.getDepth() + 1)
			{
				// Check for the most important properties: all those defined
				// in Residue::PROPERTIES
				property = entry_it->getKey();
				property.toUpper();
				if (property[0] == '!')
				{
					property.erase(0, 1);
					invert = true;
				} 
				else 
				{
					invert = false;
				}
				

				Property prop = Limits<Property>::max();
				if (property == "NON_STANDARD")
				{
					prop = Residue::PROPERTY__NON_STANDARD;
				}
				else if (property == "AMINO_ACID")
				{
					prop = Residue::PROPERTY__AMINO_ACID;
				}
				else if (property == "WATER")
				{
					prop = Residue::PROPERTY__WATER;
				}
				else if (property == "HAS_SSBOND")
				{
					prop = Residue::PROPERTY__HAS_SSBOND;
				}
				else if (property == "C_TERMINAL")
				{
					prop = Residue::PROPERTY__C_TERMINAL;
				}
				else if (property == "N_TERMINAL")
				{
					prop = Residue::PROPERTY__N_TERMINAL;
				}
				else if (property == "NUCLEOTIDE")
				{
					prop = Nucleotide::PROPERTY__NUCLEOTIDE;
				}

				if (prop == Limits<Property>::max())
				{
					// if the property was not recognized,
					// store it as a name-value pair
					if (invert)
					{
						property_man.clearProperty(property.c_str());
					} 
					else	
					{
						property_man.setProperty(property.c_str());
					}
				} 
				else 
				{
					if (invert)
					{
						property_man.clearProperty(prop);
					} 
					else 
					{
						property_man.setProperty(prop);
					}
				}
			}
		}
		
	}

	void ResourceFileFragmentStorage::parseType_(ResourceEntry&  entry, Fragment& fragment)
	{
		ResourceEntry* type_entry = entry.findChild("");

		if (type_entry!= 0)
		{
			if (type_entry->getValue() == "residue")
			{
				fragment.setProperty("TYPE", FragmentDB::TYPE__RESIDUE);
			}

			if (type_entry->getValue() == "molecule")
			{
				fragment.setProperty("TYPE", FragmentDB::TYPE__MOLECULE);
			}

			if (type_entry->getValue() == "fragment")
			{
				fragment.setProperty("TYPE", FragmentDB::TYPE__FRAGMENT);
			}
		}
	}


	void ResourceFileFragmentStorage::init()
	{
		// we are invalid until we're sure we're not...
//		valid_ = false;
		Path path;
		String filename = path.find(filename_);

		if (filename == "")
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, filename);
		}
		filename_ = filename;

		// try to open the main resource file
		ResourceFile* resource_db = new ResourceFile(filename_);

		// check for success and terminate on failure
		if (!resource_db->isValid())
		{
			delete resource_db;
			throw Exception::FileNotFound(__FILE__, __LINE__, filename_);
		}


		// copy the contents of the resource file into a tree
		if (tree_ != 0) delete tree_;
		tree_ = new ResourceEntry();			
		tree_->mergeChildrenOf(resource_db->getRoot());				

		// close the main resource file
		resource_db->close();
		delete resource_db;

		expandTree_(*tree_);

		// search for the "Fragments" entry in the main resource file
		ResourceEntry*	entry;
		entry = tree_->getRoot().findChild("Fragments");
		if (entry == 0)
		{
//			valid_ = false;
			delete tree_;

			// terminate with an exception
//			throw NoFragmentNode(__FILE__, __LINE__, filename_);
			return;
		}
		
		ResourceEntry::Iterator	frag_entry_it;
		ResourceEntry::Iterator	entry_it;
		for (frag_entry_it = ++(tree_->getEntry("/Fragments")->begin()); 
				 +frag_entry_it; ++frag_entry_it)
		{
			if (frag_entry_it->getDepth() == 2)
			{
				// create a new fragment and assign its name
				// 
				Residue* fragment = new Residue;
				fragment->setName(frag_entry_it->getKey());

				String fragment_name = (*frag_entry_it).getKey();
						
				// insert the fragment name into the corresponding lists
				Position fragment_index = addNewFragment_(fragment);
		
				// if there are no atoms in the database, something went wrong
				entry = frag_entry_it->getEntry("Atoms");
				if (entry == 0)
				{
					Log.error() << "FragmentDB: cannot find Atoms entry for " 
											<< fragment_name << std::endl;
					return;
				} 
				else	
				{
					parseAtoms_(*entry, *fragment);
				}

				// now find all the bonds for the fragment and create them
				// Fragments without bonds are legal, so we don`t complain but
				// continue
				entry = frag_entry_it->getEntry("Bonds");
				if (entry != 0)
				{
					parseBonds_(*entry, *fragment);
				}

				// now check for properties common to all variants of this
				// fragment (usually AMINO_ACID)
				// Each variant entry may also contain additional properties
				// or reset properties by specifying a "!" in front of the property
				// name
				entry = frag_entry_it->getEntry("Properties");
				if (entry != 0)
				{
					parseProperties_(*entry, *dynamic_cast<PropertyManager*>(fragment));
				}

				entry = frag_entry_it->getEntry("Connections");
				if (entry != 0)
				{
					parseConnections_(*entry, *fragment);
				}

				entry = frag_entry_it->getEntry("Type");
				if (entry != 0)
				{
					parseType_(*entry, *fragment);
				}

				// check for all aliases (given in the Names section of the db-file)
				// and insert them into the corresponding hash maps
				ResourceEntry::Iterator entry_it;
				entry = frag_entry_it->getEntry("Names");
				if (entry != 0)
				{
					String path = "/Fragments/" + fragment_name;
					for (entry_it = ++entry->begin(); +entry_it; ++entry_it)
					{
						name_to_frag_index_[entry_it->getKey()] = fragment_index;
					}
				}

				// check for possible variants of this residue type
				// (keyword Variants)
				entry = frag_entry_it->getEntry("Variants");
				if (entry != 0)
				{
					ResourceEntry::Iterator variant_it;
					Residue& original_fragment(*fragment);

					bool has_default_variant = false;

					for (variant_it = ++entry->begin(); +variant_it; ++variant_it)
					{	
						if (variant_it->getDepth() == entry->getDepth() + 1)
						{
							String variant_name = variant_it->getKey();
							Residue*	variant;
							if (variant_name == "Default")
							{
								has_default_variant = true;
								variant = new Residue(original_fragment);

								//If a default variant exists, it should take the place of the
								//basis fragment. The basis fragment itself is no longer required
								fragments_[fragment_index] = variant;
								name_to_variants_[fragment_name].push_back(fragment_index);
								name_to_frag_index_[fragment_name] = fragment_index;
							} 
							else 
							{
								variant = new Residue(original_fragment);
								variant->setName(variant_name);
								Position index = addNewFragment_(variant);

								name_to_frag_index_[variant_name] = index;
								name_to_variants_[fragment_name].push_back(index);
							}

							// Remember all variants of a certain fragment in a list.
							// This list is accessed via a hash map. It is required to 
							// determine the correct variant from given properties
							// (see getReferenceFragment(Fragment&), parseProperties_).
							for (entry_it = variant_it->begin(); +entry_it; ++entry_it)
							{
								if (entry_it->getDepth() == entry->getDepth() + 2)
								{
									const String& key = entry_it->getKey();
									if (key == "Atoms")
									{
										parseAtoms_(*entry_it, *variant);
									}
									else if (key == "Bonds")
									{
										parseBonds_(*entry_it, *variant);
									}
									else if (key == "Rename")
									{
										parseRename_(*entry_it, *variant);
									}
									else if (key == "Delete")
									{
										parseDelete_(*entry_it, *variant);
									}
									else if (key == "Connections")
									{
										parseConnections_(*entry_it, *variant);
									}
									else if (key == "Properties")
									{
										parseProperties_(*entry_it, *dynamic_cast<PropertyManager*>(variant));
									}
									else if (key == "Type")
									{
										parseType_(*entry_it, *variant);
									}
								}
							}
						}
					}

					if(has_default_variant) {
						delete fragment;
					}
				}
			}
		}

		// check for entries concerning naming standards
		entry = tree_->getEntry("/Names");
		if (entry != 0)
		{
			for (entry_it = ++entry->begin(); +entry_it; ++entry_it)
			{
				if (entry_it->getDepth() != 2) continue;
				
				// Create empty hash maps for both directions (to and from the standard).
				StringHashMap<String> map;
				standards_[entry_it->getKey() + "-" + entry_it->getValue()] = StringHashMap<String>();
				standards_[entry_it->getValue() + "-" + entry_it->getKey()] = StringHashMap<String>();
				StringHashMap<String>& name_map_to = standards_[entry_it->getKey() + "-" + entry_it->getValue()];
				StringHashMap<String>& name_map_from = standards_[entry_it->getValue() + "-" + entry_it->getKey()];
				
				// Fill those maps.
				ResourceEntry::Iterator	alias_iterator(++entry_it->begin());
				for (; +alias_iterator; ++alias_iterator)
				{
					name_map_to[alias_iterator->getKey()] = alias_iterator->getValue();
					name_map_from[alias_iterator->getValue()] = alias_iterator->getKey();
				}
			}
		}

		// check for a default naming standard
		entry = tree_->getEntry("/Defaults/Naming");
		if (entry == 0)
		{
			default_standard_ = "PDB";
		} 
		else 
		{
			default_standard_ = entry->getValue();
		}

		// OK. Everything went well, so we might consider ourselves as valid.
//		valid_ = true;

		return;
	}


}
