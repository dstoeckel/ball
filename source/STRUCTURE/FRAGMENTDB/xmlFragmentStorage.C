#include <BALL/STRUCTURE/FRAGMENTDB/xmlFragmentStorage.h>
#include <BALL/STRUCTURE/FRAGMENTDB/fragmentQuery.h>
#include <BALL/STRUCTURE/FRAGMENTDB/nameFragmentQuery.h>
#include <BALL/STRUCTURE/FRAGMENTDB/propertyFragmentQuery.h>
#include <BALL/SYSTEM/directory.h>
#include <BALL/SYSTEM/path.h>
#include <BALL/FORMAT/fragmentXMLFile.h>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#define BALL_DEBUG_XMLFRAGMENTSTORAGE

#ifdef BALL_DEBUG_XMLFRAGMENTSTORAGE
# define DEBUG(a) Log.info() << a << std::endl
#else
# define DEBUG(a)
#endif

namespace BALL {

const String XMLFragmentStorage::XML_FRAGMENT_FILE_SUFFIX = ".fdb";

	XMLFragmentStorage::XMLFragmentStorage(const String& our_naming_convention,
	                                       const String &path_to_root_directory)
	: our_convention_(our_naming_convention)
	{
		Path path;
		Directory dir(path.find(path_to_root_directory));
		if (!dir.isValid())
		{
			DEBUG("Couldn't find directory " << path_to_root_directory);
			return;
		}
		DEBUG("Loading FragmentXML from " << dir.getPath() );
		String filename;
		unsigned int fragment_number = 0;
		while (dir.getNextEntry(filename))
		{
			DEBUG("Considering " << filename);
			if (!filename.hasSuffix(XML_FRAGMENT_FILE_SUFFIX)) continue;
			DEBUG("Loading " << filename);
			FragmentXMLFile file(dir.getPath() + "/" + filename);
			std::vector<String> variants = file.getVariantNames();
			for (std::vector<String>::const_iterator it = variants.begin();
			     it != variants.end(); ++it)
			{
				Residue* res = file.residueForVariant(*it);
				if (by_name_.has(res->getFullName())) {
					Log.error() << "Fragment " << *it << " redefined in " << filename
					            << ", ignoring." << std::endl;
					continue;
				}
				buildTranslationTable_(res);
				Residue* our_res = makeCleanResidueCopy(res);
				// store
				fragments_.push_back(our_res);
				// index by full name
				by_name_[our_res->getFullName()] = fragment_number;
				// index by simple name
				by_parent_[our_res->getName()].push_back(fragment_number);
				// index by unnamed properties
				by_props_[our_res->getBitVector().getUnsignedLong()].push_back(fragment_number);
				fragment_number++;
			}
			DEBUG("Done with " << filename << " got " << fragment_number << " fragments so far.");
		}
		DEBUG("Loaded " << fragment_number << " fragments from " << path_to_root_directory << "*"
					<< XML_FRAGMENT_FILE_SUFFIX );
	}

	XMLFragmentStorage::XMLFragmentStorage(const XMLFragmentStorage &other, bool deep)
	: our_convention_(other.our_convention_)
	{
	 // TODO: copy stuff.
	}

	bool XMLFragmentStorage::query(FragmentQuery &query)
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

			if (query.getMaxResults() == 1)
				{
					if (by_name_.has(q->getFragmentName())) {
						boost::shared_ptr<Residue> frag_copy_p(
									new Residue(*fragments_[by_name_[q->getFragmentName()]])
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

					StringHashMap<list<size_t> >::ConstIterator to_find = by_parent_.find(q->getFragmentName());
					if (to_find == by_parent_.end())
						{
							// TODO: maybe we want to try exact matching here then...
							return false;
						}

					list<size_t>::const_iterator it = (*to_find).second.begin();
					const list<size_t>::const_iterator end_it = (*to_find).second.end();

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
		else if (query.selectsOn(FragmentQuery::QueryFragmentProperties))
		{
			PropertyFragmentQuery* q;
			try {
				q = boost::any_cast<PropertyFragmentQuery*>(query.getSelectorDetail(FragmentQuery::QueryFragmentProperties));
			} catch (boost::bad_any_cast& error) {
				Log.error() << "ResourceFileFragmentStorage: Unable to process query: " << query.toString() << std::endl;
				return false;
			}
			const BitVector& tplate(q->getPropertyManager().getBitVector());

			list< size_t >::iterator candidate = by_props_[tplate.getUnsignedLong()].begin();
			list< size_t >::iterator end = by_props_[tplate.getUnsignedLong()].end();

			for (; candidate != end; ++candidate)
			{
				Residue* fragment = fragments_[*candidate];
				// I'll still try the full vector, since the uLong representation
				// might have been truncated.
				if ((fragment->getBitVector() & tplate) == tplate)
				{
					// all bits in the template are also present in the fragment
					// now check the named properties. (TODO: test equality here?)
					bool namedPropsOK = true;
					NamedPropertyIterator namedProp = q->getPropertyManager().beginNamedProperty();
					NamedPropertyIterator end = q->getPropertyManager().endNamedProperty();
					for (; (namedProp != end) && namedPropsOK; ++namedProp)
					{
						if (!(fragment)->hasProperty(namedProp->getName()))
						{
							namedPropsOK = false;
						}
					}
					if (namedPropsOK)
					{
						boost::shared_ptr<Residue> frag_copy_p(
							new Residue(*fragment)
						);
						query.addResult(frag_copy_p);
					}
				}
			}
		}
	}

	String XMLFragmentStorage::getVersionTag()
	{
	// TODO: hashing
	}

	void XMLFragmentStorage::buildTranslationTable_(Residue* res)
	{
		NameMap name_in_convention;
		NamedPropertyIterator prop = res->beginNamedProperty();
		while (prop != res->endNamedProperty())
		{
			if (boost::starts_with(prop->getName(),"NAME:"))
			{
				std::vector<std::string> convention;
				boost::split(convention, prop->getName(), boost::is_any_of(":"));
				// I'm assuming names are present and valid, since I expect
				// the document to be valid to the schema.
				name_in_convention[convention[1]] = prop->getString();
			}
			++prop;
		}
		AtomIterator atom = res->beginAtom();
		while (atom != res->endAtom())
		{
			prop = atom->beginNamedProperty();
			NameMap atom_in_convention;
			while (prop != atom->endNamedProperty())
			{
				if (boost::starts_with(prop->getName(),"NAME:"))
				{
					std::vector<std::string> convention;
					boost::split(convention, prop->getName(), boost::is_any_of(":"));
					// I'm assuming names are present and valid, since I expect
					// the document to be valid to the schema.
					if (!name_in_convention.has(convention[1])) {
						Log.warn() << "Atom " << atom->getName() << " in " << res->getName()
											 << " defines name in convention " << convention[1] << ", but "
											 << res->getName() << " doesn't. Ignoring this name." << std::endl;
						++prop;
						continue;
					}
					atom_in_convention[convention[1]] = prop->getString();
				}
				++prop;
			}
			// all collected, now build table. consider all pairings:
			NameMap::Iterator i = atom_in_convention.begin();
			for (; i != atom_in_convention.end() ; ++i)
			{
				for (NameMap::Iterator j = i; j!= atom_in_convention.end() ; ++j)
				{
					if (j == i) continue;
					registerName_(j->first,name_in_convention[j->first],j->second,
					              i->first,name_in_convention[i->first],i->second);
				}
			}
			++atom;
		}
	}

	void XMLFragmentStorage::registerName_(const String& convention_a,
	                                       const String& residuename_a,
	                                       const String& atomname_a,
	                                       const String& convention_b,
	                                       const String& residuename_b,
	                                       const String& atomname_b)
	{
		// if we ever devise a cooler way of doing this, it should be
		// sufficient to change it here.
		NameMap& a_to_b(translation_tables_[convention_a+"-"+convention_b]);
		// TODO: check if using the first definition is sufficient here.
		a_to_b[residuename_a+":"+atomname_a] = residuename_b+":"+atomname_b;
		NameMap& b_to_a(translation_tables_[convention_b+"-"+convention_a]);
		b_to_a[residuename_b+":"+atomname_b] = residuename_a+":"+atomname_a;
	}

	class RenameAtoms : public UnaryProcessor<Atom>
	{
		public:
		String our_convention_;
		RenameAtoms(const String& convention) : our_convention_(convention) {}
		Processor::Result operator ()(Atom& atom)
		{
			atom.setName(atom.getProperty("NAME:"+our_convention_).getString());
			atom.clearPropertiesContaining("NAME:");
		}
	};

	Residue* XMLFragmentStorage::makeCleanResidueCopy(Residue *theirs)
	{
		Residue* ours = new Residue(*theirs);
		ours->setName(ours->getProperty("NAME:"+our_convention_).getString());
		ours->clearPropertiesContaining("NAME:");
		// and sweep up after the atoms.
		RenameAtoms broom(our_convention_);
		ours->apply(broom);
		return ours;
	}
}
