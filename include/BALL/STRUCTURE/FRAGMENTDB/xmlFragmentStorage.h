// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//

#ifndef BALL_STRUCTURE_FRAGMENTDB_XMLFRAGMENTSTORAGE_H
#define BALL_STRUCTURE_FRAGMENTDB_XMLFRAGMENTSTORAGE_H

#include<BALL/STRUCTURE/FRAGMENTDB/fragmentStorage.h>

#include<BALL/DATATYPE/stringHashMap.h>
#include<BALL/KERNEL/residue.h>

namespace BALL
{

	class XMLFragmentStorage : public FragmentStorage
	{
		public:
			static const String XML_FRAGMENT_FILE_SUFFIX;
			explicit XMLFragmentStorage(const XMLFragmentStorage& other, bool deep = true);
			explicit XMLFragmentStorage(const String& our_naming_convention = "BALL",
			                            const String& path_to_root_directory = "fragments_xml");

			bool query(FragmentQuery&);

			String getVersionTag();

		private:
			typedef StringHashMap<String> NameMap;
			void buildTranslationTable_(Residue* res);
			void registerName_(const String& convention_a,
			                   const String& residuename_a,
			                   const String& atomname_a,
			                   const String& convention_b,
			                   const String& residuename_b,
			                   const String& atomname_b);
			Residue* makeCleanResidueCopy(Residue* theirs);
			String our_convention_;
			std::vector<Residue* > fragments_;
			StringHashMap< list < size_t > > by_parent_;
			StringHashMap< size_t > by_name_;
			StringHashMap< NameMap > translation_tables_;
			HashMap< unsigned long , list< size_t > > by_props_;
			void dumpTables_();
	};

}
#endif // XMLFRAGMENTSTORAGE_H
