// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//

#ifndef BALL_STRUCTURE_FRAGMENTDB_RESOURCEFILEFRAGMENTSTORAGE_H
#define BALL_STRUCTURE_FRAGMENTDB_RESOURCEFILEFRAGMENTSTORAGE_H

#include<BALL/STRUCTURE/FRAGMENTDB/fragmentStorage.h>
#include<BALL/STRUCTURE/FRAGMENTDB/nameFragmentQuery.h>
#include<BALL/STRUCTURE/FRAGMENTDB/connection.h>

#include<BALL/FORMAT/resourceFile.h>
#include<BALL/CONCEPT/property.h>
#include<BALL/DATATYPE/stringHashMap.h>
#include<BALL/COMMON/global.h>

namespace BALL
{

	class ResourceFileFragmentStorage : public FragmentStorage
	{
		public:
			explicit ResourceFileFragmentStorage(const ResourceFileFragmentStorage& other, bool deep = true);
			explicit ResourceFileFragmentStorage(const String& path_to_root_resourcefile = "fragments/Fragments.db");

			bool query(FragmentQuery&);

			String getVersionTag();

		private:
			ResourceEntry*	tree_;


		// ----------- stuff from ex-FragmentDB.h -------------

		/** @name	Parse functions
				These functions parse different sections of the fragment DB resource
				tree and translate everything into the correct data structures.
		*/
		/// @{

		/**	Parses the Atoms entry and creates the atoms of the fragment
		*/
		void parseAtoms_(ResourceEntry& entry, Fragment& fragment);

		/**	Parses the Bonds entry and creates the bonds of the fragment
		*/
		void parseBonds_(ResourceEntry& entry, Fragment& fragment);

		/**	Parses the properties of a fragment/variant and sets the corresponding properties.
				Properties are set or reset (if the property name starts with "!") for the current 
				fragment. All properties of fragment and residue are recognized, if they
				are written exactly as in the header file (case insensitive) and set.
				Unknown properties are set as name/value pairs as bool properties and set to 
				<b>  true </b>.
		*/
		void parseProperties_(ResourceEntry& entry, PropertyManager& property_man);

		/**	Parses the Delete section.
				All atoms given in this section are removed from the fragment.
		*/
		void parseDelete_(ResourceEntry& entry, Fragment& fragment);

		/**	Parses the Rename section.
				All atoms given in this section are renamed to the given new name.
		*/
		void parseRename_(ResourceEntry& entry, Fragment& fragment);

		/** Parses the Connections section
		    There, possible inter-fragment bonds are stored.
		 */
		void parseConnections_(ResourceEntry&  entry, Fragment& fragment);

		/** Parses Type information
		 */
		void parseType_(ResourceEntry&  entry, Fragment& fragment);

		/// @}

		/** Add a new fragment pointer to the database (while parsing) */
		Position addNewFragment_(Residue* fragment);

		/**	Expands all include directives in the resource file.
				This method calls expandFirst_ until it returns true.	
		*/
		void expandTree_(ResourceEntry& root_entry);
	 
		/**	Expands the first occuring include directive.
				If no include directive is found, <b>  false </b> is returned, otherwise <b>  true </b>.
				@exception Exception::FileNotFound if the file is not found in the BALL_DATA_PATH
		*/
		bool expandFirst_(ResourceEntry& root_entry);
		
		void init();
		
				// The filename of the master fragment file.
		String 					filename_;

		// The naming standard we default to.
		String					default_standard_;

		// An array containing all allocated residues.
		std::vector<Residue*>						fragments_;

		typedef StringHashMap<String>		NameMap;

		// Maps a fragment name back to a path in the database
		NameMap													name_to_path_;

		// Maps a fragment name back to the array index in fragments_
		StringHashMap<Position>					name_to_frag_index_;

		// Maps all variants of a specific fragment back to array indices.
		StringHashMap<list<Position> >	name_to_variants_;

		// Contains the naming standards as a nested map.
		StringHashMap<NameMap>					standards_;
	};
}

#endif // RESOURCEFILEFRAGMENTSTORAGE_H
