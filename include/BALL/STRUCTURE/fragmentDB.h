// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//

#ifndef BALL_STRUCTURE_FRAGMENTDB_H
#define BALL_STRUCTURE_FRAGMENTDB_H

#ifndef BALL_DATATYPE_STRINGHASHMAP_H
#	include <BALL/DATATYPE/stringHashMap.h>
#endif

#ifndef BALL_KERNEL_RESIDUE_H
#	include <BALL/KERNEL/residue.h>
#endif

#ifndef BALL_KERNEL_MOLECULE_H
#	include <BALL/KERNEL/molecule.h>
#endif

#ifndef BALL_KERNEL_BOND_H
#	include <BALL/KERNEL/bond.h>
#endif

#ifndef BALL_STRUCTURE_RECONSTRUCTFRAGMENTPROCESSOR_H
#	include <BALL/STRUCTURE/reconstructFragmentProcessor.h>
#endif

#include <BALL/STRUCTURE/FRAGMENTDB/fragmentQuery.h>
#include <BALL/STRUCTURE/FRAGMENTDB/connection.h>

#include <vector>
#include <list>
#include <boost/smart_ptr/shared_ptr.hpp>

namespace BALL 
{

	class ResourceEntry;
	class FragmentStorage;

	/**	FragmentDB - fragment database class.
			The fragment database is used to store commonly
			used subunits as amino acids, sugars and the like.
			They are entered in a special format described below.
			The main resource file resides under 
			<a href="../../../data/fragments/Fragments.db">data/fragments/Fragments.db</a>.  \par
	\ingroup StructureMiscellaneous		
	*/
	class BALL_EXPORT FragmentDB 
	{
		public:

		BALL_CREATE_DEEP(FragmentDB)
	
		/**	@name	Enums 
		*/
		//@{

		/**	Fragment types 
		*/
		enum FragmentTypes 
		{
			TYPE__UNKNOWN	= -1,
			TYPE__FRAGMENT,
			TYPE__RESIDUE,
			TYPE__MOLECULE
		};
		//@}

		/**	@name	Type Definitions 
		*/
		//@{

		/**	Type definition for the fragment type 
		*/
		typedef short Type;
		/** A hash map used to convert one atom naming convention to another */
		typedef StringHashMap<String>	NameMap;

		//@}
		/**	@name	Exceptions
		*/
		//@{
		
		/**	No fragment node found.
				This exception is thrown by  \link init init \endlink  if the resource database
				does not contain a <tt>Fragments</tt> entry.
		*/
		class BALL_EXPORT NoFragmentNode
			:	public Exception::GeneralException
		{
			public:
			NoFragmentNode(const char* file, int line, const string& filename);
			~NoFragmentNode() throw() {}

			protected:
			string filename_;
		};

		//@}
		/**	@name Constructors and destructors 
		*/
		//@{
	
		/**	Creates a default but invalid FragmentDB instance.
		*/
		FragmentDB();

		/**	Creates a FragmentDB object and reads the contents of <tt>filename</tt>.
		 		If filename is an empty string, the default value "fragments/Fragments.db" is used.
				@exception Exception::FileNotFound if the file is not found in the BALL_DATA_PATH
		*/
		FragmentDB(const String& filename);

		/**	Copy constructor.
		*/
		FragmentDB(const FragmentDB& db, bool deep = true);	

		/// Assignment  operator 
		FragmentDB& operator = (const FragmentDB& db);

		/**	Destructor.
		*/
		virtual ~FragmentDB();

		//@}

		/**@name	Inspectors and mutators
		*/
		//@{

		/**	Checks whether a specified fragment is known to the fragment database.
		*/
		bool has(const String& fragment_name) const;
		
		/** Queries the Database backends. */
		bool query(FragmentQuery&) const;

		///
		const std::vector<Residue*>& getFragments() const;
		
		/**	Return a fragment.
		*/
		FragmentDB::Type getFragmentType(const String& fragment_name) const;

		/**	Return a list containing all variant names.
		*/
		list<String> getVariantNames(const String& name) const;

		/** Add an externally constructed store to the Queriable backends
		 */
		void addStore(boost::shared_ptr<FragmentStorage>& store);

		/**	Return a fragment.
		    @deprecated Since the new FragmentDB architecture does not allow direct access
		    to the backends raw data any more, this function cannot give you a pointer
		    owned by the FragmentDB any more.
		    Using this function is semantically equivalent to calling getFragmentCopy(),
		    you will have to clean up the Fragment* returned after you're done with it.
		    Consider using query() instead, where smart_ptrs will take care of this for you.
		*/
		BALL_DEPRECATED const Fragment* getFragment(const String& fragment_name) const;

		/**	Return a reference fragment.
				This method returns a standard template of a given fragment or a NULL pointer
				if the fragment is not known. The first criterion is the fragment name.
				If there exist multiple variants of the fragment, the correct variant is chosen 
				according to the properties set in <tt>fragment</tt>.
		*/
		boost::shared_ptr<Residue> findReferenceFragment(const Fragment& fragment) const;

		/** Returns a reference fragment.
				Returns an independent copy of the reference fragment found by #findReferenceFragment() ,
				if any.
				@deprecated Since the new FragmentDB architecture does not allow direct access
		    to the backends raw data any more, this function cannot give you a pointer
		    owned by the FragmentDB any more.
		    @see findReferenceFragment()
		    @see getFragment()
		*/
		BALL_DEPRECATED const Fragment* getReferenceFragment(const Fragment& fragment) const;

		/**	Return a residue.
		*/
		BALL_DEPRECATED const Residue* getResidue(const String& fragment_name) const;

		/**	Return a copy of a fragment.
				If a fragment with name <tt>fragment_name</tt> exists in the
				fragment database, a copy is created and returned. 
				Otherwise, a null pointer is returned. 
				Take care to destruct the copy again to avoid memory leaks.
				@return a pointer to the copied fragment or 0
				@param	fragent_name the name of the fragment in the database
		*/
		Fragment* getFragmentCopy(const String& fragment_name) const;

		/**	Return a copy of a fragment as a molecule.
				If a fragment with name <tt>fragment_name</tt> exists in the
				fragment database, a copy is created, inserted into a new molecule, and returned. 
				Otherwise, a null pointer is returned. 
				Take care to destruct the copy again to avoid memory leaks.
				@return a pointer to the copied fragment or 0
				@param	fragent_name the name of the fragment in the database
		*/
		Molecule* getMoleculeCopy(const String& fragment_name) const;

		/**	Return a copy of a residue.
				If a fragment with name <tt>fragment_name</tt> exists in the
				fragment database, a copy is created and returned as a residue. 
				Otherwise, a null pointer is returned. Take care to destruct the copy again
				to avoid memory leaks.
				@return a pointer to the copied fragment or 0
				@param	fragent_name the name of the fragment in the database
		*/
		Residue* getResidueCopy(const String& fragment_name) const;

		/**	Return the default naming standard
		*/
		const String&	getDefaultNamingStandard() const;

		/**	Return a hash map containing all naming maps.
		*/
		StringHashMap<NameMap>&	getNamingStandards();

		/**
		 * Return the naming standard given by std
		 *
		 * @return A StringHashMap that maps atom names to atom names
		 * @throw StringHashMap<String>::IllegalKey if std is not a valid naming standard
		 */
		const StringHashMap<String>& getNamingStandard(const String& std) const;

		/**
		 * Return a vector of available naming standards
		 */
		std::vector<String> getAvailableNamingStandards() const;

		//@}

		/**@name	Debugging and diagnostics 
		*/
		//@{

		/**	
		*/
		bool isValid() const;

		//@}
		/**	@name	Processors defined in the fragment DB
		*/
		//@{

		/**	Name normalization processor.
				This class is used to adopt all names in a molecular system 
				to a given naming standard (usually the PDB standard).
		*/
		class BALL_EXPORT NormalizeNamesProcessor 
			: public UnaryProcessor<Fragment>
		{
		
			public:

			/**	@name	Constructors and Destructors
			*/
			//@{

			/**	Default constructor
			*/
			NormalizeNamesProcessor();
		
		
			/**	Constructor
			*/
			NormalizeNamesProcessor(FragmentDB& db);

			/**	Destructor
			*/
			virtual ~NormalizeNamesProcessor();
		
			//@}
			/**@name	Inspectors and Mutators
			*/
			//@{
			
			/**	Bind the processor to a fragment database.
			*/
			void setFragmentDB(FragmentDB& db);

			/**
			*/
			void setNamingStandard(const String& naming_standard);
		
			/**	Retrieve the current naming standard
			*/
			const String& getNamingStandard();

			/**	Try to match a name in one of the maps
			*/
			bool matchName(String& res_name, String& atom_name, const NameMap& map) const;
	 
			//@}
			/**@name	Processor specific methods
			*/
			//@{
			
			/**	Start method
			*/
			virtual bool start();

			/**	Finish method	
			*/
			virtual bool finish();
		
			/**	Application method
			*/
			virtual Processor::Result operator () (Fragment& fragment);

			//@}

			private:
				enum CountingMode { ADD, OVERWRITE };
				String getSuffix_(const Fragment* frag) const;
				bool doMatch_(String& res_name, const String& res_name_suffix, String& atom_name, const NameMap& map) const;
				void countHits_(HashMap<NameMap*, Index>& maps, const std::list<Fragment*>& frags);
				void countHits_(HashMap<NameMap*, Index>& maps, const Fragment* frag, CountingMode mode = OVERWRITE);
				const NameMap* getBestMap_(const HashMap<NameMap*, Index>& maps) const;
				void normalizeFragments_(const NameMap* map, const std::list<Fragment*>& frags);
				void normalizeFragment_ (const NameMap* map, Fragment* frag);

			String								naming_standard_;

			FragmentDB*						fragment_db_;
		
			std::list<Fragment*>	fragments_;

		};


		/**	Bond creation processor
		*/
		class BALL_EXPORT BuildBondsProcessor 
			: public UnaryProcessor<Fragment> 
		{

			public:

			/**	@name Type definitions	
			*/
			//@{
			///
			typedef std::list<Connection> ConnectionList;
			//@}


			/** @name	Constructors and Destructors
			*/
			//@{
		
			///	
			BuildBondsProcessor();
			
			///
			BuildBondsProcessor(const FragmentDB& db);

			/// 
			virtual ~BuildBondsProcessor();
			//@}

			/**	@name	Processor-related methods 
			*/
			//@{

			///
			virtual bool finish();

			///
			virtual bool start();

			///
			virtual Processor::Result operator () (Fragment& fragment);
			//@}

			/**	@name	Accessors
			*/
			//@{
			
			/// Return the number of bonds built during the last application.
			Size getNumberOfBondsBuilt();

			/// Set the fragment database
			void setFragmentDB(const FragmentDB& fragment_db);

			//@}

			/**	@name	Bond building methods 
			*/
			//@{

			/**	Build all bonds in a fragment.
					This method builds all bonds that are contained
					in the template.
					@return the number of bonds built
			*/
			Size buildFragmentBonds(Fragment& fragment);

			/**	Build all bonds in a fragment according to a manually supplied
					template.
					This method builds all bonds that are contained
					in manually provided template.
					@return the number of bonds built
					@exception Exception::TooManyBonds if an atom would be assigned too many bonds
			*/
			Size buildFragmentBonds(Fragment& fragment, const Fragment& tplate);

			/**	Build all possible bonds between two fragments.
					This method builds all bonds that are allowed by
					the <b>Connections</b> entries in a resource database.
					@return the number of bonds built
					@exception Exception::TooManyBonds if an atom would be assigned too many bonds
			*/
			Size buildInterFragmentBonds(Fragment& first, Fragment& second) const;

			//@}

			protected:

			/**	Build a connection between two atoms, if possible
					@exception Exception::TooManyBonds if an atom would be assigned too many bonds
			*/
			bool buildConnection_(Connection& con1, Connection& con2);
		
			/**	A pointer to the fragment database 
			*/
			FragmentDB*			fragment_db_;
			
			/**	A list of all fragments.
					This list is constructed incrementally by the operator ()
					and is used by finish() to create the inter-fragment bonds
			*/
			std::list<Fragment*>	fragment_list_;

			/*_	The number of bonds built.
					This value is reset in the start method, so each application of 
					the processor, so <TT>  getNumberOfBuiltBonds </TT> always returns
					the number of bonds built in the last application.
			*/
			Size	bonds_built_;

			/*_	The list of connects between(!) fragments still to be examined
			*/
			ConnectionList	connections_;
		};

		//@}
		/**	@name	Public Variables
		*/
		//@{

		/**	The standard name normalization processor
		*/
		NormalizeNamesProcessor		normalize_names;

		/**	The standard hydrogen adder
		*/
		ReconstructFragmentProcessor	add_hydrogens;

		/**	The standard bond builder
		*/
		BuildBondsProcessor				build_bonds;

		//@}

		private:

		// The status of the FragmentDB
		bool	valid_;
		typedef std::set<boost::shared_ptr<FragmentStorage> > FragmentStoreSet;
		FragmentStoreSet stores_;

	};
  
} // namespace BALL 


#endif // BALL_STRUCTURE_FRAGMENTDB_H
