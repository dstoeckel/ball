// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//

#ifndef BALL_FORMAT_FRAGMENTXMLFILE_H
#define BALL_FORMAT_FRAGMENTXMLFILE_H

#include <BALL/FORMAT/genericMolFile.h>
#include <BALL/DATATYPE/stringHashMap.h>
#include <BALL/KERNEL/bond.h>
#include <BALL/KERNEL/PTE.h>
#include <BALL/CONCEPT/property.h>
#include <BALL/COMMON/global.h>

#include <QtXml>

namespace BALL
{
	/** Fragment XML file class.
		Support for reading fragmentXML files, one variant at a time.
		\ingroup StructureFormats
	*/
	class BALL_EXPORT FragmentXMLFile
				: public GenericMolFile
	{
		public:

		/**	@name	Constructors and Destructors
		*/
		//@{
		/// Default constructor.
		FragmentXMLFile();
		/// Detailed constructor.
		FragmentXMLFile(const String& filename, File::OpenMode open_mode = std::ios::in);
		/// Destructor.
		virtual ~FragmentXMLFile();
		//@}

		/** @name IO Operations
		 */
		//@{
		/** Read a single molecule.
		 *  Returns a molecule for each variant in the file.
		 *  Subsequent calls iterate through the variants
		 */
		virtual Molecule* read();

		/** Write a single molecule.
		 *  Unsupported, will print an error and do nothing.
		 */
		virtual bool write(const Molecule& molecule);

		/** Check whether the file is open.
		 */
		virtual bool isOpen() const;

		/** Load a System from this file.
		 *  Specialty: Tells BALLView to skip name normalization and fragment reconstruction.
		 *  @see MolecularStructure::addComposite_
		 */
		virtual GenericMolFile& operator>> (System& system);
		//@}

		/** @name Format specific functions
		 */
		//@{
		/** Performs a validation of the file against the XML Schema.
		 */
		bool validate();
		//@}

		FragmentXMLFile& operator = (const FragmentXMLFile& rhs);

		// (no need to reimplement other operators, they call our read/write from the superclass)
		private:
		/// returns the names of all known variants in the file
		std::vector<String> getVariantNames();
		/** returns a Molecule for a given Variant.
		 *  the return value is not cleaned up on destruction.
		 */
		Molecule* moleculeForVariant(const String&);
		/** returns a Residue for a given Variant.
		 *  this gives an internal pointer, which IS DELETED on destruction.
		 */
		Residue* residueForVariant(const String&);
		/// main dom interpreation routine.
		void parse();
		/// parses a property tag into a PropertyManager (Residue, Atom or Bond)
		void parseProperty(QDomElement&, PropertyManager&);
		/// parses a new atom into the Residue from an <atom> QDomElement
		void parseAtom(QDomElement&, Residue*);
		/// parses a new Bond into the Residue. Atoms need to be present already.
		void parseBond(QDomElement&, Residue*);
		/// parses residue names from a <variant> tag out into properties.
		void parseNames(QDomElement&, PropertyManager&);
		/// records a name in different conventions into named properties of the object.
		void registerNameForConventions(String& name, String& conventions, PropertyManager& object);
		/// adds all "NAME:Convention"="Name" properties to the map, keyed by Name.
		void collectNameForConventions(PropertyManager& object, StringHashMap<String>& names);
		/// 
		void parseConnections(QDomElement&, Atom*);
		
		Bond::Order parseBondOrder(const QString&);
		QDomDocument* data_;
		Position currentVariant_;
		std::vector< Residue* > residues_;
		StringHashMap< Position > residue_by_name_;
		StringHashMap< StringHashMap< String > > translation_tables_;
	};
}

#endif // BALL_FORMAT_FRAGMENTXMLFILE_H
