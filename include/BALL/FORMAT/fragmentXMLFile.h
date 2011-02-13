// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//

#ifndef BALL_FORMAT_FRAGMENTXMLFILE_H
#define BALL_FORMAT_FRAGMENTXMLFILE_H

#ifndef BALL_FORMAT_GENERICMOLFILE_H
#   include <BALL/FORMAT/genericMolFile.h>
#endif

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

		/// Read a single molecule, in this case the first variant in the file.
		Molecule* read();

		/// Read a specific molecule (variant) from the file.
		Molecule* getMoleculeForVariant(const String& variantName);

		/// alias for getMoleculeForVariant().
		Molecule* read(const String& variantName);

		bool write(const Molecule &molecule);

		bool isOpen() const;

		FragmentXMLFile& operator = (const FragmentXMLFile& rhs);

		// (no need to reimplement operators, they call read/write in the superclass)
		private:
		std::vector<String> getVariantNames();
		Atom* atomFromAtomDomForVariant(QDomElement& atomElement, const String& variant);
		bool bondFromInDomForMolecule(QDomElement& inElement, Molecule* molecule);
		QDomDocument* data_;
		Position currentVariant_;
	};
}

#endif // BALL_FORMAT_FRAGMENTXMLFILE_H
