// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//

#ifndef BALL_FORMAT_FRAGMENTXMLFILE_H
#define BALL_FORMAT_FRAGMENTXMLFILE_H

#include <BALL/FORMAT/genericMolFile.h>
#include <BALL/DATATYPE/hashMap.h>
#include <BALL/KERNEL/bond.h>
#include <BALL/KERNEL/PTE.h>

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
		virtual Molecule* read();

		virtual bool write(const Molecule &molecule);

		virtual bool isOpen() const;

		FragmentXMLFile& operator = (const FragmentXMLFile& rhs);

		// (no need to reimplement operators, they call read/write in the superclass)
		private:
		class Variant;
		class VariantedBond;
		/// translates to several Atoms
		class VariantedAtom {
			public:
				/// creates a new VariantedAtom. Only the element is fixed over its lifetime.
				VariantedAtom(Element el);
				void addName(String name, Variant& v);
				void addPosition(Vector3, Variant& v);
				String getMajorityName();
				VariantedBond* getBondTo(VariantedAtom& partner);
				friend class VariantedBond;
				friend class FragmentXMLFile;
			private:
				HashMap<String, Vector3> positionIn_;
				Element element_;
				std::vector<VariantedBond*> bonds_;
				HashMap<String, std::list<String> > occursAsIn_;
		};

		/// translates to several Bonds
		class VariantedBond {
			public:
				VariantedBond(VariantedAtom* first, VariantedAtom* second);
				void addOccurenceWithOrder(Variant& v, Bond::BondOrder order);
				Bond::BondOrder orderIn(const Variant& v);
				friend class VariantedAtom;
			protected:
				VariantedAtom* first_;
				VariantedAtom* second_;
				std::list<String> occursAsIn_[Bond::NUMBER_OF_BOND_ORDERS];
		};

		/// translates to a Molecule
		class Variant {
			public:
				Variant(String id, FragmentXMLFile& container);
				void addAtom(Atom& atom);
				const String& getId() const;
				friend class FragmentXMLFile;
			private:
				FragmentXMLFile* container_;
				String id_;
				std::list<VariantedAtom*> atoms_;
		};

		/// the index of all atoms, indicated by their names. in case of alternate names,
		/// multiple names may point to the same atom.
		HashMap<String,VariantedAtom*> atomIndex_;
		HashMap<String,Variant*> variants_;

		/// returns the names of all known variants in the file
		std::vector<String> getVariantNames();
		void parse();
		VariantedAtom* findOrCreateFromDomAtomNode(QDomElement& atomElement);
		VariantedAtom* findOrCreateFromAtom(Atom& atom);
		Bond::BondOrder fromDomBondOccursAsNode(QDomElement& element);
		Vector3 fromDomAtomNodeForVariant(QDomElement& element, const Variant& variant);
		Atom* fromVariantedAtomForVariant(VariantedAtom* va, const Variant& variant);
		Molecule* fromVariant(Variant* v);
		bool bondFromInDomForMolecule(QDomElement& inElement, Molecule* molecule);
		QDomDocument* data_;
		Position currentVariant_;
	};
}

#endif // BALL_FORMAT_FRAGMENTXMLFILE_H
