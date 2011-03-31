// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//

#ifndef BALL_FORMAT_FRAGMENTXMLFILE_H
#include <BALL/FORMAT/fragmentXMLFile.h>
#endif

#ifndef BALL_KERNEL_ATOM_H
#include <BALL/KERNEL/atom.h>
#endif

#ifndef BALL_KERNEL_MOLECULE_H
#include <BALL/KERNEL/molecule.h>
#endif

#ifndef BALL_KERNEL_SYSTEM_H
#include <BALL/KERNEL/system.h>
#endif

#ifndef BALL_KERNEL_PTE_H
#include <BALL/KERNEL/PTE.h>
#endif

#include <set>

namespace BALL
{

	/// a bond that knows which variants it occurs in
	FragmentXMLFile::VariantedBond::VariantedBond(FragmentXMLFile::VariantedAtom* first, FragmentXMLFile::VariantedAtom* second)
	{
		first_ = first;
		second_ = second;
		first->bonds_.push_back(this);
		second->bonds_.push_back(this);
	}

	void FragmentXMLFile::VariantedBond::addOccurenceWithOrder(FragmentXMLFile::Variant& v, Bond::BondOrder order)
	{
		occursAsIn_[order].push_back(v.getId());
	}

	Bond::BondOrder FragmentXMLFile::VariantedBond::orderIn(const Variant &v)
	{
		for (Size b = 0; b < (Size)Bond::NUMBER_OF_BOND_ORDERS; b++)
		{
			for (std::list<String>::const_iterator variantName = occursAsIn_[b].begin();
					 variantName != occursAsIn_[b].end();
					 ++variantName)
			{
				if (*variantName == v.getId())
				{
					return (Bond::BondOrder)b;
				}
			}
		}

		return Bond::ORDER__UNKNOWN;
	}

	FragmentXMLFile::Variant::Variant(String id, FragmentXMLFile& container) : id_(id)
	{
		container_ = &container;
		container.variants_[id_] = this;
	}

	void FragmentXMLFile::Variant::addAtom(Atom &atom) {
		FragmentXMLFile::VariantedAtom* vAtom = new FragmentXMLFile::VariantedAtom(atom.getElement());
		container_->atomIndex_[atom.getName()] = vAtom;
		atoms_.push_back(vAtom);
		vAtom->addName(atom.getName(), *this);
		// TODO: extract alternate names
	}

	const String& FragmentXMLFile::Variant::getId() const
	{
		return id_;
	}

	FragmentXMLFile::VariantedAtom::VariantedAtom(Element el)
	{
		element_ = el;
	}

	void FragmentXMLFile::VariantedAtom::addName(String name, Variant &v)
	{
		occursAsIn_[name].push_back(v.getId());
	}

	FragmentXMLFile::VariantedBond* FragmentXMLFile::VariantedAtom::getBondTo(FragmentXMLFile::VariantedAtom& second)
	{
				for (Position i = 0; i < bonds_.size(); i++)
				{
						if ((bonds_[i]->first_ == this && bonds_[i]->second_ == &second)
								|| (bonds_[i]->second_ == this && bonds_[i]->first_ == &second))
						{
							return bonds_[i];
						}
				}
				/// bond not yet found, create a new Bond.
				FragmentXMLFile::VariantedBond* newBond = new FragmentXMLFile::VariantedBond(this, &second);
				bonds_.push_back(newBond);
				return newBond;
	}

	String FragmentXMLFile::VariantedAtom::getMajorityName()
	{
		String majorityName;
		HashMap<String, std::list<String> >::const_iterator it = occursAsIn_.begin();
		Size majoritySize = 0;
		for (; it != occursAsIn_.end(); ++it)
		{
			if (it->second.size() > majoritySize)
			{
				majoritySize = it->second.size();
				majorityName = it->first;
			}
		}
		return majorityName;
	}

	FragmentXMLFile::FragmentXMLFile()
		: GenericMolFile()
	{
		// always create the DomDocument container.
		data_ = new QDomDocument();
	}

	FragmentXMLFile::FragmentXMLFile(const String &filename, File::OpenMode open_mode)
		: GenericMolFile()
	{
		data_ = new QDomDocument;
		File::name_ = filename;
		File::open_mode_ = open_mode;
		File::is_open_ = true;
		currentVariant_ = 0;
		if (open_mode & std::ios::in) {
			// load the file if we're supposed to.
			QFile file(filename.c_str());
			if (file.open(QIODevice::ReadOnly))
			{
				// doesn't really matter if it succeeds, since an empty dom is still valid in here.
				int line, col;
				QString error;
				if (!data_->setContent(&file, &error, &line, &col))
				{
					Log.error() << "Error parsing xml in " << filename << " line "<< line << " col "<<col;
					Log.error() << ":" << String(error) << std::endl;
				}
				file.close();
			}
		}
		parse();
	}

	FragmentXMLFile::~FragmentXMLFile() {
		// free the DOM
		delete data_;

		// free out intermediate data
		std::set<VariantedAtom* > atoms;
		std::set<VariantedBond* > bonds;

		// collect atoms
		for (HashMap<String,VariantedAtom*>::Iterator it = atomIndex_.begin();
		     it != atomIndex_.end();
		     ++it)
		{
			atoms.insert(it->second);
		}
		for (std::set<VariantedAtom*>::iterator it = atoms.begin();
		     it != atoms.end();
		     ++it)
		{
			// collect each atoms bonds
			for (std::vector<VariantedBond* >::iterator bond = (*it)->bonds_.begin();
			     bond != (*it)->bonds_.end();
			     ++it)
			{
				bonds.insert(*bond);
			}
			// delete the atom
			delete *it;
		}
		// delete all bonds
		for (std::set<VariantedBond*>::iterator it = bonds.begin();
		     it != bonds.end();
		     ++it)
		{
			delete *it;
		}
		// delete the variant containers
		for (HashMap<String,Variant*>::Iterator it = variants_.begin();
		     it != variants_.end();
		     ++it)
		{
			delete it->second;
		}
	}

	bool FragmentXMLFile::isOpen() const {
		// this is used by several methods of the parent classes.
		// since here, we only operate on the in-memory representation
		// always return true
		return true;
	}

	std::vector<String> FragmentXMLFile::getVariantNames() {
		// this has a definite order: the preorder of the nodes in the document.
		std::vector<String> variants;
		QDomNodeList nodes = data_->elementsByTagName("variant");
		for (int i = 0; i < nodes.count(); i++) {
			QDomElement node = nodes.item(i).toElement();
			if (node.hasAttribute("id")) {
				variants.push_back(String(node.attribute("id")));
			}
		}
		return variants;
	}


	Atom* FragmentXMLFile::fromVariantedAtomForVariant(VariantedAtom *va, const Variant &variant)
	{
		Atom* theAtom = new Atom;
		theAtom->setName(va->getMajorityName());
		theAtom->setElement(va->element_);
		theAtom->setPosition(va->positionIn_[variant.getId()]);
		return theAtom;
	}

	Vector3 FragmentXMLFile::fromDomAtomNodeForVariant(QDomElement &element, const Variant &variant)
	{
		Vector3 posVector;
		QDomNodeList positions = element.elementsByTagName("position");
		for (int i = 0; i < positions.count(); i++)
		{
			QDomElement pos = positions.item(i).toElement();
			if (pos.attribute("variant") == variant.getId().c_str())
			{
				posVector.x = pos.attribute("x").toFloat();
				posVector.y = pos.attribute("y").toFloat();
				posVector.z = pos.attribute("z").toFloat();
			}
		}
		return posVector;
	}

	Molecule* FragmentXMLFile::fromVariant(Variant *v)
	{
		Molecule* molecule = new Molecule;
		molecule->setName(v->getId());

		HashMap<VariantedAtom* , Atom* > atomForVariantedAtom; // for creating the bonds later
		std::list<VariantedAtom* >::iterator outerAtom = v->atoms_.begin();
		for (; outerAtom != v->atoms_.end(); ++outerAtom) {
			atomForVariantedAtom[*outerAtom] = fromVariantedAtomForVariant(*outerAtom, *v);
			molecule->appendChild(*atomForVariantedAtom[*outerAtom]);
		}

		std::list<VariantedAtom* >::iterator innerAtom = outerAtom = v->atoms_.begin();
		for (; outerAtom != v->atoms_.end(); ++outerAtom)
		{
			for (; innerAtom != v->atoms_.end(); ++innerAtom)
			{
				if (*innerAtom == *outerAtom)
				{
					continue;
				}
				// bond constructor takes care of registering the bond.
				// also, WhyTF do bonds have names?
				new Bond(String(""),
				         *atomForVariantedAtom[*outerAtom],
				         *atomForVariantedAtom[*innerAtom],
				         (*outerAtom)->getBondTo(**innerAtom)->orderIn(*v));
			}
		}

		return molecule;
	}

	bool FragmentXMLFile::bondFromInDomForMolecule(QDomElement& inElement, Molecule* molecule)
	{
		QDomElement asElement = inElement.parentNode().toElement();
		QDomElement bondElement = asElement.parentNode().parentNode().toElement();
		QDomNodeList partners = bondElement.elementsByTagName("partner");
		if (partners.size() < 2) return false;
		String partnerA(partners.at(0).toElement().attribute("atom"));
		String partnerB(partners.at(1).toElement().attribute("atom"));
		Atom* first = molecule->getAtom(partnerA);
		Bond* theBond = first->createBond(* molecule->getAtom(partnerB));
		theBond->setOrder(fromDomBondOccursAsNode(asElement));

		return true;
	}

	Bond::BondOrder FragmentXMLFile::fromDomBondOccursAsNode(QDomElement &element)
	{
		Bond::BondOrder retVal;
		float bondOrder = element.attribute("order").toFloat();
		// bond order: x.5 = x, aromatic; x = x
		bool isAromatic = (bondOrder - floorf(bondOrder)) > 0.4f;
		// FIXME: this matches up with BALL::Bond::Order, but is not really future-proof.
		retVal = (Bond::BondOrder)(int)floorf(bondOrder);
		if (isAromatic) retVal = Bond::ORDER__AROMATIC;

		return retVal;
	}


	bool FragmentXMLFile::write(const Molecule &molecule)
	{
		// idea: merge molecule into data, then write everything
		// since this isn't necessarily a linear, linebased format, but interleaved.
		return false;
	}

	Molecule* FragmentXMLFile::read()
	{
		std::vector<String> variants = getVariantNames();
		if (currentVariant_ >= variants.size()) return 0;
		return fromVariant(variants_[variants[currentVariant_++]]);
	}


	void FragmentXMLFile::parse() {
		// phase 1 - prepare all variants.
		QDomNodeList nodes = data_->elementsByTagName("variant");
		for (int i = 0; i < nodes.count(); i++) {
			QDomElement node = nodes.item(i).toElement();
			if (node.hasAttribute("id")) {
				String variantID = String(node.attribute("id"));
				variants_[variantID] = new FragmentXMLFile::Variant(variantID, *this);
			}
			// TODO: preserve rest of variant metadata
		}

		// phase 2 - parse the atoms
		nodes = data_->elementsByTagName("atom");
		for (int i = 0; i < nodes.count(); i++)
		{
			QDomElement atomElement = nodes.item(i).toElement();

			QDomNodeList occursInElements = atomElement.elementsByTagName("in");
			for (int v = 0; v < occursInElements.count(); v++)
			{
				QDomElement inElement = occursInElements.item(v).toElement();
				String variantWhereAtomOccurs = String(inElement.attribute("variant"));
				if (!variants_.has(variantWhereAtomOccurs))
				{
					Log.error() << "Atom " << i << " occurs in a Variant that is not defined in the file. Skipping variant." << std::endl;
					continue;
				}
				VariantedAtom* currentAtom = findOrCreateFromDomAtomNode(atomElement);
				if (!currentAtom)
				{
					Log.error() << "Can't parse atom from DOM. Skipping Atom." << std::endl;
					continue;
				}

				variants_[variantWhereAtomOccurs]->atoms_.push_back(currentAtom);
			}
		}

		// pase 3 - parse the bonds
		nodes = data_->elementsByTagName("bond");
		for (int i = 0; i < nodes.count(); i++)
		{
			QDomElement bondElement = nodes.item(i).toElement();

			QDomNodeList partnerElements = bondElement.elementsByTagName("partner");
			if (partnerElements.count() != 2)
			{
				Log.error() << "Bond with " << partnerElements.count() << " != 2 Atoms. Skipping bond." << std::endl;
				continue;
			}
			String partner1 = String(partnerElements.item(0).toElement().attribute("atom"));
			String partner2 = String(partnerElements.item(1).toElement().attribute("atom"));
			if (atomIndex_.has(partner1) && atomIndex_.has(partner2))
			{
				VariantedBond* bond = atomIndex_[partner1]->getBondTo(*atomIndex_[partner2]);

				QDomNodeList occursInElements = bondElement.elementsByTagName("in");
				for (int v = 0; v < occursInElements.count(); v++)
				{
					QDomElement inElement = occursInElements.item(v).toElement();
					String variantWhereBondOccurs = String(inElement.attribute("variant"));
					if (!variants_.has(variantWhereBondOccurs))
					{
						Log.error() << "Bond (" << partner1 << "," << partner2 <<") occurs in a Variant ("
						            << variantWhereBondOccurs <<") that is not defined in the file. Skipping variant." << std::endl;
						continue;
					}
					QDomElement asElement = inElement.parentNode().toElement();
					bond->addOccurenceWithOrder(*variants_[variantWhereBondOccurs], fromDomBondOccursAsNode(asElement));
				}
			}
			else
			{
				Log.error() << "Partner not found for Bond, "<< partner1 << " or "<< partner2 <<" missing. Skipping bond." << std::endl;
				continue;
			}
		}
	}

	FragmentXMLFile::VariantedAtom* FragmentXMLFile::findOrCreateFromDomAtomNode(QDomElement& atomElement)
	{
		QDomNodeList occursAsElements = atomElement.elementsByTagName("as");
		for (int n = 0; n < occursAsElements.count(); n++)
		{
			String atomName(occursAsElements.item(n).toElement().attribute("name"));
			if (atomIndex_.has(atomName))
			{
				Log.info() << "Atom "<< atomName << " already parsed, returning cached copy." << std::endl;
				return atomIndex_[atomName];
			}
		}

		// none found, create.
		// determine element from periodic table
		QDomNodeList elementTags = atomElement.elementsByTagName("element");
		if (elementTags.count() != 1)
		{
			Log.info() << "No ("<<elementTags.count()<< ") element Tag for Atom " << ". Skipping." << std::endl;
			return 0; // can't create this atom :-(
		}
		String element(elementTags.at(0).toElement().text());
		VariantedAtom* theAtom = new VariantedAtom(PTE[element]);
		for (int n = 0; n < occursAsElements.count(); n++)
		{
			String atomName(occursAsElements.item(n).toElement().attribute("name"));
			atomIndex_[atomName] = theAtom;
		}

		// TODO: preserve rest of Atom Metadata (?)

		return theAtom;
	}

	FragmentXMLFile::VariantedAtom* FragmentXMLFile::findOrCreateFromAtom(Atom &atom)
	{
		if (atomIndex_.has(atom.getName()))
		{
			return atomIndex_[atom.getName()];
		}

		// TODO: import multiple names (named properties?)
		// TODO: import rest of atom metadata?

		VariantedAtom* theAtom = new VariantedAtom(atom.getElement());
		atomIndex_[atom.getName()] = theAtom;
		return theAtom;
	}
} // namespace BALL
