// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//

#ifndef BALL_FORMAT_FRAGMENTXMLFILE_H
#include <BALL/FORMAT/fragmentXMLFile.h>
#endif

#ifndef BALL_KERNEL_ATOM_H
#include <BALL/KERNEL/atom.h>
#endif

#ifndef BALL_KERNEL_BOND_H
#include <BALL/KERNEL/bond.h>
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

namespace BALL
{
	FragmentXMLFile::FragmentXMLFile()
		: GenericMolFile()
	{
		// always create the DomDocument container.
		data_ = new QDomDocument();
	}

	FragmentXMLFile::FragmentXMLFile(const String &filename, File::OpenMode open_mode)
	{
		FragmentXMLFile();
		File::name_ = filename;
		File::open_mode_ = open_mode;
		if (open_mode & ios::in) {
			// load the file if we're supposed to.
			QFile file(filename.c_str());
			if (file.open(QIODevice::ReadOnly))
			{
				// doesn't really matter if it succeeds, since an empty dom is still valid in here.
				data_->setContent(&file);
				file.close();
			}
		}
	}

	FragmentXMLFile::~FragmentXMLFile() {
		delete data_;
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

	Atom* FragmentXMLFile::atomFromAtomDomForVariant(QDomElement& atomElement, const String& variant) {
		Atom* theAtom = new Atom;

		theAtom->setName(String(atomElement.attribute("id")));
		// TODO: preserve aliases
		String element = String(atomElement.elementsByTagName("element").at(0).toElement().text());
		theAtom->setElement(PTE_::getElement(element));
		QDomNodeList positions = atomElement.elementsByTagName("position");
		for (int i = 0; i < positions.count(); i++) {
			QDomElement pos = positions.item(i).toElement();
			if (pos.attribute("variant") == variant.c_str()) {
				Vector3 posVector;
				posVector.x = pos.attribute("x").toFloat();
				posVector.y = pos.attribute("y").toFloat();
				posVector.z = pos.attribute("z").toFloat();
				theAtom->setPosition(posVector);
			}
		}

		return theAtom;
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
		float bondOrder = asElement.attribute("order").toFloat();
		// bond order: x.5 = x, aromatic; x = x
		bool isAromatic = (bondOrder - floorf(bondOrder)) > 0.4f;
		// FIXME: this matches up with BALL::Bond::Order, but is not really future-proof.
		theBond->setOrder((int)floorf(bondOrder));
		if (isAromatic) theBond->setOrder(Bond::ORDER__AROMATIC);

		return true;
	}

	Molecule* FragmentXMLFile::getMoleculeForVariant(const String &variantName)
	{
		Molecule* theMolecule = new Molecule;
		QDomNodeList inList = data_->elementsByTagName("in");
		// first pass: create all the atoms.
		for (Position i = 0; i < inList.length(); i++) {
			QDomElement inElement = inList.item(i).toElement();
			if (inElement.attribute("variant") == variantName.c_str()) {
				// this thing occurs in this variant.
				// let's look at what it is. we're at $thing/occurs/as/in
				// so the information we want is in our grand-grand-parent.
				QDomElement thing = inElement.parentNode().parentNode().parentNode().toElement();
				if (thing.nodeName() == "atom") {
					Atom* atom = atomFromAtomDomForVariant(thing, variantName);
					theMolecule->append(*atom);
				}
			}
		}
		// second pass: create all bonds.
		for (Position i = 0; i < inList.length(); i++) {
			QDomElement inElement = inList.item(i).toElement();
			if (inElement.attribute("variant") == variantName.c_str()) {
				QDomElement thing = inElement.parentNode().parentNode().parentNode().toElement();
				if (thing.nodeName() == "bond") {
					bondFromInDomForMolecule(inElement, theMolecule);
				}
			}
		}

		return theMolecule;
	}

	bool FragmentXMLFile::write(const Molecule &molecule)
	{
		// idea: merge molecule into data, then write everything
		// since this isn't necessary a linear, linebased format, but interleaved.
		return false;
	}

	Molecule* FragmentXMLFile::read()
	{
		std::vector<String> variants = getVariantNames();
		if (currentVariant_ > variants.size()) return 0;

		return getMoleculeForVariant(variants[currentVariant_++]);
	}

	Molecule* FragmentXMLFile::read(const String &variantName)
	{
		return getMoleculeForVariant(variantName);
	}
}
