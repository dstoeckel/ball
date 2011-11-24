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

#ifndef BALL_SYSTEM_PATH_H
#include <BALL/SYSTEM/path.h>
#endif

#ifndef BALL_DATATYPE_STRING_H
#include <BALL/DATATYPE/string.h>
#endif

#ifndef BALL_KERNEL_RESIDUE_H
#include <BALL/KERNEL/residue.h>
#endif

#ifndef BALL_COMMON_LIMITS_H
#include <BALL/COMMON/limits.h>
#endif

#ifndef BALL_MATHS_VECTOR3_H
#include <BALL/MATHS/vector3.h>
#endif

#ifndef BALL_STRUCTURE_FRAGMENTDB_CONNECTION_H
#include <BALL/STRUCTURE/FRAGMENTDB/connection.h>
#endif

#include <set>
#include <QXmlSchema>
#include <QXmlSchemaValidator>


namespace BALL
{


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
		for (int i = 0; i < nodes.count(); i++)
		{
			QDomElement node = nodes.item(i).toElement();
			if (node.hasAttribute("id")) {
				variants.push_back(String(node.attribute("id")));
			}
		}
		return variants;
	}

	bool FragmentXMLFile::write(const Molecule &molecule)
	{
		// idea: merge molecule into data, then write everything
		// since this isn't necessarily a linear, linebased format, but interleaved.
		Log.warn() << "Writing single Molecules is not supported as of now, please write a complete System." << std::endl;
		throw (File::CannotWrite(__FILE__, __LINE__, File::name_));
		return false;
	}

	Molecule* FragmentXMLFile::read()
	{
		std::vector<String> variants = getVariantNames();
		if (currentVariant_ >= variants.size()) return 0;
		return moleculeForVariant(variants[currentVariant_++]);
	}

	GenericMolFile& FragmentXMLFile::operator >>(System& sys)
	{
		GenericMolFile::operator >>(sys);
		sys.setProperty("SKIP_NORMALIZATION", true);
		return *this;
	}

	void FragmentXMLFile::parse() {
		// phase 1 - prepare all variants.
		QDomNodeList variants = data_->elementsByTagName("variant");
		for (int i = 0; i < variants.count(); ++i)
		{
			QDomElement node = variants.item(i).toElement();
			if (node.hasAttribute("id"))
			{
				String variantID(node.attribute("id"));
				Residue* res = new Residue();
				res->setName(variantID);
				residues_.push_back(res);
				residue_by_name_.insert(variantID,residues_.size()-1);
				Log.info() << "Created residue " << variantID << std::endl;
				parseNames(node, *res);
				QDomNodeList props = node.elementsByTagName("property");
				for (int j = 0; j < props.count(); ++j) 
				{
					QDomElement propertyEl = props.at(i).toElement();
					parseProperty(propertyEl, *res);
				}
			}
		}

		// phase 2 -- parse atoms
		QDomNodeList atoms = data_->elementsByTagName("atom");
		for (int i = 0; i < atoms.count(); ++i)
		{
			QDomElement atom = atoms.item(i).toElement();
			if (atom.hasAttribute("id"))
			{
				String atomID(atom.attribute("id"));
				// TODO: only occurs/in
				QDomNodeList occurences = atom.elementsByTagName("in");
				for (int j = 0; j < occurences.count(); ++j)
				{
					QDomElement occurs_in = occurences.item(j).toElement();
					if (occurs_in.hasAttribute("variant"))
					{
						String occurence_variant(occurs_in.attribute("variant"));
						if (residue_by_name_.has(occurence_variant)) {
							if (!residues_[residue_by_name_[occurence_variant]]->getAtom(atomID))
							{
								parseAtom(atom, residues_[residue_by_name_[occurence_variant]]);
							}
						}
					}
				}
			}
		}

		// phase 3 -- parse bonds
		QDomNodeList bonds = data_->elementsByTagName("bond");
		for (int i = 0; i < bonds.count(); ++i)
		{
			QDomElement bond = bonds.item(i).toElement();
			QDomNodeList occurences = bond.elementsByTagName("in");
			for (int j = 0; j < occurences.count(); ++j)
			{
				QDomElement occurs_in = occurences.item(j).toElement();
				if (occurs_in.hasAttribute("variant"))
				{
					String occurence_variant(occurs_in.attribute("variant"));
					if (residue_by_name_.has(occurence_variant))
					{
						parseBond(bond, residues_[residue_by_name_[occurence_variant]]);
					}
				}
			}
		}
	}

	void FragmentXMLFile::parseNames(QDomElement& name_container, PropertyManager& props)
	{
		QDomNodeList names = name_container.elementsByTagName("name");
		for (int n = 0; n < names.count(); ++n)
		{
			QDomElement name = names.at(n).toElement();
			if (name.hasAttribute("convention"))
			{
				String conventionName(name.text());
				String convention(name.attribute("convention"));
				registerNameForConventions(conventionName,convention,props);
			}
		}
	}

	void FragmentXMLFile::registerNameForConventions(String& name, String& conventions, PropertyManager& object)
	{
		// split the attribute since it may actually be multiple conventions.
		std::vector<String> all_conventions;
		conventions.split(all_conventions);
		for (std::vector<String>::iterator cv = all_conventions.begin();
		     cv != all_conventions.end(); ++cv)
		{
			// save every name as a property. (we'll reconstruct the
			// compact representation with collectNameForConventions)
			object.setProperty("NAME:"+*cv,name);
		}
	}

	void FragmentXMLFile::collectNameForConventions(PropertyManager& object, StringHashMap<String>& names)
	{
		NamedPropertyIterator prop_it = object.beginNamedProperty();
		for (;prop_it != object.endNamedProperty(); ++prop_it) {
			if ((prop_it->getType() == NamedProperty::STRING) && 
			    (prop_it->getName().substr(0,4) == "NAME:"))
			{
				String convention(prop_it->getName());
				convention = convention.getField(1,":");
				const String& name(prop_it->getString());
				if (names.has(name))
				{
					names[name] += " ";
					names[name] += convention;
				}
				else
				{
					names.insert(name, convention);
				}
			}
		}
	}

	void FragmentXMLFile::parseAtom(QDomElement &domAtom, Residue *residue)
	{
		// first we create the atom
		QDomElement domElement = domAtom.elementsByTagName("element").at(0).toElement();
		String atom_id(domAtom.attribute("id"));
		if (domElement.isNull())
		{
			Log.warn() << "No element defined for atom " << atom_id << ", ignoring." << std::endl;
			return;
		}
		Atom* newAtom = new Atom();
		newAtom->setElement(PTE[String(domElement.text())]);
		newAtom->setName(atom_id);

		const char* residuename = residue->getName().c_str();

		// then figure out position...
		QDomNodeList occurences = domAtom.elementsByTagName("in");
		for (int n = 0; n < occurences.count(); ++n)
		{
			QDomElement domIn = occurences.at(n).toElement();
			if (domIn.hasAttribute("variant") && 
			    (domIn.attribute("variant") == residuename))
			{
				// so this atom occurs in this variant.
				// go up in the tag hierarchy to figure out the details.
				QDomElement domParent = domIn.parentNode().toElement();
				// it may either be the position.
				if (domParent.tagName() == "position") 
				{
					newAtom->setPosition(Vector3(
					domParent.attribute("x").toFloat(),
					domParent.attribute("y").toFloat(),
					domParent.attribute("z").toFloat()
					));
				}
				// or a plain occurence.
				else if (domParent.tagName() == "occurs")
				{
					// in which case we don't even need to do something.
				}
				else
				{
					Log.warn() << "Unexpected tag " << String(domParent.tagName()) 
					           << " on line " << domParent.lineNumber() << std::endl;
				}
			}
		}

		// parse the properties
		QDomNodeList props = domAtom.elementsByTagName("property");
		for (int j = 0; j < props.count(); ++j) 
		{
			QDomElement propertyEl = props.at(j).toElement();
			parseProperty(propertyEl, *newAtom);
		}

		// the names
		parseNames(domAtom, *newAtom);
		// and connections
		parseConnections(domAtom, newAtom);

		residue->appendChild(*newAtom);
	}

	void FragmentXMLFile::parseConnections(QDomElement &domAtom, Atom* atom)
	{
		QDomNodeList connections = domAtom.elementsByTagName("connection");
		for (int i = 0; i < connections.count(); ++i)
		{
			QDomElement connection = connections.at(i).toElement();
#define ASSERT_ATTRIBUTE(x) if (!connection.hasAttribute(x)) { \
				Log.warn() << "Connection for atom " << String(domAtom.attribute("id")) \
				           << " is missing required attribute " << x \
				           << " on line " << connection.lineNumber() \
				           << ". Ignoring connection." << std::endl; \
				continue; \
			}
			ASSERT_ATTRIBUTE("type")
			ASSERT_ATTRIBUTE("partnertype")
			ASSERT_ATTRIBUTE("order")
			ASSERT_ATTRIBUTE("distance")
			ASSERT_ATTRIBUTE("deviation")
#undef ASSERT_ATTRIBUTE
			Connection conn;
			conn.atom = atom;
			conn.type_name = String(connection.attribute("type"));
			conn.connect_to = String(connection.attribute("partnertype"));
			conn.dist = connection.attribute("distance").toFloat();
			conn.delta = connection.attribute("deviation").toFloat();
			conn.order = parseBondOrder(connection.attribute("order"));

			boost::shared_ptr<PersistentObject> connptr(new Connection(conn));
			String connID = "CONNECTION:";
			connID += String(i);
			atom->setProperty(NamedProperty(connID, connptr));
		}
	}

	void FragmentXMLFile::parseBond(QDomElement &domBond, Residue *residue)
	{
		QDomNodeList partners = domBond.elementsByTagName("partner");
		if (partners.count() != 2) {
			Log.warn() << "Disregarding Bond without exactly 2 partner atoms on line"
			           << domBond.lineNumber() << std::endl;
			return;
		}
		String partner_a(partners.at(0).toElement().attribute("atom"));
		String partner_b(partners.at(1).toElement().attribute("atom"));
		Atom* atom_a = residue->getAtom(partner_a);
		if (!atom_a)
		{
			Log.warn() << "Can't find atom "<< partner_a<< " in variant "
			           << residue->getName() << ", disregarding bond on line "
			           << domBond.lineNumber() << std::endl;
			return;
		}
		Atom* atom_b = residue->getAtom(partner_b);
		if (!atom_b)
		{
			Log.warn() << "Can't find atom "<< partner_b << " in variant "
			           << residue->getName() << ", disregarding bond on line "
			           << domBond.lineNumber() << std::endl;
			return;
		}
		Bond::Order order = Bond::ORDER__UNKNOWN;
		QDomNodeList inVariant = domBond.elementsByTagName("in");
		for (int i = 0; i < inVariant.count(); ++i)
		{
			QDomElement domIn = inVariant.at(i).toElement();
			if (domIn.attribute("variant") == residue->getName().c_str())
			{
				QDomElement domOrder = domIn.parentNode().toElement();
				if (domOrder.nodeName() == "as" && domOrder.hasAttribute("order")) {
					if (order != Bond::ORDER__UNKNOWN) {
						Log.warn() << "Bond order redefinition for " << partner_a << "-" << partner_b
						           << " in variant " << residue->getName() << " ignored "
						           << "on line " << domOrder.lineNumber();
					}
					else
					{
						order = parseBondOrder(domOrder.attribute("order"));
						if (order == Bond::ORDER__UNKNOWN)
						{
							Log.warn() << "Unknown Bond order '"<< String(domOrder.attribute("order")) 
							           << "' on line " << domOrder.lineNumber() << std::endl;
						}
					}
				}
			}
		}
		// the constructor of bond takes care of the registration.
		Bond* newBond = new Bond("",*atom_a,*atom_b,order);
	}

	Bond::Order FragmentXMLFile::parseBondOrder(const QString &theOrder)
	{
		String order(theOrder);
		if (order == "aromatic")
		{
			return Bond::ORDER__AROMATIC;
		}
		else if (order == "1")
		{
			return Bond::ORDER__SINGLE;
		}
		else if (order == "2")
		{
			return Bond::ORDER__DOUBLE;
		}
		else if (order == "3")
		{
			return Bond::ORDER__TRIPLE;
		}
		else if (order == "4")
		{
			return Bond::ORDER__QUADRUPLE;
		}
		else
		{
			return Bond::ORDER__UNKNOWN;
		}
	}

	void FragmentXMLFile::parseProperty(QDomElement &domProperty, PropertyManager &object)
	{
		if (domProperty.hasAttribute("type") && (domProperty.attribute("type") != "UNNAMED"))
		{
			if (!domProperty.hasAttribute("name")) {
				Log.warn() << "Named property of type " 
				           << String(domProperty.attribute("type")) << " has no name at line "
				           << domProperty.lineNumber() << std::endl;
				return;
			}
			else
			{
				String type(domProperty.attribute("type"));
				String value(domProperty.text());
				String name(domProperty.attribute("name"));
				if (object.hasProperty(name))
				{
					Log.warn() << "Overwriting named property " << name << " with value "
					           << value << " (was: " << object.getProperty(name).toString() << ") " 
					           << "redefined on line " << domProperty.lineNumber() << std::endl;
				}
				// relying on the BALL::String conversions for consistency here,
				// even though QString provides some of them as well
				if (type == "bool") {
					bool val = value.toBool();
					object.setProperty(name, val);
				}
				else if (type == "int")
				{
					int val = value.toInt();
					object.setProperty(name, val);
				}
				else if (type == "unsigned")
				{
					unsigned int val = value.toUnsignedInt();
					object.setProperty(name, val);
				}
				else if (type == "float")
				{
					float val = value.toFloat();
					object.setProperty(name, val);
				}
				else if (type == "double")
				{
					double val = value.toDouble();
					object.setProperty(name, val);
				}
				else if (type == "string")
				{
					object.setProperty(name, value);
				}

			}
		}
		else
		{
			QString property = domProperty.text();
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
				Log.warn() << "Unreconized unnamed property: " << String(property)
				           << " on line " << domProperty.lineNumber() << std::endl;
			}
			else
			{
				object.setProperty(prop);
			}
		}
	}

	bool FragmentXMLFile::validate() {
		Path p;
		String schemaInPath = "file://" + p.find("fragments_xml/FragmentXML.xsd");
		QUrl schemaUrl(schemaInPath.c_str());

		QXmlSchema schema;
		schema.load(schemaUrl);

		if (schema.isValid()) {
			QFile file(File::name_.c_str());
			if (file.open(QIODevice::ReadOnly))
			{
				QXmlSchemaValidator validator(schema);
				return validator.validate(&file, QUrl::fromLocalFile(file.fileName()));
			}
			else
			{
				Log.warn() << "Cannot open File for verification of Schema: " << File::name_;
				return false;
			}
		}
		else
		{
			Log.warn() << "Cannot find Schema or Schema is invalid: " << schemaInPath;
			return false;
		}
	}

	Molecule* FragmentXMLFile::moleculeForVariant(const String & variant)
	{
		Molecule* mol = new Molecule();
		Residue* res = residueForVariant(variant);
		if (res) {
			mol->append(*res);
			return mol;
		}
		else
		{
			delete mol;
			return NULL;
		}
	}

	Residue* FragmentXMLFile::residueForVariant(const String & variant)
	{
		if (residue_by_name_.has(variant)) {
			return residues_[residue_by_name_[variant]];
		}
		else
		{
			return NULL;
		}
	}

} // namespace BALL
