// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//

#include <BALL/FORMAT/fragmentXMLFile.h>
#include <BALL/KERNEL/atom.h>
#include <BALL/KERNEL/bond.h>
#include <BALL/KERNEL/molecule.h>
#include <BALL/KERNEL/system.h>
#include <QtXml>

namespace BALL
{
	FragmentXMLFile::FragmentXMLFile()
		: GenericMolFile()
	{
		// always create the DomDocument container.
		data = new QDomDocument();
	}

	FragmentXMLFile::FragmentXMLFile(const String &filename, File::OpenMode open_mode)
	{
		FragmentXMLFile();
		if (open_mode == ios::in) {
			QFile file(filename.c_str());
			if (file.open(QIODevice::OpenMode::ReadOnly))
			{
				data->setContent(file);
			}
		}
	}

	FragmentXMLFile::~FragmentXMLFile() {
		delete data;
	}

	FragmentXMLFile::isOpen() {
		// this is used by several methods of the parent classes.
		// since here, we only operate on the in-memory representation
		// always return true
		return true;
	}

	FragmentXMLFile::write(const Molecule &molecule)
	{

	}

	Molecule* FragmentXMLFile::read()
	{

	}

	Molecule* FragmentXMLFile::read(const String &variantName)
	{

	}
}
