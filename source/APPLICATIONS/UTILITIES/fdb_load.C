// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
//

#include<iostream>

#include <BALL/FORMAT/fragmentXMLFile.h>
#include <BALL/KERNEL/forEach.h>
#include <BALL/DATATYPE/string.h>
#include <BALL/KERNEL/system.h>
#include <BALL/KERNEL/molecule.h>
#include <BALL/KERNEL/atom.h>
#include <BALL/KERNEL/bond.h>
#include <BALL/KERNEL/bondIterator.h>

using namespace BALL;

int main(int argc, char** argv)
{
	if (argc < 2) {
		std::cerr << "Usage :" << argv[0] << " <fragmentXML input>" << std::endl;
		return 1;
	}

	String infilename(argv[1]);
	FragmentXMLFile fragments(infilename);
	System testsystem;

	fragments >> testsystem;

	std::cout << "Contents:" << std::endl << "------" << std::endl;

	MoleculeIterator mol;
	BALL_FOREACH_MOLECULE(testsystem, mol)
	{
		std::cout << "Fragment:" << mol->getName() << std::endl;
		AtomIterator atom;
		BALL_FOREACH_ATOM(*mol, atom)
		{
			std::cout << "`- Atom:" << atom->getName() << std::endl;
			AtomBondIterator bond;
			BALL_FOREACH_ATOM_BOND(*atom, bond)
			{
				std::cout << "   `- Bond to:" << bond->getPartner(*atom)->getName() << std::endl;
			}
		}
	}
}
