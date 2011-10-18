
#include <BALL/STRUCTURE/atomBiMap.h>
#include <BALL/MATHS/vector3.h>
#include <BALL/KERNEL/atom.h>
#include <BALL/KERNEL/atomContainer.h>
#include <BALL/KERNEL/extractors.h>
#include <BALL/KERNEL/residue.h>
#include <BALL/DATATYPE/hashGrid.h>

namespace BALL {

	/* Calculate the root mean squared deviation */
	double AtomBiMap::calculateRMSD() const
	{
		double sum_of_squares = 0.0;
		
		if (!empty())
		{
			left_const_iterator it = left.begin();

			for(; it != left.end(); ++it) {
				Vector3& r(it->first->getPosition());
				sum_of_squares += r.getSquareDistance(it->second->getPosition());
			}

			// calculate mean square deviation
			sum_of_squares = sqrt(sum_of_squares / (double)size());
		}

		// return RMSD
		return sum_of_squares;
	}

	Atom* AtomBiMap::leftFromRight(Atom* atom) {
		right_const_iterator it = right.find(atom);
		if (it == right.end())
		{
			return 0;
		}
		else
		{
			return it->second;
		}
	}

	Atom* AtomBiMap::rightFromLeft(Atom* atom) {
		left_const_iterator it = left.find(atom);
		if (it == left.end())
		{
			return 0;
		}
		else
		{
			return it->second;
		}
	}

}
