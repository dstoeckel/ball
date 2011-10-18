// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//

#include <BALL/STRUCTURE/atomBijection.h>

#include <BALL/STRUCTURE/geometricProperties.h>
#include <BALL/KERNEL/PTE.h>
#include <BALL/KERNEL/extractors.h>
#include <BALL/KERNEL/residue.h>
#include <BALL/DATATYPE/hashGrid.h>

using namespace std;

namespace BALL
{

	AtomBijection::AtomBijection(AtomContainer& A, AtomContainer& B, bool limit_to_selection)
		: std::vector<std::pair<Atom*, Atom*> >()
	{
		assignByName(A, B, limit_to_selection);
	}

	/* Calculate the root mean squared deviation */
	double AtomBijection::calculateRMSD() const
	{
		double sum_of_squares = 0.0;

		if (!empty())
		{
			for (Size i = 0; i < size(); ++i)
			{
				Vector3& r(operator [] (i).first->getPosition());
				sum_of_squares += r.getSquareDistance(operator [] (i).second->getPosition());
			}

			// calculate mean square deviation
			sum_of_squares = sqrt(sum_of_squares / (double)size());
		}

		// return RMSD
		return sum_of_squares;
	}
} // namespace BALL
