// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//

#ifndef BALL_STRUCTURE_FRAGMENTDB_FRAGMENTSTORAGE_H
#define BALL_STRUCTURE_FRAGMENTDB_FRAGMENTSTORAGE_H

#include<BALL/DATATYPE/string.h>
#include<bitset>

namespace BALL
{
	class FragmentQuery;

	class FragmentStorage
	{
	public:
		/** Retrieve the Storages version tag.
		 *  The version tag should change if the underlying data changes.
		 *  Other layers of the database hierarchy may use it as an indicator
		 *  on when to update indices or other support structures.
		 */
		virtual String getVersionTag() = 0;

		/** Submit a query to the Storage.
		 *  If the query (or part of it) can be handled by the Storage,
		 *  it is updated with the results as pertains to this Storage.
		 *  In that case, "true" is returned. If the Storage can't do
		 *  anything to that query, it returns "false".
		 */
		virtual bool query(FragmentQuery&) = 0;
	};

}

#endif // FRAGMENTSTORAGE_H
