// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//

#ifndef BALL_STRUCTURE_FRAGMENTDB_PROPERTYFRAGMENTQUERY_H
#define BALL_STRUCTURE_FRAGMENTDB_PROPERTYFRAGMENTQUERY_H

#include<BALL/STRUCTURE/FRAGMENTDB/fragmentQuery.h>
#include<BALL/CONCEPT/property.h>

namespace BALL
{

class PropertyFragmentQuery : public FragmentQuery
{
	public:
		PropertyFragmentQuery(const PropertyManager& properties , unsigned int query_limit = 1);
		
		bool selectsOn(QuerySelector);
		void* getSelectorDetail(QuerySelector);
		String toString();

		const PropertyManager& getPropertyManager();
	private:
		PropertyManager pm_;
};

}

#endif // PROPERTYFRAGMENTQUERY_H
