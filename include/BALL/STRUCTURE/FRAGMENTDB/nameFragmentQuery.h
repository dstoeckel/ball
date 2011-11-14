// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//

#ifndef BALL_STRUCTURE_FRAGMENTDB_NAMEFRAGMENTQUERY_H
#define BALL_STRUCTURE_FRAGMENTDB_NAMEFRAGMENTQUERY_H

#include<BALL/STRUCTURE/FRAGMENTDB/fragmentQuery.h>

namespace BALL
{

class NameFragmentQuery : public FragmentQuery
{
	public:
		NameFragmentQuery(const String& fragment_name, const String& naming_standard = "PDB", unsigned int query_limit = 1);
		
		bool selectsOn(QuerySelector);
		boost::any getSelectorDetail(QuerySelector);
		String toString();

		String getFragmentName();
		String getNamingStandard();
	private:
		String fragment_name_;
		String naming_standard_;
};

}

#endif // NAMEFRAGMENTQUERY_H
