// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//

#include<BALL/STRUCTURE/FRAGMENTDB/nameFragmentQuery.h>

namespace BALL
{
	NameFragmentQuery::NameFragmentQuery(const String& fragment_name, const String& naming_standard, unsigned int max_results)
		: FragmentQuery(max_results),
			fragment_name_(fragment_name),
			naming_standard_(naming_standard)
	{
	}

	bool NameFragmentQuery::selectsOn(QuerySelector s)
	{
		return (s == QueryFragmentName);
	}
	
	boost::any NameFragmentQuery::getSelectorDetail(QuerySelector s) {
		switch (s) {
			case QueryFragmentName:
				return boost::any(this);
			default:
				return NULL;
		}
	}
	
	String NameFragmentQuery::getFragmentName()
	{
		return fragment_name_;
	}
	
	String NameFragmentQuery::getNamingStandard()
	{
		return naming_standard_;
	}
	
	String NameFragmentQuery::toString()
	{
		String name;
		name += "NameFragmentQuery: ";
		name += naming_standard_;
		name += "/";
		name += fragment_name_;
		name += "; ";
		name += ((getMaxResults() == 0)? "unlimited results" : "max "+getMaxResults());
		return name;
	}
}
