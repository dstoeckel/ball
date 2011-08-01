
#include <BALL/STRUCTURE/FRAGMENTDB/nameMapQuery.h>

namespace BALL {

	NameMapQuery::NameMapQuery(const String &part_of_name)
		: FragmentQuery(0),
		  name_part_(part_of_name)
	{
	}

	String NameMapQuery::toString()
	{
		return "Blah"; // FIXME
	}

	void* NameMapQuery::getSelectorDetail(QuerySelector q)
	{
		if (selectsOn(q)) 
		{
			return this;
		}
	}

	bool NameMapQuery::selectsOn(QuerySelector q)
	{
		return (q == QueryNameMap);
	}

	const StringHashMap<NameMapQuery::NameMap*>& NameMapQuery::getMaps()
	{
		return maps_;
	}

	void NameMapQuery::addMap(const String & translate, NameMapQuery::NameMap * map)
	{
		maps_.insert(translate, map);
	}

	const String& NameMapQuery::getMapName()
	{
		return name_part_;
	}
}
