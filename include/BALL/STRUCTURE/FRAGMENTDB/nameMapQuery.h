// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//

#ifndef BALL_STRUCTURE_FRAGMENTDB_NAMEMAPQUERY_H
#define BALL_STRUCTURE_FRAGMENTDB_NAMEMAPQUERY_H

#include<BALL/STRUCTURE/FRAGMENTDB/fragmentQuery.h>
#include<BALL/DATATYPE/stringHashMap.h>

namespace BALL
{

class NameMapQuery : public FragmentQuery
{
	public:
		typedef StringHashMap<String>		NameMap;

		NameMapQuery(const String& part_of_name);

		void* getSelectorDetail(QuerySelector);
		bool selectsOn(QuerySelector);
		void addMap(const String&, NameMap*);
		String toString();
		const StringHashMap<NameMap*>& getMaps();
		const String& getMapName();

	private:
		String name_part_;
		StringHashMap<NameMap*> maps_;
};

}

#endif // NAMEMAPQUERY_H
