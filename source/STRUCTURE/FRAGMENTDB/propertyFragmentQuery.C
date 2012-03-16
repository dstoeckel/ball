#include <BALL/STRUCTURE/FRAGMENTDB/propertyFragmentQuery.h>

namespace BALL {

	PropertyFragmentQuery::PropertyFragmentQuery(const PropertyManager &properties, unsigned int query_limit) 
		:	FragmentQuery(query_limit),
			pm_(properties)
	{
	}

	String PropertyFragmentQuery::toString() {
		String result;
		result += "Property Query: Unnamed: ";
		result += pm_.getBitVector().getUnsignedLong();
		result += " Named: (";
		NamedPropertyIterator it = pm_.beginNamedProperty();
		for (;it != pm_.endNamedProperty(); ++it)
		{
			result += it->getName();
		}
		result += ")";
	}

	const PropertyManager& PropertyFragmentQuery::getPropertyManager()
	{
		return pm_;
	}

	bool PropertyFragmentQuery::selectsOn(QuerySelector q)
	{
		return (q == QueryFragmentProperties);
	}

	void* PropertyFragmentQuery::getSelectorDetail(QuerySelector q)
	{
		if (selectsOn(q))
		{
			return this;
		}
		else
		{
			return NULL;
		}
	}
}
