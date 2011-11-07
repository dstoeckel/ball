// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//

#ifndef BALL_STRUCTURE_FRAGMENTDB_FRAGMENTQUERY_H
#define BALL_STRUCTURE_FRAGMENTDB_FRAGMENTQUERY_H

#include<BALL/DATATYPE/string.h>
#include<BALL/KERNEL/fragment.h>

#include<boost/smart_ptr.hpp>
#include<boost/any.hpp>
#include<set>

namespace BALL
{

	/** Abstract Class for querying Molecular (Fragment) Databases, such as the FragmentDB.
	 *  A FragmentQuery holds some deatiled information about the Query, such as features
	 *  to query on. These can be obtained via the #selectsOn() and #getSelectorDetail() mechanics.
	 */
	class FragmentQuery
	{
	public:
		/** Enum containing the Query Types known to the library.
		 *  If a FragmentQuery #selectsOn() one of these, it means that #getSelectorDetail()
		 *  will return a pointer to the corresponding detail Class.
		 */
		enum QuerySelector {
			QueryFragmentName = 0,     //< will return a pointer to a #NameFragmentQuery
			QueryFragmentProperties = 1,
//			QueryFragmentVariants, //< will return a pointer to a #VariantFragmentQuery
			QueryConnectionTable = 2,
			QueryNameMap = 3,
			MAX_QUERYSELECTORS
		};

		typedef std::set<boost::shared_ptr<Residue> > ResultSet;

		FragmentQuery();
		FragmentQuery(unsigned int maximum_results) : max_results_(maximum_results) {}

		/** Get a human-readable representation of the Query. */
		virtual String toString() = 0;

		/** Retrieve the Results of the Query -- if any -- after it has been resolved. */
		virtual const ResultSet& getResults() { return results_; }

		/** Add a Fragment to the result set of the Query. */
		virtual void addResult(const boost::shared_ptr<Residue> &p) { results_.insert(p); }

		/** Check whether this Query has detailed information about a particular Fragment Feature. 
		 *  A 'true' here implies that #getSelectorDeatil() will return the proper child type
		 *  according to the #QuerySelector.
		 */
		virtual bool selectsOn(QuerySelector) = 0;
		
		/** Get a detailed information record out of this query.
		 *  If the query doesn't #selectsOn() the #QuerySelector, this will return NULL.
		 *  Think SQLs "WHERE" clauses.
		 */
		virtual boost::any getSelectorDetail(QuerySelector) = 0;

		/** The maximum number of results to query for.
		 *  Think SQLs "LIMIT" clause.
		 */
		unsigned int getMaxResults() { return max_results_; }
	protected:
		ResultSet results_;
		unsigned int max_results_;
	};

}

#endif // FRAGMENTQUERY_H
