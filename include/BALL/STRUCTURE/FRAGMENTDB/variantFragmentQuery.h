// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//

#ifndef BALL_STRUCTURE_FRAGMENTDB_VARIANTFRAGMENTQUERY_H
#define BALL_STRUCTURE_FRAGMENTDB_VARIANTFRAGMENTQUERY_H

#include<BALL/STRUCTURE/FRAGMENTDB/fragmentQuery.h>

namespace BALL {

	/** A FragmentQuery for enabling or limiting the number of variants returend in a query.
	 *  This Query only makes sense in conjunction with other selectors.
	 *  If another Query should return multiple Variants, it should selectOn(VariantQuery).
	 *  (And all that is entailed by that, @see FragmentQuery.)
	 */
	class VariantFragmentQuery : public FragmentQuery
	{
		public:
			VariantFragmentQuery();
			VariantFragmentQuery(unsigned int max_number_of_variants);

			bool selectsOn(QuerySelector);
			void* getSelectorDetail(QuerySelector);

			/** Maximum number of Variants to find.
			 *  Zero (0) means unlimited.
			 */
			unsigned int getMaxNumberOfVariants();
		private:
			unsigned int max_variants_;
	};
}

#endif // VARIANTFRAGMENTQUERY_H
