// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:

#ifndef BALL_STRUCTURE_FORCEFIELDEVALUATION_H
#define BALL_STRUCTURE_FORCEFIELDEVALUATION_H

#ifndef BALL_STRUCTURE_DOCKING_ENERGETICEVALUATION_H
# include <BALL/STRUCTURE/DOCKING/energeticEvaluation.h>
#endif

#ifndef BALL_MOLMEC_COMMON_FORCEFIELD_H
# include <BALL/MOLMEC/COMMON/forceField.h>
#endif

namespace BALL
{
		/** Base class for energetic evaluators of docking results using
		    a force field as scoring function.
				\ingroup Docking
		 */
		class BALL_EXPORT ForceFieldEvaluation 
			: public EnergeticEvaluation
		{
			public:

				/// Default constructor.
				ForceFieldEvaluation()
					;
				
				/// 
				ForceFieldEvaluation(ForceField& ff)
					;
				
				/// 
				virtual ~ForceFieldEvaluation()
					;
				
				/** Operations
				*/
				void setForceField(ForceField& ff)
					;

				/// 
				void setOptions(const Options& options)
					;

				/// 
				ForceField& getForceField()
					;

				/// 
				const ForceField& getForceField() const
					;

				/// 
				Options& getOptions()
					;

				/// 
				const Options& getOptions() const
					;
				
				/// 
				virtual std::vector<ConformationSet::Conformation> operator () (ConformationSet& conformations)
					throw(Exception::TooManyErrors);

			protected:

				ForceField* ff_;
				Options     options_;
				bool 				delete_force_field_;
		};
}
#endif