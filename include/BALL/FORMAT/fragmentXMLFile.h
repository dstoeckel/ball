// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//

#ifndef BALL_FORMAT_FRAGMENTXMLFILE_H
#define BALL_FORMAT_FRAGMENTXMLFILE_H

#ifndef BALL_FORMAT_GENERICMOLFILE_H
#   include <BALL/FORMAT/genericMolFile.h>
#endif

namespace BALL {
        /** Fragment XML file class.
            Support for reading fragmentXML files, one variant at a time.
            \ingroup StructureFormats
        */
        class BALL_EXPORT FragmentXMLFile
                  : public GenericMolFile
        {

        };
}

#endif // BALL_FORMAT_FRAGMENTXMLFILE_H
