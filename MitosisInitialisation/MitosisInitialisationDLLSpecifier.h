#ifndef MITOSISINITIALISATION_EXPORT_H 
#define MITOSISINITIALISATION_EXPORT_H

    #if defined(_WIN32)
      #ifdef MitosisInitialisation_EXPORTS
          #define MITOSISINITIALISATION_EXPORT __declspec(dllexport)
          #define MITOSISINITIALISATION_EXPIMP_TEMPLATE
      #else
          #define MITOSISINITIALISATION_EXPORT __declspec(dllimport)
          #define MITOSISINITIALISATION_EXPIMP_TEMPLATE extern
      #endif
    #else
         #define MITOSISINITIALISATION_EXPORT
         #define MITOSISINITIALISATION_EXPIMP_TEMPLATE
    #endif

#endif

