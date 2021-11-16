
#ifndef MITOSISPOSTINITIALISE_EXPORT_H 
#define MITOSISPOSTINITIALISE_EXPORT_H

    #if defined(_WIN32)
      #ifdef MitosisPostInitialise_EXPORTS
          #define MITOSISPOSTINITIALISE_EXPORT __declspec(dllexport)
          #define MITOSISPOSTINITIALISE_EXPIMP_TEMPLATE
      #else
          #define MITOSISPOSTINITIALISE_EXPORT __declspec(dllimport)
          #define MITOSISPOSTINITIALISE_EXPIMP_TEMPLATE extern
      #endif
    #else
         #define MITOSISPOSTINITIALISE_EXPORT
         #define MITOSISPOSTINITIALISE_EXPIMP_TEMPLATE
    #endif

#endif
