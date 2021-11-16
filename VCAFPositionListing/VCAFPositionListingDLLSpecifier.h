#ifndef VCAFPOSITIONLISTING_EXPORT_H
#define VCAFPOSITIONLISTING_EXPORT_H

    #if defined(_WIN32)
      #ifdef VCAFPositionListing_EXPORTS
          #define VCAFPOSITIONLISTING_EXPORT __declspec(dllexport)
          #define VCAFPOSITIONLISTING_EXPIMP_TEMPLATE
      #else
          #define VCAFPOSITIONLISTING_EXPORT __declspec(dllimport)
          #define VCAFPOSITIONLISTING_EXPIMP_TEMPLATE extern
      #endif
    #else
         #define VCAFPOSITIONLISTING_EXPORT
         #define VCAFPOSITIONLISTING_EXPIMP_TEMPLATE
    #endif

#endif
