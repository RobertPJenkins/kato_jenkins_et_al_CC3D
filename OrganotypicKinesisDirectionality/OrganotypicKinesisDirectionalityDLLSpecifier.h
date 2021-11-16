#ifndef ORGANOTYPICKINESISDIRECTIONALITY_EXPORT_H
#define ORGANOTYPICKINESISDIRECTIONALITY_EXPORT_H

    #if defined(_WIN32)
      #ifdef OrganotypicKinesisDirectionalityShared_EXPORTS
          #define ORGANOTYPICKINESISDIRECTIONALITY_EXPORT __declspec(dllexport)
          #define ORGANOTYPICKINESISDIRECTIONALITY_EXPIMP_TEMPLATE
      #else
          #define ORGANOTYPICKINESISDIRECTIONALITY_EXPORT __declspec(dllimport)
          #define ORGANOTYPICKINESISDIRECTIONALITY_EXPIMP_TEMPLATE extern
      #endif
    #else
         #define ORGANOTYPICKINESISDIRECTIONALITY_EXPORT
         #define ORGANOTYPICKINESISDIRECTIONALITY_EXPIMP_TEMPLATE
    #endif

#endif
