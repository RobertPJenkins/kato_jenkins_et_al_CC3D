#ifndef ECMPENETRATION_EXPORT_H
#define ECMPENETRATION_EXPORT_H

    #if defined(_WIN32)
      #ifdef ECMPenetrationShared_EXPORTS
          #define ECMPENETRATION_EXPORT __declspec(dllexport)
          #define ECMPENETRATION_EXPIMP_TEMPLATE
      #else
          #define ECMPENETRATION_EXPORT __declspec(dllimport)
          #define ECMPENETRATION_EXPIMP_TEMPLATE extern
      #endif
    #else
         #define ECMPENETRATION_EXPORT
         #define ECMPENETRATION_EXPIMP_TEMPLATE
    #endif

#endif
