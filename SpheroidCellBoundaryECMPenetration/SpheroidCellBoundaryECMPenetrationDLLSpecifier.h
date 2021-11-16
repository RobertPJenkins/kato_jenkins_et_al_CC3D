#ifndef SPHEROIDCELLBOUNDARYECMPENETRATION_EXPORT_H 
#define SPHEROIDCELLBOUNDARYECMPENETRATION_EXPORT_H

    #if defined(_WIN32)
      #ifdef SpheroidCellBoundaryECMPenetration_EXPORTS
          #define SPHEROIDCELLBOUNDARYECMPENETRATION_EXPORT __declspec(dllexport)
          #define SPHEROIDCELLBOUNDARYECMPENETRATION_EXPIMP_TEMPLATE
      #else
          #define SPHEROIDCELLBOUNDARYECMPENETRATION_EXPORT __declspec(dllimport)
          #define SPHEROIDCELLBOUNDARYECMPENETRATION_EXPIMP_TEMPLATE extern
      #endif
    #else
         #define SPHEROIDCELLBOUNDARYECMPENETRATION_EXPORT
         #define SPHEROIDCELLBOUNDARYECMPENETRATION_EXPIMP_TEMPLATE
    #endif

#endif