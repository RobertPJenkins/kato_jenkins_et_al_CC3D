#ifndef SPHEROIDADHESIONECMCONCENTRATION_EXPORT_H 
#define SPHEROIDADHESIONECMCONCENTRATION_EXPORT_H

    #if defined(_WIN32)
      #ifdef SpheroidAdhesionECMConcentration_EXPORTS
          #define SPHEROIDADHESIONECMCONCENTRATION_EXPORT __declspec(dllexport)
          #define SPHEROIDADHESIONECMCONCENTRATION_EXPIMP_TEMPLATE
      #else
          #define SPHEROIDADHESIONECMCONCENTRATION_EXPORT __declspec(dllimport)
          #define SPHEROIDADHESIONECMCONCENTRATION_EXPIMP_TEMPLATE extern
      #endif
    #else
         #define SPHEROIDADHESIONECMCONCENTRATION_EXPORT
         #define SPHEROIDADHESIONECMCONCENTRATION_EXPIMP_TEMPLATE
    #endif

#endif