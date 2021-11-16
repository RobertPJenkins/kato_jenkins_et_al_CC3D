#ifndef SPHEROIDCENTROID_EXPORT_H
#define SPHEROIDCENTROID_EXPORT_H

    #if defined(_WIN32)
      #ifdef SpheroidCentroidShared_EXPORTS
          #define SPHEROIDCENTROID_EXPORT __declspec(dllexport)
          #define SPHEROIDCENTROID_EXPIMP_TEMPLATE
      #else
          #define SPHEROIDCENTROID_EXPORT __declspec(dllimport)
          #define SPHEROIDCENTROID_EXPIMP_TEMPLATE extern
      #endif
    #else
         #define SPHEROIDCENTROID_EXPORT
         #define SPHEROIDCENTROID_EXPIMP_TEMPLATE
    #endif

#endif
