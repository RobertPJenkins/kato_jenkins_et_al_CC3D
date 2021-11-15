#ifndef DIRECTIONALBIAS_EXPORT_H
#define DIRECTIONALBIAS_EXPORT_H

    #if defined(_WIN32)
      #ifdef DirectionalBiasShared_EXPORTS
          #define DIRECTIONALBIAS_EXPORT __declspec(dllexport)
          #define DIRECTIONALBIAS_EXPIMP_TEMPLATE
      #else
          #define DIRECTIONALBIAS_EXPORT __declspec(dllimport)
          #define DIRECTIONALBIAS_EXPIMP_TEMPLATE extern
      #endif
    #else
         #define DIRECTIONALBIAS_EXPORT
         #define DIRECTIONALBIAS_EXPIMP_TEMPLATE
    #endif

#endif
