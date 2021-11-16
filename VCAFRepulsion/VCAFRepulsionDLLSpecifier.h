#ifndef VCAFREPULSION_EXPORT_H
#define VCAFREPULSION_EXPORT_H

    #if defined(_WIN32)
      #ifdef VCAFRepulsionShared_EXPORTS
          #define VCAFREPULSION_EXPORT __declspec(dllexport)
          #define VCAFREPULSION_EXPIMP_TEMPLATE
      #else
          #define VCAFREPULSION_EXPORT __declspec(dllimport)
          #define VCAFREPULSION_EXPIMP_TEMPLATE extern
      #endif
    #else
         #define VCAFREPULSION_EXPORT
         #define VCAFREPULSION_EXPIMP_TEMPLATE
    #endif

#endif
