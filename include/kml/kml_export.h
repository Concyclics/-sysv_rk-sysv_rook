
#ifndef KML_EXPORT_H
#define KML_EXPORT_H

#ifdef KML_STATIC_DEFINE
#  define KML_EXPORT
#  define KML_NO_EXPORT
#else
#  ifndef KML_EXPORT
#    ifdef Service_so_EXPORTS
        /* We are building this library */
#      define KML_EXPORT __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define KML_EXPORT __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef KML_NO_EXPORT
#    define KML_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef KML_DEPRECATED
#  define KML_DEPRECATED 
#endif

#ifndef KML_DEPRECATED_EXPORT
#  define KML_DEPRECATED_EXPORT KML_EXPORT KML_DEPRECATED
#endif

#ifndef KML_DEPRECATED_NO_EXPORT
#  define KML_DEPRECATED_NO_EXPORT KML_NO_EXPORT KML_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef KML_NO_DEPRECATED
#    define KML_NO_DEPRECATED
#  endif
#endif

#endif /* KML_EXPORT_H */
