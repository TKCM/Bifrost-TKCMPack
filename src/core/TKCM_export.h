#ifndef TKCM_EXPORT_H
#define TKCM_EXPORT_H

#ifndef AMINO_PARSING_HEADERS
    #if defined(_WIN32)
        #define TKCM_EXPORT __declspec(dllexport)
        #define TKCM_IMPORT __declspec(dllimport)
    #elif defined(__GNUC__)
        #define TKCM_EXPORT __attribute__ ((visibility("default")))
        #define TKCM_IMPORT __attribute__ ((visibility("default")))
    #else
        #define TKCM_EXPORT
        #define TKCM_IMPORT
    #endif

    #if defined(TKCM_BUILD_NODEDEF_DLL)
        #define TKCM_DECL TKCM_EXPORT
    #else
        #define TKCM_DECL TKCM_IMPORT
    #endif

    #define TKCM_DEFAULT_VALUE_DECL(NAME) \
        AMINO_DEFAULT_VAL_SHARED_DECLARE(TKCMInterface::NAME, AMINO_CONCAT(TKCMInterface__, NAME), TKCM_DECL)

#else
    #define TKCM_DECL
    #define TKCM_DEFAULT_VALUE_DECL( NAME )
#endif

#endif // TKCM_EXPORT_H