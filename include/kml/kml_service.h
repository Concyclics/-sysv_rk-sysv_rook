/*******************************************************************************
 * Copyright (c) Huawei Technologies Co., Ltd. 2020-2020. All rights reserved.
 * Description: Part of KML library
 * Author: KML
 * Create: 2020
 ******************************************************************************/
#ifndef LIBKSERVICE_KML_SERVICE_H
#define LIBKSERVICE_KML_SERVICE_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "kml_export.h"

#ifdef kml_complex_float

#elif defined(__cplusplus)

#include <complex>
#define kml_complex_float std::complex<float>

#elif __STDC_VERSION__ >= 199901L && !defined(__STDC_NO_COMPLEX__)

#define kml_complex_float float _Complex

#else

typedef struct kml_complex_float {
    float real;
    float imag;
} kml_complex_float;

#define kml_complex_float kml_complex_float

#endif

#ifdef kml_complex_double

#elif defined(__cplusplus)

#include <complex>
#define kml_complex_double std::complex<double>

#elif __STDC_VERSION__ >= 199901L && !defined(__STDC_NO_COMPLEX__)

#define kml_complex_double double _Complex

#else

typedef struct kml_complex_double {
    double real;
    double imag;
} kml_complex_double;

#define kml_complex_double kml_complex_double

#endif

#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus

template<class T, class U>
constexpr T KmlMin(T a, U b)
{
    static_assert(std::is_same<T, U>::value, "KmlMin doesn't allow implicit casts on arguments.");
    if constexpr (std::is_floating_point<T>::value) {
        if (std::isnan(a)) {
            return a;
        }
    }
    return (a < b) ? a : b;
}

template<class T, class U>
constexpr T KmlMax(T a, U b)
{
    static_assert(std::is_same<T, U>::value, "KmlMax doesn't allow implicit casts on arguments.");
    if constexpr (std::is_floating_point<T>::value) {
        if (std::isnan(a)) {
            return a;
        }
    }
    return (a > b) ? a : b;
}

extern "C" {
#endif

#define KML_DEFAULT_ALLOCATE_ALIGN 128

#ifndef KML_ATTR_MALLOC
#define KML_ATTR_MALLOC __attribute__((malloc))
#endif

#ifndef KML_ATTR_WEAK
#define KML_ATTR_WEAK __attribute__((weak))
#endif

#ifndef KML_ATTR_PUBLIC
#define KML_ATTR_PUBLIC __attribute__((visibility("default")))
#endif

#ifndef KML_ATTR_ALLOC_SIZE
#define KML_ATTR_ALLOC_SIZE(...) __attribute__((alloc_size(__VA_ARGS__)))
#endif

#ifndef KML_ATTR_UNUSED
#define KML_ATTR_UNUSED __attribute__((unused))
#endif

#ifndef KML_ATTR_NORETURN
#define KML_ATTR_NORETURN __attribute__((noreturn))
#endif

#ifndef KML_ATTR_ALWAYS_INLINE
#define KML_ATTR_ALWAYS_INLINE __attribute__((always_inline))
#endif

#ifndef KML_THREAD_LOCAL
#define KML_THREAD_LOCAL __thread
#endif

#ifndef KML_ATTR_CONSTRUCTOR
#define KML_ATTR_CONSTRUCTOR(x) __attribute__((constructor((x))))
#endif

#ifndef KML_FALLTHROUGH
#define KML_FALLTHROUGH __attribute__((fallthrough))
#endif

#ifndef L1_CACHELINE_SIZE
#define L1_CACHELINE_SIZE 64
#endif

#ifndef L2_CACHELINE_SIZE
#define L2_CACHELINE_SIZE 64
#endif

#ifndef L3_CACHELINE_SIZE
#define L3_CACHELINE_SIZE 128
#endif

#ifndef KML_WAIL
#define KML_WAIL(why, ...)                                                   \
    do {                                                                     \
        KmlCannotContinue("%s:%i: " why, __FILE__, __LINE__, ##__VA_ARGS__); \
    } while (0)
#endif

#ifndef KML_WAIL_NO_MEMORY
#define KML_WAIL_NO_MEMORY KML_WAIL("OUT OF MEMORY")
#endif

#ifndef LIKELY
#define LIKELY(x) __builtin_expect(!!(x), 1)
#endif

#ifndef UNLIKELY
#define UNLIKELY(x) __builtin_expect(!!(x), 0)
#endif

#ifndef SELDOM
#define SELDOM UNLIKELY
#endif

#ifndef OFTEN
#define OFTEN LIKELY
#endif

#ifndef KML_ATTR_ALIGNED
#define KML_ATTR_ALIGNED(x) __attribute__((aligned((x))))
#endif

#ifndef KML_ATTR_RESTRICT
#define KML_ATTR_RESTRICT __restrict__
#endif

#define KML_VERSION_PART_LEN 100
struct KMLVersion {
    char component[KML_VERSION_PART_LEN];
    char version[KML_VERSION_PART_LEN];
    char supportPackage[KML_VERSION_PART_LEN];
    char compileTime[KML_VERSION_PART_LEN];
};
typedef struct KMLVersion KSERVICEVersion;

extern KML_EXPORT int KSERVICEGetVersion(KSERVICEVersion *ver);

enum KmlSetting {
    KML_VERBOSE,
    KML_NUM_THREADS,
    KML_THREAD_USING_POLICY,
    KML_CFMA_ACCELERATION,
    KML_CFMA_DRIVER_LOADED,
    KML_BLAS_SKIP_ZEROS,
    KML_GEMM_TO_GEMM3M,
    KML_BITWISE_COMPATIBILITY,
    KML_FAST_IGNORE_DENORMS,
    KML_ALLOCATE_POLICY,
    KML_ALLOCATE_HUGEPAGE_SIZE,
    SETTINGS_FLAG_COUNT
};

enum SettingsStates { SETTINGS_ENTRY_NOT_INITED, SETTINGS_ENTRY_INITED };

union KmlSettingValue {
    size_t sizetValue;
    double doubleValue;
    char *stringValue;
    int intValue;
};

struct KmlSettingEntry {
    union KmlSettingValue value;
    enum SettingsStates entryState;
};

enum KmlThreadUsingPolicyOption {
    KML_USE_ALL_THREADS = 0,
    KML_THRESHOLD_ENABLED = 1,
    KML_CROSSPOINTS_ENABLED = 2,
};

enum KmlBlasSkipZerosOption {
    KML_BLAS_NO_ZERO_SKIP = 0,
    KML_BLAS_NO_INTRUSIVE_ZERO_SKIP = 1,
    KML_BLAS_ENFORCED_ZERO_SKIP = 2,
};

enum KmlGemmToGemM3MOption {
    KML_GEMM_TO_GEMM_3M_NEVER = 0,
    KML_GEMM_TO_GEMM_3M_OPTIMIZED = 1,
    KML_GEMM_TO_GEMM_3M_ALWAYS = 2,
};

enum KmlAllocatePolicy {
    KML_HUGEPAGES_NEVER = 0,
    KML_HUGEPAGES_ALWAYS = 1,
    KML_HUGEPAGES_AS_HEAP = 2,
};

struct KmlHugePagesInfo {
    size_t size;
    size_t mmaped;
    char *heaptop;
    char *heapbase;
};

// defined in memory.c
extern struct KmlHugePagesInfo g_kmlHugepages;

#ifndef TO_UPPER_CONSTANT
#define TO_UPPER_CONSTANT (~32)
#endif

KML_ATTR_ALWAYS_INLINE static inline char ToUpper(const char a)
{
    return (char)((unsigned char)a & (unsigned char)TO_UPPER_CONSTANT);
}

KML_ATTR_ALWAYS_INLINE static inline int KmlLsame(const char a, const char b)
{
    return ToUpper(a) == ToUpper(b);
}

#ifndef __cplusplus

#define AUTOTYPE_MIN(a, b)               \
    ({                                   \
        __auto_type const _a = (a);      \
        __auto_type const _b = (b);      \
        (_a != _a || _a < _b) ? _a : _b; \
    })
#define KmlMin(a, b)                                                                                                   \
    ({                                                                                                                 \
        _Static_assert(__builtin_types_compatible_p(typeof((a)), typeof((b))), "KmlMin doesn't allow implicit casts"); \
        AUTOTYPE_MIN((a), (b));                                                                                        \
    })

#define AUTOTYPE_MAX(a, b)               \
    ({                                   \
        __auto_type const _a = (a);      \
        __auto_type const _b = (b);      \
        (_a != _a || _a > _b) ? _a : _b; \
    })
#define KmlMax(a, b)                                                                                                   \
    ({                                                                                                                 \
        _Static_assert(__builtin_types_compatible_p(typeof((a)), typeof((b))), "KmlMax doesn't allow implicit casts"); \
        AUTOTYPE_MAX((a), (b));                                                                                        \
    })

#endif

KML_ATTR_ALWAYS_INLINE static inline float Int2FloatCeil(const int i)
{
    float rv = i;
    if (i > (int)rv) {
        int *irv = (int *)&rv;
        ++*(irv);
    }
    return rv;
}

KML_ATTR_ALWAYS_INLINE static inline int NonzeroDenom(int x)
{
    return x ? x : 1;
}

KML_ATTR_ALWAYS_INLINE static inline long long int KmlAlign(long long int val, long long int align)
{
    return (val + align - 1) & (~(align - 1));
}

// atomic.c
// Do not be confused by "atomic" types. Use atomic operations.
// For non-atomic operations, like initialization, use "atoms"
// of atomic variables: atomicValue.a = 0;

typedef struct {
    int a;
} atomic_int;
typedef struct {
    int64_t a;
} atomic_int64;
typedef struct {
    uint64_t a;
} atomic_uint64;

extern KML_EXPORT int AtomicLoad(atomic_int *p);
extern KML_EXPORT int AtomicLoadRelaxed(atomic_int *p);
extern KML_EXPORT int64_t AtomicLoadInt64(atomic_int64 *p);
extern KML_EXPORT void AtomicStore(atomic_int *p, int value);
extern KML_EXPORT void AtomicStoreInt64(atomic_int64 *p, int64_t value);
extern KML_EXPORT int AtomicInc(atomic_int *p);
extern KML_EXPORT int64_t AtomicIncInt64(atomic_int64 *p);
extern KML_EXPORT int AtomicDec(atomic_int *p);
extern KML_EXPORT int64_t AtomicDecInt64(atomic_int64 *p);
extern KML_EXPORT int AtomicCmpxchg(atomic_int *p, int *expected, int desired);
extern KML_EXPORT int AtomicCmpxchgInt64(atomic_int64 *p, int64_t *expected, int64_t desired);

typedef struct barrier {
    KML_ATTR_ALIGNED(L1_CACHELINE_SIZE) atomic_int64 counter;
    KML_ATTR_ALIGNED(L1_CACHELINE_SIZE) atomic_int64 fence;
    int nthr, nb, nInit;
} barrier_t;

extern KML_EXPORT void BarrierInit(barrier_t *b, int nthrTeam);
extern KML_EXPORT void BarrierInit2(barrier_t *b, int nthrTeam, int nb, int nInit);
extern KML_EXPORT void BarrierWait(barrier_t *b, int ithrTeam);

// Array Division
typedef struct KmlArrayDivisionInfo {
    int ithr;
    int nthr;
    int size;
    int align;
} KmlArrayDivisionInfo;

extern KML_EXPORT void KmlArrayDivisionInfoInit(KmlArrayDivisionInfo *info, int ithr, int nthr, int size, int align);
extern KML_EXPORT void KmlArrayDivision(const KmlArrayDivisionInfo *arraysInfo, int *chunk, int *ind);

// settings.c
extern const union KmlSettingValue *KmlGetFlag(enum KmlSetting flag);
extern union KmlSettingValue KmlGetFlagOrDefault(enum KmlSetting flag);
extern KML_EXPORT enum KmlThreadUsingPolicyOption KmlGetThreadUsingPolicy();

// memory.c
extern KML_EXPORT void *KML_ATTR_MALLOC KML_ATTR_ALLOC_SIZE(1) KmlAllocate(size_t sz, size_t align);
extern KML_EXPORT void *KML_ATTR_MALLOC KML_ATTR_ALLOC_SIZE(1, 2) KmlCalloc(size_t num, size_t size, size_t align);
extern KML_EXPORT void KmlTypedFree(char *p);
extern KML_EXPORT void KmlFree(void *p);

// timer.c
extern KML_EXPORT double dsecnd(void);
extern KML_EXPORT float second(void);
extern KML_EXPORT float second_(void);
extern KML_EXPORT double dsecnd_(void);

// parallel.c
extern KML_EXPORT int KmlGetNumProcs(void);
extern KML_EXPORT int KmlGetMaxThreads(void);
extern KML_EXPORT void KmlSetNumThreads(int numThreads);

// utils.c
extern KML_EXPORT void KmlCannotContinue(const char *msg, ...);
extern KML_EXPORT void KmlInform(const char *msg, ...);
extern KML_EXPORT int Lsame(const char *a, const char *b);
extern KML_EXPORT int lsame_(const char *a, const char *b);
extern KML_EXPORT int Xerbla(const char *func, const int *param, int funcLength);
extern KML_EXPORT int xerbla_(const char *func, const int *param, int funcLength);
extern KML_EXPORT double Dlamch(const char *c);
extern KML_EXPORT double dlamch_(const char *c);
extern KML_EXPORT float Slamch(const char *c);
extern KML_EXPORT float slamch_(const char *c);

// version.c
extern KML_EXPORT const char *KmlGetBuildConfig(void);

#define KML_WORK_MEM_ERROR (-1313)
#define KML_NULL_PTR_ERROR (100000)

#define KML_CHECK_PARAM_FOR_NULL_PTR(param, paramNum, funcName, info, ...) \
    do {                                                                   \
        if (UNLIKELY((param) == NULL)) {                                   \
            if ((info) != NULL) {                                          \
                *(info) = -(paramNum);                                     \
            }                                                              \
            int funcLen = (sizeof(funcName) / sizeof(char)) - 1;           \
            int xerblaInfo = KML_NULL_PTR_ERROR + (paramNum);              \
            Xerbla((funcName), &xerblaInfo, funcLen);                      \
            return __VA_ARGS__;                                            \
        }                                                                  \
    } while (0)

#ifdef __cplusplus
}
#endif

#endif // #ifndef LIBKSERVICE_KMLSERVICE_H
