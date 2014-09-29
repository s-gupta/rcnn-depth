/*
 * Lambda abstractions.
 */
#ifndef LANG__LAMBDA_HH
#define LANG__LAMBDA_HH

#define LambdaXY(X,Y) X##Y
#define LambdaMakeNameXY(FX,LINE) LambdaXY(FX,LINE)
#define LambdaMakeName(FX) LambdaMakeNameXY(FX,__LINE__)

/* #define _ , */

#define lambda(args, ret_type, body) \
class LambdaMakeName(__lambda__) {\
public: \
   ret_type operator() args { body; } \
}

#define lambda_named(name, args, ret_type, body) \
class LambdaMakeName(name) { \
public: \
   ret_type operator() args { body; } \
}

#define lambda_derived(parents, args, ret_type, body) \
class LambdaMakeName(__lambda_derived__) parents {\
public: \
   ret_type operator() args { body; } \
}

#define lambda_named_derived(name, parents, args, ret_type, body) \
class LambdaMakeName(name) parents { \
public: \
   ret_type operator() args { body; } \
}

#endif
