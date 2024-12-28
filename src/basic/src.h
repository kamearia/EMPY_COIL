#pragma once
/**/
#if defined (SSIL_LINUX) || defined (FOCUS_LINUX)
#define CPP_SRC
#endif

#ifdef CPP_SRC
#define C_HEADER "C"
#define C_HEADER_IN_CPP
#define C_START "C" {
#define C_END }
#else
#define C_HEADER
#define C_HEADER_IN_CPP "C"
#define C_START  
#define C_END 
#endif

