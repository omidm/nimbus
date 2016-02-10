//#####################################################################
// Copyright 2004-2011, Eran Guendelman, Geoffrey Irving, Andrew Selle, Rahul Sheth.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace PROCESS_UTILITIES
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>
#if defined(WIN32)
#include <windows.h>
#include <psapi.h>
#include <tchar.h>
#include <StackWalker.h>
#pragma comment(lib,"psapi")
#elif defined(__linux__) || defined(__APPLE__)
#include <csignal>
#ifndef USE_OPENGLES
#include <execinfo.h>
#endif
#include <fenv.h>
#include <sys/resource.h>
#include <sys/time.h>
#endif
#ifdef __APPLE__
#include <xmmintrin.h>
#endif

namespace PhysBAM{
namespace PROCESS_UTILITIES{
//###################################################################
// Win32 Specific Function Definitions
//###################################################################
#if defined(WIN32)

unsigned int Memory_Usage()
{
    HANDLE process=GetCurrentProcess();
    PROCESS_MEMORY_COUNTERS counters;
    GetProcessMemoryInfo(process,&counters,sizeof(counters));
    return (unsigned int)counters.WorkingSetSize;
}

void Set_Floating_Point_Exception_Handling(const bool enable,const bool division_by_zero,const bool invalid_operation,
                                           const bool overflow,const bool underflow,const bool inexact_result)
{
    static bool have_original_state=false;
    static unsigned int saved_state=0;
    if(!have_original_state){
        _controlfp_s(&saved_state,0,0);
        have_original_state=true;}
    if(enable){
        unsigned int disabledExceptions=0;
        if(!division_by_zero) disabledExceptions|=_EM_ZERODIVIDE;
        if(!invalid_operation) disabledExceptions|=_EM_INVALID;
        if(!overflow) disabledExceptions|=_EM_OVERFLOW;
        if(!underflow) disabledExceptions|=_EM_UNDERFLOW;
        if(!inexact_result) disabledExceptions|=_EM_INEXACT;
        _clearfp();
        if(_controlfp_s(NULL,disabledExceptions,_MCW_EM)) PHYSBAM_FATAL_ERROR("Could not disable selected FPE exceptions");}
    else{
        if(!have_original_state || _controlfp_s(NULL,saved_state,_MCW_EM)) PHYSBAM_FATAL_ERROR("Could not restore default FPE exceptions");
        _clearfp();}
}

// From: Unhandled exceptions in VC8 and above for x86 and x64
// http://blog.kalmbach-software.de/2008/04/02/unhandled-exceptions-in-vc8-and-above-for-x86-and-x64/
#if defined _M_X64 || defined _M_IX86 //No support for IA64
LPTOP_LEVEL_EXCEPTION_FILTER WINAPI DummySetUnhandledExceptionFilter(LPTOP_LEVEL_EXCEPTION_FILTER lpTopLevelExceptionFilter) {return NULL;}

//Try to get exception before Windows Error Reporting (WER)
BOOL PreventSetUnhandledExceptionFilter(){
  HMODULE hKernel32 = LoadLibrary(_T("kernel32.dll"));
  if (hKernel32 == NULL) return FALSE;
  void *pOrgEntry = GetProcAddress(hKernel32, "SetUnhandledExceptionFilter");
  if(pOrgEntry == NULL) return FALSE;
  DWORD dwOldProtect = 0;
  SIZE_T jmpSize = 5;
#ifdef _M_X64
  jmpSize = 13;
#endif
  BOOL bProt = VirtualProtect(pOrgEntry, jmpSize, PAGE_EXECUTE_READWRITE, &dwOldProtect);
  BYTE newJump[20];
  void *pNewFunc = &DummySetUnhandledExceptionFilter;
#ifdef _M_IX86
  DWORD dwOrgEntryAddr = (DWORD) pOrgEntry;
  dwOrgEntryAddr += jmpSize; // add 5 for 5 op-codes for jmp rel32
  DWORD dwNewEntryAddr = (DWORD) pNewFunc;
  DWORD dwRelativeAddr = dwNewEntryAddr - dwOrgEntryAddr;
  // JMP rel32: Jump near, relative, displacement relative to next instruction.
  newJump[0] = 0xE9;  // JMP rel32
  memcpy(&newJump[1], &dwRelativeAddr, sizeof(pNewFunc));
#elif _M_X64
  newJump[0] = 0x49;  // MOV R15, ...
  newJump[1] = 0xBF;  // ...
  memcpy(&newJump[2], &pNewFunc, sizeof (pNewFunc));
  //pCur += sizeof (ULONG_PTR);
  newJump[10] = 0x41;  // JMP R15, ...
  newJump[11] = 0xFF;  // ...
  newJump[12] = 0xE7;  // ...
#endif
  SIZE_T bytesWritten;
  BOOL bRet = WriteProcessMemory(GetCurrentProcess(), pOrgEntry, newJump, jmpSize, &bytesWritten);
  if (bProt != FALSE){
    DWORD dwBuf;
    VirtualProtect(pOrgEntry, jmpSize, dwOldProtect, &dwBuf);}
  return bRet;
}

//override the built in output function to print to console
class StackWalkerToConsole : public StackWalker
{
public:
    StackWalkerToConsole() : StackWalker(0xF,
        NULL,GetCurrentProcessId(),GetCurrentProcess()){}
protected:
  virtual void OnOutput(LPCSTR szText){
    printf("%s", szText);}
};

static LONG __stdcall CrashHandlerExceptionFilter(EXCEPTION_POINTERS* pExPtrs)
{
#ifdef _M_IX86
  if (pExPtrs->ExceptionRecord->ExceptionCode == EXCEPTION_STACK_OVERFLOW){
    static char MyStack[1024*128];
    __asm mov eax,offset MyStack[1024*128];
    __asm mov esp,eax;
  }
#endif
  StackWalkerToConsole sw;
  sw.ShowCallstack(GetCurrentThread(), pExPtrs->ContextRecord);
  TCHAR lString[500];
  _stprintf_s(lString,
     _T("*** Unhandled Exception!\n")
     _T("   ExpCode: 0x%8.8X\n")
     _T("   ExpFlags: %d\n")
     _T("   ExpAddress: 0x%8.8X"),
     pExPtrs->ExceptionRecord->ExceptionCode,
     pExPtrs->ExceptionRecord->ExceptionFlags,
     pExPtrs->ExceptionRecord->ExceptionAddress);
  FatalAppExit(-1, lString);
  return EXCEPTION_CONTINUE_SEARCH;
}

//you can only call prevent unhandled once
static BOOL s_bUnhandledExceptionFilterSet=FALSE;

void Backtrace()
{
    LOG::cerr<<"=================== BEGIN STACK BACKTRACE ==================="<<std::endl;
    StackWalkerToConsole sw;
    sw.ShowCallstack();
    LOG::cerr<<"==================== END STACK BACKTRACE ===================="<<std::endl;
}

void Block_Interrupt_Signal(const bool block)
{
    PHYSBAM_WARNING("PROCESS_UTILITIES::Block_Interrupt_Signal undefined for windows");
}

void Set_Backtrace(const bool enable)
{
    if(enable){
        // set global exception handler (for handling all unhandled exceptions)
        SetUnhandledExceptionFilter(CrashHandlerExceptionFilter);
        if(s_bUnhandledExceptionFilterSet==FALSE){
    #if defined _M_X64 || defined _M_IX86
         //PreventSetUnhandledExceptionFilter(); //this is for disabling all default exception handlers
    #endif
         s_bUnhandledExceptionFilterSet=TRUE;}}
    else{
        SetUnhandledExceptionFilter(0);
    }
}

#else //IA64 not supported
void Backtrace()
{
    PHYSBAM_WARNING("PROCESS_UTILITIES::Backtrace undefined for windows");
}

void Block_Interrupt_Signal(const bool block)
{
    PHYSBAM_WARNING("PROCESS_UTILITIES::Block_Interrupt_Signal undefined for windows");
}

void Set_Backtrace(const bool enable)
{
    PHYSBAM_WARNING("PROCESS_UTILITIES::Set_Debug_Backtrace undefined for windows");
}
#endif

void PB_Sleep(const unsigned int milliseconds)
{
    Sleep(milliseconds);
}

//###################################################################
// Linux Specific Function Definitions
//###################################################################
#elif defined(__linux__) || defined(__APPLE__)

unsigned int Memory_Usage()
{
    struct rusage usage;
    getrusage(RUSAGE_SELF,&usage);
    return usage.ru_idrss+usage.ru_isrss;
}

static void Floating_Point_Exception_Handler(int sig_number,siginfo_t* info,void *data)
{
    if(sig_number!=SIGFPE) PHYSBAM_FATAL_ERROR();
    LOG::cerr<<"** ERROR: SIGNAL "<<"SIGFPE ("<<sig_number<<") **"<<std::endl;
    LOG::cerr<<"Floating point exception: reason "<<info->si_code<<" = \""<<
        (info->si_code==FPE_INTDIV?"integer divide by zero":info->si_code==FPE_INTOVF?"integer overflow":
        info->si_code==FPE_FLTDIV?"FP divide by zero":info->si_code==FPE_FLTOVF?"FP overflow":
        info->si_code==FPE_FLTUND?"FP underflow":info->si_code==FPE_FLTRES?"FP inexact result":
        info->si_code==FPE_FLTINV?"FP invalid operation":info->si_code==FPE_FLTSUB?"subscript out of range":"unknown")
        << "\", from address 0x"<<std::hex<<(unsigned long)info->si_addr<<std::endl;
    Backtrace();
    LOG::Finish_Logging();
    exit(sig_number);
}

#ifdef __APPLE__

// feenableexcept and fedisableexcept are not defined, so define them

static int Flags_To_Mask(int flags)
{
    int result=0;
    if(flags&FE_INVALID) result|=_MM_MASK_INVALID;
    if(flags&FE_DIVBYZERO) result|=_MM_MASK_DIV_ZERO;
    if(flags&FE_OVERFLOW) result|=_MM_MASK_OVERFLOW;
    if(flags&FE_UNDERFLOW) result|=_MM_MASK_UNDERFLOW;
    if(flags&FE_INEXACT) result|=_MM_MASK_INEXACT;
    return result;
}

static void fedisableexcept(int flags)
{
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() | Flags_To_Mask(flags));
}

static void feenableexcept(int flags)
{
    _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~Flags_To_Mask(flags));
}

#endif

void Set_Floating_Point_Exception_Handling(const bool enable,const bool division_by_zero,const bool invalid_operation,const bool overflow,const bool underflow,const bool inexact_result)
{
    static bool have_original_action=false;
    static struct sigaction original_action;
    if(!have_original_action){ // initialize with original action
        sigaction(SIGFPE,0,&original_action);}
    if(enable){
        int exceptions=0;
        if(division_by_zero) exceptions|=FE_DIVBYZERO;
        if(invalid_operation) exceptions|=FE_INVALID;
        if(overflow) exceptions|=FE_OVERFLOW;
        if(underflow) exceptions|=FE_UNDERFLOW;
        if(inexact_result) exceptions|=FE_INEXACT;
        // avoid catching delayed exceptions caused by external code
        fedisableexcept(FE_ALL_EXCEPT);
        feclearexcept(exceptions);
        // install new handler
        struct sigaction action;
        action.sa_flags=SA_SIGINFO;
        action.sa_sigaction=Floating_Point_Exception_Handler;
        sigemptyset(&action.sa_mask);
        if(sigaction(SIGFPE,&action,0)) PHYSBAM_FATAL_ERROR("Could not register FPE signal handler");
        feenableexcept(exceptions);}
    else{
        if(sigaction(SIGFPE,&original_action,0)) PHYSBAM_FATAL_ERROR("Could not restore FPE signal handler");
        fedisableexcept(FE_ALL_EXCEPT);}
}

void Backtrace()
{
    const int stack_entries=50;
    void *stack_array[stack_entries];
#ifdef USE_OPENGLES  //TODO: Add Statck trace ability for andoid
    size_t size=0;
    char **strings=0;
#else
    size_t size=backtrace(stack_array,stack_entries);
    char **strings=backtrace_symbols(stack_array,size);
#endif
    LOG::cerr<<"=================== BEGIN STACK BACKTRACE ==================="<<std::endl;
    int pid=getpid();
    static const unsigned int buf_size=2048;
    static char buffer[buf_size];
    sprintf(buffer,"/tmp/physbam-%d.stack",pid);
    {FILE* fp=fopen(buffer,"w");
    if(fp){
        for(size_t i=0;i<size;i++) fprintf(fp,"%s\n",strings[i]);
        fclose(fp);}}
    {sprintf(buffer,"c++filt < /tmp/physbam-%d.stack",pid);
    FILE* fp=popen(buffer,"r");
    if(fp){
        int depth=0;
        while(!feof(fp) && fgets(buffer,buf_size-1,fp)){
#ifdef __linux__
            // The format should be "library.so(function) ...", so parse up to '(' to skip the library
            char* paren=index(buffer,'(');
            if(paren) LOG::cerr<<"#"<<std::setiosflags(std::ios::left)<<std::setw(5)<<depth<<" "<<paren;
#else
            // We're probably on a Mac, and the format is already reasonable.  Just print the whole line.
            LOG::cerr<<buffer;
#endif
            depth+=1;
        }
        fclose(fp);}}
    LOG::cerr<<"==================== END STACK BACKTRACE ===================="<<std::endl;
    free(strings);
}

static bool caught_interrupt_signal=false;

static void Interrupt_Signal_Handler(int signal_id)
{
    caught_interrupt_signal=true;
}

void Block_Interrupt_Signal(const bool block)
{
    static bool have_original_action=false;
    static struct sigaction original_action;
    if(block){
        if(have_original_action) PHYSBAM_FATAL_ERROR("Nested call to Block_Interrupt_Signal(true).");
        struct sigaction action;
        action.sa_flags=0;
        action.sa_handler=Interrupt_Signal_Handler;
        sigemptyset(&action.sa_mask);
        if(sigaction(SIGINT,&action,&original_action)) PHYSBAM_FATAL_ERROR("Could not block interrupt signal.");
        have_original_action=true;}
    else{
        if(!have_original_action) PHYSBAM_FATAL_ERROR("Call to Block_Interrupt_Signal(false) before Block_Interrupt_Signal(true).");
        if(sigaction(SIGINT,&original_action,0)) PHYSBAM_FATAL_ERROR("Could not unblock interrupt signal.");
        if(caught_interrupt_signal){
            LOG::cerr<<"Caught delayed interrupt signal."<<std::endl;
            raise(SIGINT);}
        have_original_action=false;}
}

static int physbam_catch_signals[]={SIGINT,SIGABRT,SIGSEGV,SIGBUS,SIGTERM,SIGHUP,SIGUSR2,0};
static const char* physbam_catch_signal_names[]={"SIGINT","SIGABRT","SIGSEGV","SIGBUS","SIGTERM","SIGHUP","SIGUSR2",0};

void Backtrace_And_Abort(int signal_number,siginfo_t* info,void *data)
{
    std::stringstream ss;ss<<std::flush;LOG::filecout(ss.str());
    LOG::cerr<<"\n";
    Backtrace();
    const char** names=physbam_catch_signal_names;const char *signal_name=0;
    for(int *i=physbam_catch_signals;*i!=0;i++,names++) if(signal_number==*i) signal_name=*names;
    LOG::cerr<<"\n*** ERROR: SIGNAL "<<(signal_name?signal_name:"UNKNOWN")<<" ("<<signal_number<<")\n"<<std::endl;
    if(signal_number!=SIGUSR2){
        LOG::Finish_Logging();
        exit(signal_number);}
}

void Set_Backtrace(const bool enable)
{
    if(enable){
        struct sigaction action;
        action.sa_flags=SA_SIGINFO;
        action.sa_sigaction=Backtrace_And_Abort;
        sigemptyset(&action.sa_mask);
        for(int *i=physbam_catch_signals;*i!=0;i++) sigaddset(&action.sa_mask,*i);
        for(int *i=physbam_catch_signals;*i!=0;i++) if(sigaction(*i,&action,0)) PHYSBAM_FATAL_ERROR("Failed to install backtrace handler.");}
    else for(int *i=physbam_catch_signals;*i!=0;i++) signal(*i,SIG_DFL);
}

void PB_Sleep(const unsigned int milliseconds)
{
    sleep(milliseconds);
}

//###################################################################
// Default (Unimplemented) Function Definitions
//###################################################################
#else

unsigned int Memory_Usage()
{
    LOG::cerr<<"PROCESS_UTILITIES::Memory_Usage undefined"<<std::endl;
    return 0;
}

void Set_Floating_Point_Exception_Handling(const bool enable,const bool division_by_zero,const bool invalid_operation,
                                           const bool overflow,const bool underflow,const bool inexact_result)
{
    LOG::cerr<<"PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling undefined"<<std::endl;
}

void Backtrace()
{
    LOG::cerr<<"PROCESS_UTILITIES::Backtrace undefined"<<std::endl;
}

void Block_Interrupt_Signal(const bool block)
{
    LOG::cerr<<"PROCESS_UTILITIES::Block_Interrupt_Signal undefined"<<std::endl;
}

void Set_Backtrace(const bool enable)
{
    LOG::cerr<<"PROCESS_UTILITIES::Set_Debug_Backtrace undefined"<<std::endl;
}

void PB_Sleep(const unsigned int milliseconds)
{
    PHYSBAM_WARNING("PROCESS_UTILITIES::Sleep undefined");
}

#endif

//###################################################################
}
}
