// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME nanoDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "JToyMCTrack.h"
#include "JHistos.h"
#include "JInputs.h"
#include "JEventLists.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_JToyMCTrack(void *p = 0);
   static void *newArray_JToyMCTrack(Long_t size, void *p);
   static void delete_JToyMCTrack(void *p);
   static void deleteArray_JToyMCTrack(void *p);
   static void destruct_JToyMCTrack(void *p);
   static void streamer_JToyMCTrack(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::JToyMCTrack*)
   {
      ::JToyMCTrack *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::JToyMCTrack >(0);
      static ::ROOT::TGenericClassInfo 
         instance("JToyMCTrack", ::JToyMCTrack::Class_Version(), "JToyMCTrack.h", 7,
                  typeid(::JToyMCTrack), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::JToyMCTrack::Dictionary, isa_proxy, 16,
                  sizeof(::JToyMCTrack) );
      instance.SetNew(&new_JToyMCTrack);
      instance.SetNewArray(&newArray_JToyMCTrack);
      instance.SetDelete(&delete_JToyMCTrack);
      instance.SetDeleteArray(&deleteArray_JToyMCTrack);
      instance.SetDestructor(&destruct_JToyMCTrack);
      instance.SetStreamerFunc(&streamer_JToyMCTrack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::JToyMCTrack*)
   {
      return GenerateInitInstanceLocal((::JToyMCTrack*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::JToyMCTrack*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static TClass *JHistos_Dictionary();
   static void JHistos_TClassManip(TClass*);
   static void *new_JHistos(void *p = 0);
   static void *newArray_JHistos(Long_t size, void *p);
   static void delete_JHistos(void *p);
   static void deleteArray_JHistos(void *p);
   static void destruct_JHistos(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::JHistos*)
   {
      ::JHistos *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::JHistos));
      static ::ROOT::TGenericClassInfo 
         instance("JHistos", "JHistos.h", 11,
                  typeid(::JHistos), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &JHistos_Dictionary, isa_proxy, 0,
                  sizeof(::JHistos) );
      instance.SetNew(&new_JHistos);
      instance.SetNewArray(&newArray_JHistos);
      instance.SetDelete(&delete_JHistos);
      instance.SetDeleteArray(&deleteArray_JHistos);
      instance.SetDestructor(&destruct_JHistos);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::JHistos*)
   {
      return GenerateInitInstanceLocal((::JHistos*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::JHistos*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *JHistos_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::JHistos*)0x0)->GetClass();
      JHistos_TClassManip(theClass);
   return theClass;
   }

   static void JHistos_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *JInputs_Dictionary();
   static void JInputs_TClassManip(TClass*);
   static void *new_JInputs(void *p = 0);
   static void *newArray_JInputs(Long_t size, void *p);
   static void delete_JInputs(void *p);
   static void deleteArray_JInputs(void *p);
   static void destruct_JInputs(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::JInputs*)
   {
      ::JInputs *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::JInputs));
      static ::ROOT::TGenericClassInfo 
         instance("JInputs", "JInputs.h", 37,
                  typeid(::JInputs), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &JInputs_Dictionary, isa_proxy, 0,
                  sizeof(::JInputs) );
      instance.SetNew(&new_JInputs);
      instance.SetNewArray(&newArray_JInputs);
      instance.SetDelete(&delete_JInputs);
      instance.SetDeleteArray(&deleteArray_JInputs);
      instance.SetDestructor(&destruct_JInputs);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::JInputs*)
   {
      return GenerateInitInstanceLocal((::JInputs*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::JInputs*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *JInputs_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::JInputs*)0x0)->GetClass();
      JInputs_TClassManip(theClass);
   return theClass;
   }

   static void JInputs_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *JEventLists_Dictionary();
   static void JEventLists_TClassManip(TClass*);
   static void *new_JEventLists(void *p = 0);
   static void *newArray_JEventLists(Long_t size, void *p);
   static void delete_JEventLists(void *p);
   static void deleteArray_JEventLists(void *p);
   static void destruct_JEventLists(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::JEventLists*)
   {
      ::JEventLists *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::JEventLists));
      static ::ROOT::TGenericClassInfo 
         instance("JEventLists", "JEventLists.h", 6,
                  typeid(::JEventLists), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &JEventLists_Dictionary, isa_proxy, 0,
                  sizeof(::JEventLists) );
      instance.SetNew(&new_JEventLists);
      instance.SetNewArray(&newArray_JEventLists);
      instance.SetDelete(&delete_JEventLists);
      instance.SetDeleteArray(&deleteArray_JEventLists);
      instance.SetDestructor(&destruct_JEventLists);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::JEventLists*)
   {
      return GenerateInitInstanceLocal((::JEventLists*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::JEventLists*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *JEventLists_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::JEventLists*)0x0)->GetClass();
      JEventLists_TClassManip(theClass);
   return theClass;
   }

   static void JEventLists_TClassManip(TClass* ){
   }

} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr JToyMCTrack::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *JToyMCTrack::Class_Name()
{
   return "JToyMCTrack";
}

//______________________________________________________________________________
const char *JToyMCTrack::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::JToyMCTrack*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int JToyMCTrack::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::JToyMCTrack*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *JToyMCTrack::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::JToyMCTrack*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *JToyMCTrack::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::JToyMCTrack*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void JToyMCTrack::Streamer(TBuffer &R__b)
{
   // Stream an object of class JToyMCTrack.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      lVec.Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, JToyMCTrack::IsA());
   } else {
      R__c = R__b.WriteVersion(JToyMCTrack::IsA(), kTRUE);
      TObject::Streamer(R__b);
      lVec.Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_JToyMCTrack(void *p) {
      return  p ? new(p) ::JToyMCTrack : new ::JToyMCTrack;
   }
   static void *newArray_JToyMCTrack(Long_t nElements, void *p) {
      return p ? new(p) ::JToyMCTrack[nElements] : new ::JToyMCTrack[nElements];
   }
   // Wrapper around operator delete
   static void delete_JToyMCTrack(void *p) {
      delete ((::JToyMCTrack*)p);
   }
   static void deleteArray_JToyMCTrack(void *p) {
      delete [] ((::JToyMCTrack*)p);
   }
   static void destruct_JToyMCTrack(void *p) {
      typedef ::JToyMCTrack current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_JToyMCTrack(TBuffer &buf, void *obj) {
      ((::JToyMCTrack*)obj)->::JToyMCTrack::Streamer(buf);
   }
} // end of namespace ROOT for class ::JToyMCTrack

namespace ROOT {
   // Wrappers around operator new
   static void *new_JHistos(void *p) {
      return  p ? new(p) ::JHistos : new ::JHistos;
   }
   static void *newArray_JHistos(Long_t nElements, void *p) {
      return p ? new(p) ::JHistos[nElements] : new ::JHistos[nElements];
   }
   // Wrapper around operator delete
   static void delete_JHistos(void *p) {
      delete ((::JHistos*)p);
   }
   static void deleteArray_JHistos(void *p) {
      delete [] ((::JHistos*)p);
   }
   static void destruct_JHistos(void *p) {
      typedef ::JHistos current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::JHistos

namespace ROOT {
   // Wrappers around operator new
   static void *new_JInputs(void *p) {
      return  p ? new(p) ::JInputs : new ::JInputs;
   }
   static void *newArray_JInputs(Long_t nElements, void *p) {
      return p ? new(p) ::JInputs[nElements] : new ::JInputs[nElements];
   }
   // Wrapper around operator delete
   static void delete_JInputs(void *p) {
      delete ((::JInputs*)p);
   }
   static void deleteArray_JInputs(void *p) {
      delete [] ((::JInputs*)p);
   }
   static void destruct_JInputs(void *p) {
      typedef ::JInputs current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::JInputs

namespace ROOT {
   // Wrappers around operator new
   static void *new_JEventLists(void *p) {
      return  p ? new(p) ::JEventLists : new ::JEventLists;
   }
   static void *newArray_JEventLists(Long_t nElements, void *p) {
      return p ? new(p) ::JEventLists[nElements] : new ::JEventLists[nElements];
   }
   // Wrapper around operator delete
   static void delete_JEventLists(void *p) {
      delete ((::JEventLists*)p);
   }
   static void deleteArray_JEventLists(void *p) {
      delete [] ((::JEventLists*)p);
   }
   static void destruct_JEventLists(void *p) {
      typedef ::JEventLists current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::JEventLists

namespace {
  void TriggerDictionaryInitialization_nanoDict_Impl() {
    static const char* headers[] = {
"JToyMCTrack.h",
"JHistos.h",
"JInputs.h",
"JEventLists.h",
0
    };
    static const char* includePaths[] = {
"/home/heimarry/alice/sw/slc7_x86-64/ROOT/v6-10-08-1/include",
"/home/heimarry/Documents/simulaatioita/ToyFlow/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "nanoDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$JToyMCTrack.h")))  JToyMCTrack;
class __attribute__((annotate("$clingAutoload$JHistos.h")))  JHistos;
class __attribute__((annotate("$clingAutoload$JInputs.h")))  JInputs;
class __attribute__((annotate("$clingAutoload$JEventLists.h")))  JEventLists;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "nanoDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif
#ifndef toka
  #define toka 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "JToyMCTrack.h"
#include "JHistos.h"
#include "JInputs.h"
#include "JEventLists.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"JEventLists", payloadCode, "@",
"JHistos", payloadCode, "@",
"JInputs", payloadCode, "@",
"JToyMCTrack", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("nanoDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_nanoDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_nanoDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_nanoDict() {
  TriggerDictionaryInitialization_nanoDict_Impl();
}
