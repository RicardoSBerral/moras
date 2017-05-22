//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
#ifndef __LANGUAGE_H
#define __LANGUAGE_H
//---------------------------------------------------------------------------
#include <string>
#include <map>
//---------------------------------------------------------------------------

class StringRepository
{
   protected:
     std::map<const std::string, std::string> strings;
     static StringRepository *sr;

   protected:
     StringRepository() {}

   public:
     static StringRepository *GetStringRepository();

     static std::string GetString(std::string Key);

};
//---------------------------------------------------------------------------
class StringRepositorySpanish : public StringRepository
{
  public:
    StringRepositorySpanish();
};
//---------------------------------------------------------------------------
class StringRepositoryEnglish : public StringRepository
{
  public:
    StringRepositoryEnglish();
};
//---------------------------------------------------------------------------
#endif
