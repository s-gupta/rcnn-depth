// Copyright (C) 2002 David R. Martin <dmartin@eecs.berkeley.edu>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of the
// License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
// 02111-1307, USA, or see http://www.gnu.org/copyleft/gpl.html.
    
#ifndef REGEX_HH
#define REGEX_HH
    
// A C++ wrapper for the GNU regex routines.
    
#include <regex.h>

namespace Util
{
  class Regex  
  {

    public:

      static const int cflagsDefault = REG_EXTENDED | REG_NEWLINE;
      
          // Construct an uncompiled RE.
          Regex();
      
          // Construct a compiled RE.
          Regex(const char *re, int cflags = cflagsDefault);
      
          // Destructor.
          ~Regex();
      
          // Compile or recompile an RE.
      void compile(const char *re, int cflags = cflagsDefault);
      
          // Substitutions.  Modifies argument.  Returns modified argument.
          // Caller must free text.  Argument must be free()able.
          // The RE must be compiled.  
      char *subst(char *&text, const char *rep, bool global = true);
      
          // Returns matched().
          // The RE must be compiled.  
          bool match(const char *text, int eflags = 0);
      
          // Return the number of sub-matches.  
          // There is always >= 1 sub-match.  
          // The RE must be matched.
      unsigned count() const;
      
          // The array 'match' should be at least the length of the original
          // text to be safe.
          // The i'th sub-match is written into 'match'.
          // Returns match.
          // The RE must be matched.
      char *get(unsigned i, char *match) const;
      
          // The starting/ending index of the i'th match, and its length.
          // The RE must be matched.
      unsigned start(unsigned i) const;
      unsigned end(unsigned i) const;
      unsigned length(unsigned i) const;
      
          // Is this RE matched or compiled?
          bool matched()const {
          return _matched;
      }
      bool compiled() const {
          return _compiled;
      }
      
          // Convenience routines.
      
          // Compile and match.
       bool match(const char *text, const char *re, int eflags = 0) {
          compile(re);
          return match(text, eflags);
      }
      
          // One-time matching with no sub-matches.
      static bool matches(const char *text, const char *re, int eflags = 0) {
          Regex regex(re);
          return regex.match(text, eflags);
      }
      
          // One-time substitutions.
      static char *subst(char *&text, const char *re, const char *rep, bool global = true) 
      {
          Regex regex(re);
          return regex.subst(text, rep, global);
      } 
    
    protected:

      void _init();
      void _zero();
      void _delete();
      void _checkIndex(const char *origin, unsigned i) const;
      static const unsigned _nmatchDefault = 8;
      bool _compiled;
      bool _matched;
      unsigned _nmatch;         // Initialized by construction.
      regmatch_t *_pmatch;
      char *_regex;             // Initialized by compiling.
      regex_t _preg;
      char *_text;              // Initialized by matching.
      unsigned _textLen;
      unsigned _numMatches;
  };

} //namespace Util

#endif                          // __Regex_h__
