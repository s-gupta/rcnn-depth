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
    
#include <stdlib.h>
#include <regex.h>
#include <string.h>
#include <assert.h>
#include "regex.hh"
#include "exception.hh"
#include "string.hh"

namespace Util
{
  void  Regex::_init() 
  {
      _nmatch = _nmatchDefault;
      _pmatch = new regmatch_t[_nmatch];
  } void  Regex::_zero() 
  {
      _compiled = false;
      _matched = false;
      _regex = NULL;
      _text = NULL;
      _pmatch = NULL;
  } void  Regex::_delete() 
  {
      if (_compiled) {
          regfree(&_preg);
      }
      if (_pmatch != NULL) {
          delete[]_pmatch;
      }
      if (_text != NULL) {
          free(_text);
      }
      if (_regex != NULL) {
          free(_regex);
      }
  }

  Regex::Regex() 
  {
      _zero();
      _init();
  }

  Regex::~Regex() 
  {
      _delete();
  }
  Regex::Regex(const char *re, int cflags) 
  {
      _zero();
      try {
          _init();
          compile(re, cflags);
      } catch(Exception & e) {
          _delete();
          throw;
  } } void  Regex::compile(const char *re, int cflags) 
  {
      _matched = false;
      if (_compiled) {
          regfree(&_preg);
          _compiled = false;
          if (_text != NULL) {
              free(_text);
              _text = NULL;
          }
          if (_regex != NULL) {
              free(_regex);
              _regex = NULL;
          }
      }
      int rval = regcomp(&_preg, re, cflags);
      _compiled = true;
      if (rval != 0) {
          static const int errbufsize = 1024;
          static char errbuf[errbufsize];
          regerror(rval, &_preg, errbuf, errbufsize);
          throw Exception(String("Regex::Regex: %s", errbuf));
      }
      _regex = strdup(re);
      _text = NULL;
  } bool  Regex::match(const char *text, int eflags) 
  {
      if (!_compiled) {
          throw Exception("Regex::match: Regex not compiled!");
      }
      
          // Try matching the RE to text.
          _pmatch[_nmatch - 1].rm_so = -1;
      int rval = regexec(&_preg, text, _nmatch, _pmatch, eflags);
      _matched = (rval == 0);
      
          // Return right away if there was no match.
          if (!_matched) {
          return false;
      }
      
          // Make sure we have enough space for all the matches.
          if (_pmatch[_nmatch - 1].rm_so != -1) {
          fprintf(stderr, "resize: %d-->%d\n", _nmatch, _nmatch * 2);
          delete[]_pmatch;
          _nmatch *= 2;
          _pmatch = new regmatch_t[_nmatch];
          return match(text);
      }
      
          // Save a copy of the text.
          if (_text != NULL) {
          free(_text);
      }
      _text = strdup(text);
      _textLen = strlen(_text);
      
          // Count the number of matches.
          _numMatches = 0;
      for (unsigned i = 0; i < _nmatch; i++) {
          if (_pmatch[i].rm_so == -1) {
              break;
          }
          _numMatches++;
      }
      assert(_numMatches > 0);
      return true;
  }
  void  Regex::_checkIndex(const char *origin, unsigned i) const 
  {
      if (!_matched) {
          throw
              Exception(String
                        ("Regex::%s: /%s/ has not been matched!", origin,
                         _regex));
      }
      if (i >= _nmatch || _pmatch[i].rm_so == -1) {
          throw
              Exception(String
                        ("Regex::%s: Request for sub-match %u of /%s/ exceeds number "
                          "of sub-matches (%u).", origin, i, _regex,
                         _numMatches));
      }
      assert(_text != NULL);
      assert(_pmatch[i].rm_so >= 0);
      assert(_pmatch[i].rm_eo >= _pmatch[i].rm_so);
      assert((unsigned) _pmatch[i].rm_so <= _textLen);
      assert((unsigned) _pmatch[i].rm_eo <= _textLen);
  } char * Regex::get(unsigned i, char *match) const 
  {
      _checkIndex("get", i);
      unsigned matchLen = _pmatch[i].rm_eo - _pmatch[i].rm_so;
      strncpy(match, _text + _pmatch[i].rm_so, matchLen);
      match[matchLen] = '\0';
      return match;
  }
  unsigned  Regex::count() const 
  {
      _checkIndex("count", 0);
      return _numMatches;
  }
  unsigned  Regex::start(unsigned i) const 
  {
      _checkIndex("start", i);
      return _pmatch[i].rm_so;
  }
  unsigned  Regex::end(unsigned i) const 
  {
      _checkIndex("end", i);
      return _pmatch[i].rm_eo;
  }
  unsigned  Regex::length(unsigned i) const 
  {
      _checkIndex("length", i);
      return _pmatch[i].rm_eo - _pmatch[i].rm_so;
  }
  char * Regex::subst(char *&s, const char *rep, bool global) 
  {
      unsigned repLen = strlen(rep);
      unsigned srcLen = strlen(s);
      unsigned last = 0;         // Loop invariant: We have processed s[0..last).
      String dest;
      for (unsigned i = 0;; i++) {
          if (i > 0 && !global) {
              break;
          }
          bool matched = match(s + last);
          if (!matched) {
              break;
          }
          unsigned a = start(0);
          unsigned b = end(0);
          
              // Copy src string up to the matched text.
              dest.append(a, s + last);
          
              // Skip the matched text.
              last += b;
          
              // Copy the replacement text instead of the matched text.
              dest.append(repLen, rep);
      } 
          // Copy the remainder of the string after the last match.
          dest.append(srcLen - last, s + last);
      free(s);
      s = strdup(dest.text());
      return s;
  }

} //namespace Util

