
%{

// Copyright (C) 2002 David R. Martin <dmartin@eecs.berkeley.edu>
// Copyright (C) 2002 Charless C. Fowlkes <fowlkes@eecs.berkeley.edu>
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

// Simple lexer that (1) strips out comments, (2) parses strings
// including escape caracters, (3) produces a stream of the remaining
// whiespace-separated tokens.  The comment character is '#'.

#include "string.hh"

// This is used in the mini-parser for strings (state STRING).
static Util::String str;

// Callbacks.
extern void __cfl_SyntaxError (const char* msg);
extern void __cfl_FoundToken (const char* s);
extern void __cfl_FoundString (const char* s);
extern void __cfl_AllDone ();

// Track line and column numbers.
// LEX's yylineno doesn't work at all.
int __cfl_lineNo = 1;
int __cfl_columnNo = 0;
static void eat () {
    for (int i = 0; i < yyleng; i++) {
	if (yytext[i] == '\n') {
	    __cfl_lineNo++;
	    __cfl_columnNo = 0;
	} else {
	    __cfl_columnNo++;
	}
    }
}

%}

WS		[ \t\r\n]+
TOKEN		[^ \t\r\n#\"]+

%x STRING 
%option noyywrap 

%%

#.*$			eat(); // discard comments
{TOKEN} 		__cfl_FoundToken (__cfl_text); eat();

\"			str.clear(); BEGIN(STRING); eat();
<STRING>\"		__cfl_FoundString (str.text()); BEGIN(INITIAL); eat();
<STRING>\\0		str.append ('\0'); eat();
<STRING>\\a		str.append ('\a'); eat();
<STRING>\\b		str.append ('\b'); eat();
<STRING>\\t		str.append ('\t'); eat();
<STRING>\\n		str.append ('\n'); eat();
<STRING>\\v		str.append ('\v'); eat();
<STRING>\\f		str.append ('\f'); eat();
<STRING>\\r		str.append ('\r'); eat();
<STRING>\\.		str.append (__cfl_text[1]); eat();
<STRING>\r?\n		|
<STRING><<EOF>>		__cfl_SyntaxError ("Unterminated string"); eat();
<STRING>.		str.append(__cfl_text[0]); eat();

{WS}			eat(); // discard whitespace
.			__cfl_SyntaxError ("Unknown token"); eat();

<<EOF>>			__cfl_AllDone(); return 0;

%%


 
