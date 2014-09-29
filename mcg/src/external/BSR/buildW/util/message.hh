// Copyright (C) 2002 Charless C. Fowlkes <fowlkes@eecs.berkeley.edu>
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

#ifndef MESSAGE_HH
#define MESSAGE_HH

#include "timer.hh"
#include "string.hh"

namespace Util {
  class Message 
  {
    public:
       static void registerConfig();
       static void init();

       static void error(const String mess);
       static void error(const char* mess);
       static void debug(const String mess, const int level = 0);
       static void debug(const char* mess, const int level = 0);
       static void startBlock(const char* message, const int level = 0);
       static void startBlock(const int totalSteps, const char* message, const int level = 0);
       static void stepBlock(const int level = 0);
       static void endBlock(const int level = 0);

    private:
       static void indent(const int inset);

       static const int NUM_DOTS=20;
       static const int MAX_DEPTH=500;
       static const int INDENT_SIZE=4;

       static int debuglevel;
       static int frame;                        //current frame depth
       static char* message[MAX_DEPTH];         //the messages associated with each frame
       static int numSteps[MAX_DEPTH];          //total number of steps for given frame
       static int currentStep[MAX_DEPTH];       //current step
       static int stepSize[MAX_DEPTH];          //how many steps between updates
       static Timer timers[MAX_DEPTH];          //timers for each block
  };
} // namespace Util

#endif // __Message_hh__

