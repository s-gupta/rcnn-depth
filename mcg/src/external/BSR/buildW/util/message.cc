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

#include <stdio.h>
#include <iostream>
#include <string.h>
#include "message.hh"
#include "util.hh"
#include "timer.hh"
#include "string.hh"
#include "configure.hh"

namespace Util {

  int Message::debuglevel = 1;
  int Message::frame = 0;
  char* Message::message[Message::MAX_DEPTH];
  int Message::numSteps[Message::MAX_DEPTH];
  int Message::currentStep[Message::MAX_DEPTH];
  int Message::stepSize[Message::MAX_DEPTH];
  Timer Message::timers[Message::MAX_DEPTH];

  void Message::registerConfig()
  {
    static bool called = false;
    if (called)
    {
      return;    
    }
    Configure::registerInt("Message::debuglevel",0,"Debug message level 0-5, bigger is more verbose");
    called = true;
  }

  void Message::init()
  {
    static bool called = false;
    if (called)
    {
      return;    
    }
    debuglevel = Configure::getInt("Message::debuglevel");
    called = true;
  }

  void Message::error(const char* mess)
  {
    std::cerr << "!!ERROR!! : " << mess << std::endl;
  }

  void Message::error(const String mess)
  {
    error(mess.text());
  }

  void Message::debug(const char* mess, const int level)
  {
    if (level <= debuglevel)
    {
      indent(frame);
      std::cerr << mess << std::endl;
    }
  }

  void Message::debug(const String mess, const int level)
  {
    debug(mess.text(),level);
  }

  void Message::startBlock(const char* mess, const int level)
  {
    if (level <= debuglevel)
    {
      frame++;
      assert(frame < MAX_DEPTH);

      int len = strlen(mess);
      message[frame] = new char[len+1];
      strcpy(message[frame],mess);
   
      numSteps[frame] = 0;
      stepSize[frame] = 0;
      currentStep[frame] = 0;

      indent(frame-1);
      std::cerr << "entering " << message[frame] << std::endl;

      timers[frame].reset();
      timers[frame].start();
    }
  }

  void Message::startBlock(const int totalSteps, const char* mess, const int level)
  {
    if (level <= debuglevel)
    {
      frame++;
      assert(frame < MAX_DEPTH);

      int len = strlen(mess);
      message[frame] = new char[len+1];
      strcpy(message[frame],mess);
   
      numSteps[frame] = totalSteps;
      stepSize[frame] = (int) Util::max((float)totalSteps / NUM_DOTS,2.0f);
      currentStep[frame] = 0;

      indent(frame-1);
      std::cerr << message[frame] << "  [";

      timers[frame].reset();
      timers[frame].start();
    }
  }

  void Message::stepBlock(const int level)
  {
    if ((level <= debuglevel) && (stepSize[frame] != 0))
    {
      if ((currentStep[frame] % stepSize[frame]) == 1)
      {
        std::cerr << ".";    
      }
      currentStep[frame]++;
    }
  }

  void Message::endBlock(const int level)
  {
    if (level <= debuglevel)
    {

      timers[frame].stop();
      if (numSteps[frame] > 0)
      {
        std::cerr << "] done (" << Timer::formatTime(timers[frame].cpu()) << " elapsed) " << std::endl;
      }
      else
      {
        assert(frame > 0);
        indent(frame-1);
        std::cerr << "leaving " << message[frame] << " (" 
                  << Timer::formatTime(timers[frame].cpu()) << " elapsed) " << std::endl;
      }

      delete[] message[frame];
      frame--;
      assert(frame >= 0);
    }
  }

  void Message::indent(const int inset)
  {
    for (int i = 0; i < INDENT_SIZE*inset; i++)
    {
      std::cerr << " ";    
    }
  }

} // namespace Util

