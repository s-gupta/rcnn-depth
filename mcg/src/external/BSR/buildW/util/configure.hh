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

#ifndef CONFIGURE_HH
#define CONFIGURE_HH

// TODO: print null strings in comments

#include <stdio.h>
#include "types.hh"

  // This module maintains a global map of the form module::option -->
  // (value), where the values are typed.  The map is created when modules
  // register options.  The map can be modified by command-line options
  // and configuration files.  Configuration files are read in the order
  // they appear on the command-line.  Command-line paramters override
  // configuration file parameters.

  // Command-line options are of the form "module::option=value" or
  // "-module::option value".  The "--" argument terminates command-line
  // processing.

  // The option "Config=<file>" is reserved, and causes parameter values
  // to be read from <file>.  In the configuration file, '#' is the
  // comment character.  The configuration file consists of a sequence
  // of whitespace-separated "module::option value" pairs.  The
  // "module::option=value" syntax is not supported in the configuration
  // file.

  // In all method calls, module may be "" or NULL.  Either signifies
  // an empty module name.
  class Configure {
    public:
      //static const char *const configKey = "config";  //moved to Configure.cc

      // Process options on command-line.
      // The index at which processing stopped is returned.
      static Util::uint init(int argc, const char **argv);

      // Initialize an individual option.
      static void init(const char *key, const char *value);

      // Print a usage message showing default values and option descriptions.
      static void usage();

      // Print a table of parameters and values.
      static void show();

      // Register options along with default values.
      // For ENUM, values is a NULL-terminated list of possible values
      // for the option.
      static void registerBool(const char *key, bool defaultValue,
                               const char *description);
      static void registerInt(const char *key, Util::int64 defaultValue,
                              const char *description, bool printHex =
                              false);
      static void registerFloat(const char *key, double defaultValue,
                                const char *description);
      static void registerString(const char *key, const char *defaultValue,
                                 const char *description);
      static void registerEnum(const char *key, const char *const *values,
          Util::uint defaultValue, const char *description);

      // Retrieve option values.  The string returned by getString()
      // is created with strdup().  Caller is responsible for freeing it.
      static bool getBool(const char *key);
      static Util::int64 getInt(const char *key);
      static double getFloat(const char *key);
      static char *getString(const char *key);
      static Util::uint getEnum(const char *key);

  };

#endif                          // __Configure_h__
